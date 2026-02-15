#!/usr/bin/env python3
"""
CGCP_pro.py: Causative Genotype Combination Pattern Analysis
Author: XiaoDong Zheng
License: MIT
Description: A high-performance tool to identify multi-locus genotype combinations 
             potentially causative for complex diseases. 
             Optimized for HWE file parsing and automated output naming.
"""

import pandas as pd
import numpy as np
from itertools import combinations
import time
from multiprocessing import Pool, cpu_count
import sys

# Ensure tqdm is installed for progress bar support
try:
    from tqdm import tqdm
except ImportError:
    class tqdm:
        def __init__(self, total, **kwargs): self.n = 0
        def update(self, n): self.n += n
        def close(self): pass

def process_combination_batch(args):
    """
    Worker function to process a chunk of SNP combinations.
    """
    (combo_batch, geno_df_values, snp_list, case_mask, ctrl_mask, 
     freq_map, lamda, low_boundary, up_boundary, cutoff_depth) = args
    
    batch_results = []
    
    for combo_indices in combo_batch:
        combo_snps = [snp_list[i] for i in combo_indices]
        
        # Pre-check: skip if any SNP in the combination is missing from the frequency map
        if any(snp not in freq_map for snp in combo_snps):
            continue

        # Slicing the genotype matrix for current SNP combination
        sub_matrix = geno_df_values[:, combo_indices]
        
        # Filter out samples with missing genotypes ('0' or '00')
        mask = ~(np.any(sub_matrix == '0', axis=1) | np.any(sub_matrix == '00', axis=1))
        valid_matrix = sub_matrix[mask]
        valid_case_mask = case_mask[mask]
        valid_ctrl_mask = ctrl_mask[mask]
        
        if len(valid_matrix) == 0:
            continue

        # Vectorized string concatenation for multi-locus genotypes
        combined_genos = np.array(["".join(row) for row in valid_matrix])
        
        # Frequency counts using NumPy unique
        u_cases, c_counts = np.unique(combined_genos[valid_case_mask], return_counts=True)
        u_ctrls, ct_counts = np.unique(combined_genos[valid_ctrl_mask], return_counts=True)
        ctrl_dict = dict(zip(u_ctrls, ct_counts))
        
        for genotype_str, count in zip(u_cases, c_counts):
            # CGCP Criteria: Case frequency >= cutoff_depth AND Control frequency <= 1
            if count >= cutoff_depth and ctrl_dict.get(genotype_str, 0) <= 1:
                
                # Population frequency calculation (HapFrq)
                # Formula: (FreqCase + FreqCtrl * lambda) / (1 + lambda)
                pop_freq = 1.0
                for i, snp in enumerate(combo_snps):
                    pair = genotype_str[i*2 : i*2+2]
                    # Safe dictionary access using nested .get() to prevent KeyError
                    f_case = freq_map[snp].get('AFF', {}).get(pair, 0)
                    f_ctrl = freq_map[snp].get('UNAFF', {}).get(pair, 0)
                    pop_freq *= (f_case + f_ctrl * lamda) / (1 + lamda)
                
                # Boundary filtering based on disease prevalence
                if low_boundary < pop_freq < up_boundary:
                    res_line = (f"SNPs: {' '.join(combo_snps)}\tGeno:{genotype_str}\t"
                                f"CaseNum:{count}\tCtrlNum:{ctrl_dict.get(genotype_str,0)}\t"
                                f"PopFreq:{pop_freq:.6f}")
                    batch_results.append(res_line)
    
    return batch_results

def main():
    if len(sys.argv) < 5:
        print("\nUsage: python3 CGCP_pro.py <prefix> <prevalence> <min_depth> <comb_size> [cpus]")
        print("Example: python3 CGCP_pro.py my_data 0.0047 20 3 36\n")
        sys.exit(1)

    prefix = sys.argv[1]
    prevalence = float(sys.argv[2])
    cutoff_depth = int(sys.argv[3])
    comb_size = int(sys.argv[4])
    num_processors = int(sys.argv[5]) if len(sys.argv) > 5 else max(1, cpu_count() - 2)

    print(f"\n{'='*60}")
    print(f" CGCP Analysis Pro - Performance Engine")
    print(f"{'='*60}")
    print(f"Start Time:  {time.ctime()}")
    print(f"Data Prefix: {prefix}")
    print(f"Prevalence:  {prevalence}")
    print(f"Min Depth:   {cutoff_depth}")
    print(f"Combo Size:  {comb_size}")
    print(f"Processors:  {num_processors}")

    # 1. Load MAP file (SNP Identifiers)
    map_df = pd.read_csv(f"{prefix}.map", sep=r'\s+', header=None)
    snp_list = map_df[1].tolist()
    
    # 2. Load HWE file and Parse Frequencies
    lamda = (1 - prevalence) / prevalence
    freq_map = {}
    print("Parsing HWE frequencies (Filtering 'ALL' rows to avoid conflicts)...")
    with open(f"{prefix}.hwe", 'r') as f:
        for line in f:
            p = line.strip().split()
            if not p or p[0] == "CHR": continue
            
            # Based on standard PLINK HWE format: SNP(index 1), TEST(index 2)
            snp, status, a1, a2, geno_raw = p[1], p[2], p[3], p[4], p[5]
            
            # CRITICAL: Only process AFF and UNAFF, skip ALL to prevent dict overwrite
            if status not in ['AFF', 'UNAFF']:
                continue
                
            if snp not in freq_map:
                freq_map[snp] = {'AFF': {}, 'UNAFF': {}}
            
            try:
                nums = [int(x) for x in geno_raw.split('/')]
                tot = sum(nums)
                if tot > 0:
                    freq_map[snp][status][f"{a1}{a1}"] = nums[0]/tot
                    freq_map[snp][status][f"{a1}{a2}"] = nums[1]/tot
                    freq_map[snp][status][f"{a2}{a2}"] = nums[2]/tot
            except:
                continue

    # 3. Load PED file and Vectorize
    print("Loading Genotypes (PED)...")
    ped_data = pd.read_csv(f"{prefix}.ped", sep=r'\s+', header=None)
    case_mask = (ped_data[5].values == 2)
    ctrl_mask = (ped_data[5].values == 1)
    
    # Concatenate allele columns into genotype matrix
    geno_raw_matrix = ped_data.iloc[:, 6:].values
    combined_geno_list = []
    for i in range(0, geno_raw_matrix.shape[1], 2):
        combined_geno_list.append(np.char.add(geno_raw_matrix[:, i].astype(str), 
                                              geno_raw_matrix[:, i+1].astype(str)))
    geno_values = np.array(combined_geno_list).T
    
    # 4. Prepare Combinations and Task Chunks
    all_combos = list(combinations(range(len(snp_list)), comb_size))
    total_combos = len(all_combos)
    num_chunks = num_processors * 4 
    chunk_size = (total_combos // num_chunks) + 1
    chunks = [all_combos[i:i + chunk_size] for i in range(0, total_combos, chunk_size)]
    
    # 5. Multiprocessing Execution
    print(f"Testing {total_combos} combinations...")
    fixed_params = (geno_values, snp_list, case_mask, ctrl_mask, freq_map, 
                    lamda, 0.01 * prevalence, prevalence, cutoff_depth)
    task_args = [(c,) + fixed_params for c in chunks]

    final_results = []
    with Pool(processes=num_processors) as pool:
        pbar = tqdm(total=total_combos, desc="Analysis Progress", unit="comb")
        for batch_res in pool.imap_unordered(process_combination_batch, task_args):
            final_results.append(batch_res)
            # Safe progress update
            pbar.update(len(chunks[0]) if pbar.n + len(chunks[0]) <= total_combos else total_combos - pbar.n)
        pbar.close()

    # 6. Automatic Output Naming
    output_file = f"result_depth{cutoff_depth}_size{comb_size}.txt"
    with open(output_file, "w") as f_out:
        total_found = 0
        for batch in final_results:
            for line in batch:
                f_out.write(line + "\n")
                total_found += 1
    
    print(f"\n{'='*60}")
    print(f" Analysis Completed: {time.ctime()}")
    print(f" Total CGCPs Found:  {total_found}")
    print(f" Results Saved to:   {output_file}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
