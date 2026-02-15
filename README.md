## CGCPs: Causative Genotype Combination Patterns Analysis
CGCPs is a bioinformatic tool designed with an exhaustive process. It identifies highly specific multi-locus genotype combination patterns that are potentially causative for complex diseases.

## CGCPs steps
# Data preperation
In order to facilitate the operation and reduce the amount of computation, before started CGCP，you can select loci that you intrested in. 
Prepared genotype/sequencing data to plain text PLINK format(.ped/.map/.hwe). Test files are showed in "Test documents" folder named "select-interested-loci(.ped/.map/.hwe)".

# Clone the repository
git clone https://github.com/XiaoDongZheng1234/CGCPs.git

# Install dependencies
pip install pandas numpy tqdm

# Run the analysis with 36 cores
python3 CGCP_pro.py <data_prefix> <prevalence> <min_depth> <comb_size>

# Example:
python3 CGCP_pro.py select-interested-loci 0.0047 20 3
`select-intrested-loci`: The script automatically searches for files with the suffixes `.map`, `.ped`, and `.hwe`.
`0.0047`: This is the prevalence of psoriasis in the Han Chinese population. The script uses this value to calculate the frequency across the entire population.
`20`: Indicates that a specific genotype combination occurs in at least 20 patients (Depth).
`3`: Indicates that pairwise analysis of 3 SNPs is performed.

## Key Features
# High Performance: 
Re-engineered from Perl to Multiprocessing Python, achieving up to 50x speedup on multi-core servers.
# Big Data Ready: 
Optimized for UK Biobank scale datasets using NumPy vectorized operations.
# Real-time Monitoring: 
Integrated tqdm progress bars for tracking millions of combinations.
# Strict Filtering: 
Built-in prevalence-based frequency weighting and strict case-control filtering logic.
