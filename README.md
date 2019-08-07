# CGCP-Perl
“Causual genotype combination patterns (CGCPs)” method is an exhaustive process developed by ourselves with a script written in Perl language that allows you to identify causual genotype combination patterns in complicated disease.
Two core algorithm are executed by CGCP. One is that CGCP program can screen the variants combinations existed only in patients while not in control. Another is that after eliminating the patterns existing in the control cohort, the remaining patterns would be examined by frequency according to the prevalence rate as a boundary condition to limit the scale of calculation. Combined with the prevalence of a particular complex disease, CGCP can be used to simulate the frequency of genotypes in the population by magnificating times for controls from experiment. The sum of frequency of disease-specific genotype combination in the simulated population that clearly do not meet the prevalence of disease are removed by CGCP. Once the procedure of CGCP is complete, specific genotype combinations are selected.
## How it works
### CGCPs steps
<br> In order to facilitate the operation and reduce the amount of computation, before started CGCP，you can select loci that you intrested in. 
* Preperation of raw data of genome
>> PLINK Binary files (BED/BIM/FAM) from genome sequencing raw data
* Extract SNPs (single nucleotide polymorphisms) in genome regions of interest
>> Executed by PLINK 
* Select independent SNPs
>> Considering Linkage disequilibrium, pruning redundancy SNPs to select independent SNPs is necessary, executed by PLINK
>> <br> Details of commend on these repository "pruning_SNPs" 
* Association study
>> Using Perl language to excute command of PLINK-association
>> <br> Details of commend on these repository "plink-bed.pl"
>> <br> After quality control, a batch of SNPs were screened out
* Annotate
>> Made every SNPs annotate into different genes.注释方法及命令需要修改
* Search CGCPs
>> This is the core step. 
>> <br> In the search process, all the genotype combination patterns would be searched in each
>> <br> individual among the cases and controls. Only candidate CGCPs were picked up, while the
>> <br> others, if they occurred over once in the controls, would be eliminated from the process.
>> <br> Details of commend on these repository "CGCPs.pl.txt"
* Merging the data of CGCPs 
>> Patients who had same CGCPs will be merged and statisticed.
>> <br> Details of commend on these repository "statistic.pl.txt"

## Runing test
Ruuning test need to follow above steps. Test documents were performed by Perl 5.01 and PLINK 1.07 and we upload all files on these repository

## Download/Installation
* Perl, a family of two high-level, general-purpose, interpreted, dynamic programming languages. Need to install above version 5.01.
<br> Copyright (C) 1989 Free Software Foundation, Inc. 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
* <br> PLINK, an open-source C/C++ WGAS tool set. You need to install above 1.07.
<br> Copyright (C) 2006 Shaun Purcell, GNU General Public License, v2  http://pngu.mgh.harvard.edu/purcell/plink/ 

## Notice
根据不同的数据类型，注释文件有所不同

## License
This project is licensed under the MIT License - see the LICENSE.md file for details

## Question or Feedback
Email us: zhengxiaodong.1234@163.com，we're listening!
