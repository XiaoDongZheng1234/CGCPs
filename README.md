# CGCP-Perl
"<font color=Red>Causual genotype combination patterns (CGCPs)</font>" method is an exhaustive process developed by ourselves with a script written in Perl language that allows you to identify causual genotype combination patterns in complicated disease.<br><br>Two core algorithm are executed by CGCP. One is that CGCP program can screen the variants combinations existed only in patients while not in control. Another is that after eliminating the patterns existing in the control cohort, the remaining patterns would be examined by frequency according to the prevalence rate as a boundary condition to limit the scale of calculation. Combined with the prevalence of a particular complex disease, CGCP can be used to simulate the frequency of genotypes in the population by magnificating times for controls from experiment. The sum of frequency of disease-specific genotype combination in the simulated population that clearly do not meet the prevalence of disease are removed by CGCP.<br> <br>Once the procedure of CGCP is complete, specific genotype combinations are selected.

## How it works
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### CGCPs steps
* **Data preperation**
<br> In order to facilitate the operation and reduce the amount of computation, before started CGCP，you can select loci that you intrested in. 
>> Prepared genotype/sequencing data to plain text PLINK format(.ped/.map/.hwe). Test files are showed in "Test documents" folder named "select-interested-loci(.ped/.map/.hwe)".
* **Running CGCP**
<br> The ruuning command and its usage warinings are as follow. You can also read the detail in this repository ("CGCP.pl.txt"). After running CGCP with the test file, the generated file is named "CGCP-test.txt".Test documents were performed by Perl 5.01.
```perl
#!/usr/bin/perl
use strict;
use 5.010;
#use warnings;
# run command: nohup perl CGCP.pl.txt select-intrested-loci 0.0047 10 3 > myfile.txt &
#"select-interested-loci": plink ".ped/.map/.hwe" file profix name; "0.0047": prevalence rate of specific ; "10": cutoff depth; "3": number of selected snps for combination analysis; $low_boundary=0.01*$prevalence; $up_boundary=$prevalence;

my $nowtime=localtime;
say "$nowtime";
#reading data
#Read data of variants
open IA, "$ARGV[0].map";
my %hash;
my @snp = ();
my @sample=();
while(<IA>){
chomp;
my $line1=$_;
my @line1=split/\s+/,$line1;
push (@snp,$line1[1]);
}
my @snpnum=(0..$#snp);
my %phenotype;
my %genotype;
close IA;

#Read data of genotype and build hash
open IB, "$ARGV[0].ped";
while(<IB>){
chomp;
my $line2=$_;
my @line2=split/\s+/,$line2;
push (@sample,$line2[0]);
$phenotype{$line2[0]}="$line2[5]";
	foreach(@snpnum){
	my $num=$_;
	my $aa=6+$num*2;
	my $bb=7+$num*2;
	$genotype{$line2[0].$snp[$_]}="$line2[$aa]$line2[$bb]";
	}
}
close IB;

#Read "/.hwe" file and simulate the genotype frequency of each variant in the population；
open IO, "$ARGV[0].hwe";
my %freqcase;
my %freqcontrol;
while(<IO>){
chomp;
my $line=$_;
my @line=split/\s+/,$line;
	if($line[3] =~/^AFF/){
	my @numaff=split/\//,$line[6];
	my $numaff=$numaff[0]+$numaff[1]+$numaff[2];
		if($numaff > 0){
		$freqcase{$line[2].$line[4].$line[4]}=$numaff[0]/$numaff;
		$freqcase{$line[2].$line[4].$line[5]}=$numaff[1]/$numaff;
		$freqcase{$line[2].$line[5].$line[5]}=$numaff[2]/$numaff;
		}
	}elsif($line[3] =~/UNAFF/){
	my @numunaff=split/\//,$line[6];
	my $numunaff=$numunaff[0]+$numunaff[1]+$numunaff[2];
		if($numunaff > 0){
		$freqcontrol{$line[2].$line[4].$line[4]}=$numunaff[0]/$numunaff;
		$freqcontrol{$line[2].$line[4].$line[5]}=$numunaff[1]/$numunaff;
		$freqcontrol{$line[2].$line[5].$line[5]}=$numunaff[2]/$numunaff;
		}
	}else{
	next;
	}
}



#Define the prevalence (P),Lamda, and the number of variants that make up genotype combinations
my $prevalence=$ARGV[1];
my $low_boundary=0.01*$prevalence;
my $up_boundary=$prevalence;
my $lamda=(1-$prevalence)/$prevalence;
my $deep=$ARGV[2];
my $size=$ARGV[3];

###########Master Program################
foreach(@snpnum){
&check($_,$size);
}

###########Subprogram################
sub check{
my $value = shift;
my $max = (split (/\t/,$value))[-1];
my $size = shift;
my @tmpList;
        foreach(@snpnum){
                if($max < $_){
                        push(@tmpList,$_);
                }
	
        }

        foreach(@tmpList){
        my $tmpValue = "$value"."\t$_";
	my @tmpValue = split /\t/,$tmpValue;

        
                if($#tmpValue + 1 == $size){
		my @snplistnum=@tmpValue;
		my $snpselect;
		my %haplotype;
		my %haplotypenum;
		my %haplotypenumc;
		my %haplotypesam;
		my %haplotypesamc;
			foreach(@snplistnum){
			$snpselect.="$snp[$_]\t";
			}

#Built the haplotype of case；
			foreach my $samplecase(@sample){
				if ($phenotype{$samplecase} eq 2){
				my $genotypecase=undef;
					foreach(@snplistnum){
					$genotypecase.=$genotype{$samplecase.$snp[$_]};
					}
					if($genotypecase !~0){
					$haplotype{$genotypecase}="1";
					$haplotypenum{$genotypecase}+=1;
					$haplotypesam{$genotypecase}.="$samplecase".":";########################
					}
				}
			}

#Built the haplotype of control；
			foreach my $samplecontrol(@sample){
				if ($phenotype{$samplecontrol} eq 1){
				my $genotypecontrol=undef;
					foreach(@snplistnum){
					$genotypecontrol.=$genotype{$samplecontrol.$snp[$_]};
					}
					if($genotypecontrol !~0){
					$haplotypenumc{$genotypecontrol}+=1;
					$haplotypesamc{$genotypecontrol}.="$samplecontrol".":";######################
					}	
				}
			}

#Remove the haplotype present both in case and control
			foreach my $key(sort keys %haplotype){
				delete $haplotype{$key};
			}

#review if the frequency of haplotype between 0.01P and P；
		my @snpselect=split/\t/,$snpselect;
		my %hap=();
			foreach my $hap (sort keys %haplotype){
			my @hap=split//,$hap;
			my $hapfrq=1;
				foreach my $nn(0..$#snpselect){
				my $aa=2*$nn;
				my $bb=2*$nn+1;
				my $hapfrqadd=($freqcase{$snpselect[$nn].$hap[$aa].$hap[$bb]}+$freqcontrol{$snpselect[$nn].$hap[$aa].$hap[$bb]}*$lamda)/(1+$lamda);
				$hapfrq*=$hapfrqadd;
				}
			$hap{$hap}=$hapfrq;
				unless($hapfrq>$low_boundary and $hapfrq<$up_boundary){
				delete $haplotype{$hap};
				}
				if($haplotypenum{$hap} < $deep){
				delete $haplotype{$hap};
				}
			}

#Output variants list and haplotype that meet the above conditions
		my @haplotype =	keys %haplotype;	
			if(%haplotype){ 
			print "\n$snpselect\t";
				foreach(@haplotype){
					if(exists $haplotypenumc{$_}){
					print "$_\t$haplotypenum{$_}\t$haplotypesam{$_}\t$haplotypenumc{$_}\t$haplotypesamc{$_}\t$hap{$_}\t";##################
					}else{
					print "$_\t$haplotypenum{$_}\t$haplotypesam{$_}\t0\tundef\t$hap{$_}\t";###################
					}	
				}
			}
                next;
                }


        &check($tmpValue,$size);
        }
}
my $outtime=localtime;
say "$outtime";//perl
```
* **Variants annotation and enrichment for the genotype combinations**
>> Merged these genotype combinations selected by CGCP into gene or loci combinations if it has the same gene symbol and physical location.


## Download/Installation
* Perl, a family of two high-level, general-purpose, interpreted, dynamic programming languages. Need to install above version 5.01.
<br> Copyright (C) 1989 Free Software Foundation, Inc. 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


## License
This project is licensed under the MIT License - see the LICENSE.md file for details

## Question or Feedback
Email us: zhengxiaodong.1234@163.com，we're listening!
