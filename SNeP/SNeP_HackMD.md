---
title: 'SNeP: trends in recent effective population size trajectories using genome-wide SNP data'
disqus: hackmd
---

SNeP: trends in recent effective population size trajectories using genome-wide SNP data
===
Read more about SNeP from Mario Barbato, here:
https://www.frontiersin.org/articles/10.3389/fgene.2015.00109/full

## Objective
We're interested in looking at effective population size trends in "contemporary" demographic history (most recent ~1000 generations) that would be complementary to PSMC and likely more of interest to our state agency folks than PSMC. SNeP and IBDNe are two methods to do this, but SNeP has proven easier to use because IBDNe requires a full pedigree and linkage map which we don't have at this time. 

## Background
Effective population size (Ne) is a key population genetic parameter that describes the amount of genetic drift in a population. Methods to estimate Ne from linkage disequilibrium (LD) were developed ~40 years ago but depend on the availability of large amounts of genetic marker data that only the most recent advances in DNA technology have made available. One solution to overcome the limitation of an incomplete pedigree is to estimate the recent trend in Ne using genomic data. Here we introduce SNeP, a multithreaded tool to perform the estimate of Ne using LD using the standard PLINK input file format (.ped and.map files) or by using LD values calculated using other software. Through SNeP the user can apply several corrections to take account of sample size, mutation, phasing, and recombination rate. 

## Table of Contents

[TOC] 

Installation for first time users
---
1. Load Dependencies
You can download the SNeP binaries here. I have found that this is not straightforward at all: SNeP1.1 is for Linux
https://sourceforge.net/projects/snepnetrends/files/
Or use GitHub
git clone https://git.code.sf.net/p/snepnetrends/code snepnetrends-code #this repository is empty!
```
module load anaconda/
module load gcc/4.9.2
```
Some of our packages are installed into a conda environment at: /home/tl50a/.conda/envs/NGSStep7:

We had a **major** dependency issue which was resolved by Chris Hull (MGHPCC)
./SNeP1.1: /lib64/libc.so.6: version `GLIBC_2.14' not found (required by ./SNeP1.1)
./SNeP1.1: /lib64/libc.so.6: version `GLIBC_2.17' not found (required by ./SNeP1.1)

Running SNeP: 
---
We have a work around that needs to be followed very carefully here: 
1. load singularity
module load singularity/singularity-3.5.3
2. Navigate to ./SNeP1.1 and Run
chmod u+x ./SNeP1.1
3. Instead of ./SNep1.1, use:
```
singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif ./SNeP1.1
```
This method used a singularity container to run the SNeP1.1 binary, despite the new linux system not supporting the glibc libraries that were required. You can add options or filenames on to the end of the command, as normal e.g. 
```
singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif ./SNeP1.1 -help
```
SNep works from the command line using flags as options. All of the options are discussed in the manual, and we include those that seemed relevant in the command below. SNeP's default parameters should fit MOST requirements. 

If there are no problems with the input files, SNeP will quickly (<1h) perform the analysis and create a filepath/filename.NeAll output file and a filepath/filenameSNeP.log logging file.

The ***input files*** for SNeP are map/ped files based on genome-wide SNP data. SNeP does not use whole genomes e.g. in bam format, just a VCF converted to ped/map format using vcftools or plink.) 

We've already done the format conversion from VCF to ped/map, so we will proceed with those files located here: 
/project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.ped 6SV_unfiltered_SNPs.map

These ped/map files include genome-wide SNPs that are minimally filtered (e.g. without LD-pruning). However, all of the SNPs included have been filtered for depth of coverage and other quality parameters which you can read about on my WGS Pipeline HackMD. OK carrying on:

4. I recommend using an interactive session: 
```
bsub -q interactive -R rusage[mem=24000] -n 1 -W 300 -Is bash
```
5. Running all lynx together n=61)
```
bsub -q long -R rusage[mem=2000] -n 8 -W 24:00 "singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif ./SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/ -threads 8" 
```
### Results
This worked well and produced an output with the suffix .NeAll. We will pull this into R and run some visualizations, then compare by population. But first, we'll do a little sensitivity analysis to better understand the performance of SNeP and make a few decisions regarding parameterization.

## Sensitivity Testing
Now that we have gotten SNeP to run properly, we need to ask a couple questions. 

1. Is SNeP sensitive to the number of chromosomes?
2. Is SNeP sensitive to the number of samples/population?
3. Is SNeP sensitive to parameterization?

We have run a full suite or chromomsomes (all, some, one). We've found that SNeP is fairly sensitive to the number of SNPs and number of chromosomes you feed into it, but once you reach a certain threshold (in our case N SNPs distributed across n chromosomes), the Ne trend produced via the .NeAll output begins to look the same. 

**This is the full list of all 66 chromosome names:**
NC_044303.1,NC_044304.1,NC_044305.1,NC_044306.1,NC_044307.1,NC_044308.1,NC_044309.1,NC_044310.1,NC_044311.1,NC_044312.1,NC_044313.1,NC_044314.1,NC_044315.1,NC_044316.1,NC_044317.1,NC_044318.1,NC_044319.1,NC_044320.1,NC_044321.1,NW_022059692.1,**NW_022059693.1**,NW_022059694.1,NW_022059695.1,NW_022059696.1,NW_022059697.1,NW_022059698.1,NW_022059699.1,NW_022059700.1,NW_022059701.1,NW_022059702.1,NW_022059703.1,NW_022059704.1,NW_022059705.1,NW_022059706.1,NW_022059707.1,NW_022059708.1,NW_022059709.1,NW_022059710.1,NW_022059711.1,NW_022059712.1,NW_022059713.1,NW_022059714.1,NW_022059715.1,NW_022059716.1,NW_022059717.1,NW_022059718.1,NW_022059719.1,NW_022059720.1,NW_022059721.1,NW_022059722.1,NW_022059723.1,NW_022059724.1,NW_022059725.1,NW_022059726.1,NW_022059727.1,NW_022059728.1,NW_022059729.1,NW_022059730.1,NW_022059731.1,NW_022059732.1,NW_022059733.1,NW_022059734.1,NW_022059735.1,NW_022059736.1,NW_022059737.1,NW_022059738.1, NW_022059694.1

**This is a list of all individuals**: 
L155.variant4.variant2  LIC11.variant4.variant2 LIC20.variant4.variant2 LIC23.variant4.variant2 LIC24.variant4.variant2 LIC27B.variant5.variant2        LIC28.variant5.variant2 LIC31.variant5.variant2 LIC32.variant5.variant2 LIC36.variant5.variant2 LIC46.variant5.variant2 LIC47.variant5.variant2 LIC48.variant6.variant2 LIC54.variant6.variant2 LIC57.variant6.variant2 LIC60.variant6.variant2 LIC8.variant4.variant2  LIC9.variant4.variant2  LIT2.variant6.variant2  LIT5.variant6.variant2  LRK10.variant6.variant2 LRK11.variant7.variant2 LRK12.variant7.variant2 LRK13.variant7.variant2 LRK17.variant7.variant2 LRK22.variant7.variant2 a109.variant.variant.variant2   a182.variant.variant.variant2   a202.variant.variant.variant2   a33.variant.variant.variant2    a475.variant.variant.variant2   a494.variant.variant.variant2   a507.variant.variant.variant2   a697.variant    a772.variant    a794.variant        a857.variant          b90.variant2.variant.variant2   a772.variant a803.variant    a818.variant b114.variant    b124.variant2.variant.variant2 b188.variant2.variant.variant2  b23.variant2.variant.variant2   b276.variant2.variant.variant2  b554.variant2.variant.variant2 c165.variant2.variant.variant2  c323.variant2.variant2  c548.variant2.variant2  cb15.variant2.variant2  cb42.variant2.variant2  cb7.variant2.variant2         l09_003.variant3.variant2              l09_015.variant3.variant2

### What does an individual chromosome tell us? 
```
singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/allsamples_NW_022059698.1 -chr NW_022059698.1 -maf 0 -itemsTH 1 -mindist 1 -samplesize 54
```
SNeP hit some issues at chr21 (NW_022059693.1) but we were able to get some advice from the author (Mario Barbato) that resolved it. 

### Results
Actually you can glean quite a lot from a single chromosome, but I found that the Ne trend line different quite a lot from chromosome to chromosome, and didn't converge for some, where SNP density wasn't high enough for SNeP. 

Chrom: 

Chrom: 

SNeP is geared toward genome-wide SNPs, after all, so we will input SNP data from all 66 chromosomes. Sadly, CHR Y is not included in this analysis because CHR Y wasn't assembled until pretty recently. 

### Sensitivity to parameterization
We tested sensitivity to five flags of interest, run on just 2 chromosomes NC_044303.1,NC_044304.1 for expediency. 
1. Default
```
bsub -q long -R rusage[mem=2000] -n 8 -W 24:00 -R span\[hosts=1\] "singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/sensitivity_analysis/7_sensitivity_samplesize_itemsTH_mindist_maxsnp -chr NC_044303.1 -samplesize 61 -itemsTH 1 -mindist 1 -maxsnp 300000 -threads 8" 
```
2. With sample-size
sample size is 61 for all 
3. With MAF rate
* MAF rate (default 0.5)
trying MAF 0
4. Adjusted n SNPs/bin
try 1
* itemsTH items/bin (binning rate)
5. Adjusted mindist
try 1
* mindist distance between SNPs (default is 5000, but we will set mindist to >LD because some of our chromosomes are quite small. 5000 seems excessive)
6. Adjusted maxSNP
* maxSNP (default 100k/CHR)

Combination: 
-samplesize 61

-samplesize 61 -maf 0 -itemsTH 1 -mindist 1
### Results
* MAF-rate
* itemsTH
* mindist
* sample size

Population-level SNeP trends
---
Now that we've really nailed down our parameterization, we need to subset our VCF into populations, convert those VCFs to ped/map format, run SNeP for each population, and plot the different trend lines together in R. 

### 1. Subsetting VCF into populations 
Each population-level VCF and associated ped/map file will be located in the Rgenomics/mLynCan4_v1.p folder along with the "mother" 6SV VCF

We've resolved from much trial and error that the best way of subsetting a VCF is using vcftools

**Our output files will include:** 
6SV_southSLR.vcf.gz .map/.ped
6SV_northSLR.vcf.gz .map/.ped
6SV_NFLD.vcf.gz .map/.ped
6SV_bobcat.vcf.gz .map/.ped

### NFLD --done
```
bsub -q short -W 1:00 -R rusage[mem=8000] -n 2 "vcftools --indv a475.variant.variant.variant2 --indv a494.variant.variant.variant2 --indv b276.variant2.variant.variant2 --gzvcf /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz --recode --out /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_NFLD.vcf.gz"

#The output then needs to be gzipped with:
bsub -q short -W 1:00 -R rusage[mem=48000] -n 1 -R span\[hosts=1\] gzip 6SV_NFLD.vcf

#The gzvcf then needs to be recoded to ped/map format using PLINK
plink --vcf  6SV_NFLD.vcf.gz --recode --allow-extra-chr
```
### southSLR --done
```
bsub -q short -W 1:00 -R rusage[mem=8000] -n 2 "vcftools --remove-indv a475.variant.variant.variant2 --remove-indv a494.variant.variant.variant2 --remove-indv b276.variant2.variant.variant2 --remove-indv a202.variant.variant.variant2 --remove-indv a33.variant.variant.variant2 --remove-indv a772.variant --remove-indv a803.variant --remove-indv a818.variant --remove-indv b114.variant    --remove-indv b124.variant2.variant.variant2 --remove-indv b188.variant2.variant.variant2 --remove-indv b554.variant2.variant.variant2 --remove-indv c165.variant2.variant.variant2  --remove-indv c323.variant2.variant2 --remove-indv b13.variant2.variant2 --remove-indv f264.variant2.variant2  --remove-indv f457.variant3.variant2  --remove-indv fha_024.variant3.variant2 --remove-indv fha_042.variant3.variant2 --remove-indv fha_043.variant3.variant2 --remove-indv l09_007.variant3.variant2 --gzvcf /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz --recode --out /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_southSLR"

#The output then needs to be gzipped with:
bsub -q short -W 1:00 -R rusage[mem=48000] -n 1 -R span\[hosts=1\] gzip 6SV_southSLR.vcf

#The gzvcf then needs to be recoded to ped/map format using PLINK
plink --vcf  6SV_southSLR.vcf.gz --recode --allow-extra-chr
```
### mainelynx n=26
```
bsub -q short -W 1:00 -R rusage[mem=8000] -n 2 "vcftools --indv L155.variant4.variant2  --indv LIC11.variant4.variant2 --indv LIC20.variant4.variant2 --indv LIC23.variant4.variant2 --indv LIC24.variant4.variant2 --indv LIC27B.variant5.variant2        --indv LIC28.variant5.variant2 --indv LIC31.variant5.variant2 --indv LIC32.variant5.variant2 --indv LIC36.variant5.variant2 --indv LIC46.variant5.variant2 --indv LIC47.variant5.variant2 --indv LIC48.variant6.variant2 --indv LIC54.variant6.variant2 --indv LIC57.variant6.variant2 --indv LIC60.variant6.variant2 --indv LIC8.variant4.variant2  --indv LIC9.variant4.variant2  --indv LIT2.variant6.variant2  --indv LIT5.variant6.variant2  --indv LRK10.variant6.variant2 --indv LRK11.variant7.variant2 --indv LRK12.variant7.variant2 --indv LRK13.variant7.variant2 --indv LRK17.variant7.variant2 --indv LRK22.variant7.variant2 --gzvcf /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz --recode --out /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/mainelynx/6SV_mainelynx"

#The output then needs to be gzipped with:
bsub -q short -W 1:00 -R rusage[mem=48000] -n 1 -R span\[hosts=1\] gzip 6SV_mainelynx.vcf

#The gzvcf then needs to be recoded to ped/map format using PLINK
plink --vcf  6SV_mainelynx.vcf.gz --recode --allow-extra-chr
```
### northSLR --done
```
bsub -q short -W 1:00 -R rusage[mem=8000] -n 2 "vcftools --indv a202.variant.variant.variant2 --indv a33.variant.variant.variant2 --indv a772.variant --indv a803.variant --indv a818.variant --indv b114.variant --indv b124.variant2.variant.variant2 --indv b188.variant2.variant.variant2 --indv b554.variant2.variant.variant2 --indv c165.variant2.variant.variant2  --indv c323.variant2.variant2 --indv b13.variant2.variant2 --indv f264.variant2.variant2  --indv f457.variant3.variant2  --indv fha_024.variant3.variant2 --indv fha_042.variant3.variant2 --indv fha_043.variant3.variant2 --indv l09_007.variant3.variant2 --gzvcf /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz --recode --out /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_northSLR"

#The output then needs to be gzipped with:
bsub -q short -W 1:00 -R rusage[mem=48000] -n 1 -R span\[hosts=1\] gzip 6SV_northSLR.vcf

#The gzvcf then needs to be recoded to ped/map format using PLINK
plink --vcf  6SV_northSLR.vcf.gz --recode --allow-extra-chr
```
### bobcat --done
```
bsub -q short -W 1:00 -R rusage[mem=8000] -n 2 "vcftools --indv BOBCAT1.variant2.variant2 --indv BOBCAT2.variant2.variant2 --indv BOBCAT3.variant2.variant2 --indv BOBCAT4.variant2.variant2 --indv BOBCAT5.variant2.variant2 --indv BOBCAT6.variant2.variant2 --gzvcf /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_all/allsamples.vcf.gz --recode --out /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_all/bobcats.vcf.gz"

#The output then needs to be gzipped with:
bsub -q short -W 1:00 -R rusage[mem=48000] -n 1 -R span\[hosts=1\] gzip bobcats.vcf

#The gzvcf then needs to be recoded to ped/map format using PLINK
plink --vcf  bobcats.vcf.gz --recode --allow-extra-chr
```
### 1. A wrapper script for running SNeP chromosome by chromosome
We wrote a wrapper script in each folder /northSLR /southSLR /NFLD /mainelynx. We will then plot and compare trends
```
bsub -q long -R rusage[mem=80000] -n 1 -W 24:00 -R span\[hosts=1\] "./wrapper_mainelynx_chr_by_chr.sh" 

```
### 2. Running SNeP on bobcat 
We will run SNeP on the full suite of chromosomes using the parameters we settled on above (see Sensitivity Analysis)

bsub -q long -R rusage[mem=80000] -n 1 -W 12:00 -R span\[hosts=1\] "singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_all/bobcats.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_all/bobcats.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/bobcats_NC_044303.1 -chr scaffold_11_arrow_ctg1 -samplesize 6 -threads 8" 

### 3. Running SNeP on southSLR
bsub -q short -R rusage[mem=80000] -n 1 -W 1:00 -R span\[hosts=1\] "singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/southSLR/6SV_southSLR.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/southSLR/6SV_southSLR.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/southSLR_NC_044303.1 -chr NC_044303.1 -samplesize 40 -itemsTH 1 -threads 8" #the itemsTH flag needs to be included, always, because 500+ items per bin is not always possible

### 4. Running SNeP on northSLR
bsub -q long -R rusage[mem=80000] -n 1 -W 24:00 -R span\[hosts=1\] "singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/northSLR/6SV_northSLR.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/northSLR/6SV_northSLR.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/northSLR_NC_044303.1 -chr NC_044304.1 -maxsnp 10000 -samplesize 18 -itemsTH 1 -threads 8"

### 5. Running SNeP on NFLD
bsub -q short -R rusage[mem=80000] -n 1 -W 1:00 -R span\[hosts=1\] "singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/NFLD/6SV_NFLD.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/NFLD/6SV_NFLD.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/NFLD_NC_044303.1 -chr NC_044303.1 -samplesize 3 -itemsTH 1 -threads 8"

Plotting Trend lines in R
---

I could potentially copy the R code here, but I think I'll just publish it to RPubs and link it. 

Results
---
Here is likely our best population trend plot for lynx southSLR (south of the St Lawrence River), and Newfoundland. 
Missing northSLR because it is still running
![](https://i.imgur.com/bYswF1v.jpg)

Zoom into the most recent 100 generations
![](https://i.imgur.com/yMEiv4t.jpg)


Notes 
---
### LinkNe
SNeP and LinkNe have similar methodology. They both use linkage desequilibrium to calculate Ne in the "recent"past. You can read more about LinkNe here, but I will not be using it because it requires a linkage map inferred from pedigree data that we do not have for lynx at this time: 
https://github.com/chollenbeck/LinkNe






















# Traditional Ne Estimates
# Theta
https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.2945

# NeEstimatorv2 https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12157
 ‘java –jar./NeGUI.jar’. An example of an input data file is provided as well as a help file in.pdf and.html formats.

# Ne_StrataG.R
Script for estimating Effective Population Size (Ne) using Linkage Desequilibrium in StrataG R package.
Convert a VCF/Genind to GENEPOP format to run analyses on NeEstimator.
https://github.com/jdalapicolla/Ne_StrataG.R/blob/master/Ne_Estimation.R

# Run this tomorrow
Estimate pairwise linkage disequilibrium (LD) between all SNPs in each population. (Plink)
Use the LD estimates to estimate effective population size (Ne) in each population. (R)
Compare estimates of the census (Nc) and effective (Ne) population sizes.

################ SCRATCH PAD (do not delete please)

## Diversity Calculations -- these needs to be moved to the Heterozygosity HackMD
One of the stats we're interested in calculating is a comparison to Godoy's heterozygous SNP/Mb figure, which was originally presented by Cho et al. 2013 https://www.nature.com/articles/ncomms3433?ref=driverlayer.com#accession-codes

In the original figure, the rate of heterozygous SNPs was calculated across Panthera species. the heterozygous SNV rate (y axis) was calculated by dividing the total number of heterozygous SNVs by genome size. We do the following to acquire the total count of heterozygous SNVs, then divide by the genome size (2.41)
![](https://i.imgur.com/jSfzIYp.png)

1. Start an interactive session
bsub -q interactive -R rusage[mem=24000] -n 1 -W 30 -Is bash
2. Navigate to our unfiltered VCF
3. Use the following command to count the number of homozygous major(0|0),homozygous minor(1|1) and heterozygous(1|0,0|1) alleles for each position of a chromosome, across all chromosomes. 
zcat 6SV.vcf.gz|awk -v OFS="\t" '$0 !~ "^#" {hom_ref = 0; hom_alt = 0; het_01 = 0; het_10=0; for(i=10;i<=NF;i++) { if($i ~ /0\|0/) hom_ref++; else if($i ~ /1\|1/) hom_alt++;  else if($i ~ /0\|1/) het_01++; else if($i ~ /1\|0/) het_10++; } print $1, $2, hom_ref, hom_alt, het_01, het_10}'
> hetcount

Not sure if that was correct because it looks like that heterozygosity estimate was really high (0.0013) compared to other species. Do we also have an estimate from VGP? Check that as well 

### Calculating site by site heterozygosity using -freqx
We will also eventually use this in our outlier analysis. Similarly to how we calculated FST

### Load Dependencies
module load anaconda2/4.4.0
module load plink/1.90b6.9

### PLINK
plink --make-bed --file /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs --freqx --out output_freqx2 --allow-extra-chr 

Now we can take the output and hopefully calculate average heterozygosity by position (to visualize LROH) as in the Godoy paper

rerunning all together. debating whether we should separate out populations? using plink?

## Comparing heterozygosity rates among species
We tried to replicate what was presented in Godoy and Cho, but they are not very specific about what method was used to calculate average het rate. I did my best to replicate their methods, and the het rate estimated 0.013 seems high compared to other species. Let's try something that will be more robust using jellyfish and genomescope. This is what we use for the VGP assemblies, so I know it will be a fair method of comparison for the Canada lynx reference genome and other reference genomes (mainly felids). 

There is a great beginners tutorial on GenomeScope here: 
https://github.com/schatzlab/genomescope

We start on the cluster by downloading the necessary dependencies

Working directory is: 
/project/uma_lisa_komoroske/Tanya/scripts/genomescope

### Load Dependencies
module load anaconda2/4.4.0 
module load jellyfish/2.2.3 

Jellyfish will accept .fasta files, which is good news for us!

1. Counting kmers. We use 31-mers here  
bsub -q long -R rusage[mem=16000] -n 4 -W 6:00 "jellyfish count -m 31 -o mLynCan4v1p -c 3 -s 10000000 -t 32 /project/uma_lisa_komoroske/Tanya/download/refs/mLynCan4_v1.p/GCF_007474595.1_mLynCan4_v1.p_genomic.fa" #LSF mem limit with 8G, trying 16G, 

bsub -q short -R rusage[mem=16000] -n 2 -W 1:30 "jellyfish count -m 31 -o pantig -c 3 -s 10000000 -t 32 /project/uma_lisa_komoroske/Tanya/download/refs/felidae/GCA_000464555.1_PanTig1.0_genomic.fa"


2. All of the output files named mLynCan4v1p need to be merged using: 
jellyfish merge -o mLynCan4v1p.jf mLynCan4v1p\_*

3. Then export the kmer count histogram
bsub -q short -R rusage[mem=16000] -n 2 -W 1:30 "jellyfish histo -t 10 mLynCan4v1p.jf > mLynCan4v1p.histo"

We've downloaded all of the felid genomes to /project/uma_lisa_komoroske/Tanya/download/refs/felidae
GCA_000181335.4_Felis_catus_9.0_genomic.fna.gz
GCA_000464555.1_PanTig1.0_genomic.fna.gz
GCA_001857705.1_PanPar1.0_genomic.fna.gz
GCA_003327715.1_PumCon1.0_genomic.fna.gz
GCA_003709585.1_Aci_jub_2_genomic.fna.gz
GCA_004023805.1_PanOnc_v1_BIUU_genomic.fna.gz
GCA_004023925.1_FelNig_v1_BIUU_genomic.fna.gz
GCA_005406085.1_Prionailurus_bengalensis_euptilurus_v01_genomic.fna.gz
GCA_007474595.1_mLynCan4_v1.p_genomic.fna.gz
GCA_008795835.1_PanLeo1.0_genomic.fna.gz
GCA_900661375.1_LYPA1.0_genomic.fna.gz


Let's try Jellyfish on each of them 
If you use JELLYFISH in your research, please cite:

Guillaume Marcais and Carl Kingsford, A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics (2011) 27(6): 764-770 (first published online January 7, 2011) doi:10.1093/bioinformatics/btr011

Jellyfish/Genomescope does not appear to be working on the reference genomes, so I'm going to try something else

1. calculate the mean heterozygosity in R (by population)
2. forget about comparisons. If anything try to compare to iberian lynx metric from paper. What method?

I tried running SNeP1.1 with the following argument:
./SNeP1.1 -ped ~/file.ped
My .map file has several chromosomes but I get the following message from the program:
Number of Chromosomes detected: 1
Is there an error?
Also, it did end up creating a .NeAll file. However, the Ne values are too low. Are the Ne values scaled by a factor?


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `Ne` `SNeP` `Effective Population Size` `Demographic Trends` `Contemporary Demographic Reconstruction`
