---
title: 'EEMS'
disqus: hackmd
---

EEMS estimated effetive migration surface using genome-wide SNPs
===
This is an implementation of the EEMS method for analyzing and visualizing spatial population structure from geo-referenced genetic samples (Canada lynx). EEMS uses effective migration to model the relationship between genetics and geography, and outputs an estimated effective migration surface (EEMS). It's basically a visual representation of population structure that can highlight potential regions of higher-than average and lower-than-average historic gene flow. In other words, it **highlights features of the landscape where genetic similarity decays faster than expected by null expectations under isolation-by-distance**. 

Please consider reading http://www.nature.com/ng/journal/v48/n1/full/ng.3464.html
and https://osf.io/xbf5n/ 

## Table of Contents

[TOC]

## Load Dependencies
EEMS was loaded as a module by MGHPCC because dependencies including plinkio were troublesome to install without permissions.
```
module load anaconda2/4.4.0
module load eems/201809
module load plink/1.07
``` 
## Input Files
### datapath.diffs
.diffs is a matrix of average pairwise differences (computed with bed2diffs). bed2diffs reads genotypes in PLINK binary format and computes the matrix of average pairwise differences as an input for EEMS

The VCF we want to use is: /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz
which includes all individuals in our three lynx populations (NFLD, NSLR, SSLR). We need to make sure the VCF is available in binary (bin/fam) format from PLINK.
```
plink --vcf 6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz --allow-extra-chr '0'
```


Note: I am considering using SNPSift to create a VCF file that is only intergenic (neutral) SNPs for this...

### Workaround for plink -allow-extra-chr
plinkio is having trouble with an "unrecognized chromosome code" because we're using a non-human genome with more chromosomes than expected and chromosome names that differ from the expected names in hg19 (must start with a letter). You can expect to have this issue with downstream analyses (e.g. ADMIXTURE) if you're using a VGP genome or a draft genome with scaffolds rather than chromosomes.
The solution is this: 
#### 1. use plink 1.90 to convert VCF to ped/map format
```
plink --vcf  6SV_unfiltered_SNPs.vcf.gz --recode --allow-extra-chr
```
#### 2. remove plink and module load plink/1.07
#### 3. use plink 1.07 to convert ped/map format to binary (bed/bim/fam). 
```
plink 6SV_unfiltered_SNPs --make-bed
```
#### 4. proceed with bed2diffs using the 6SV_unfiltered_SNPs bed/bim/fam files
 

### bed2diffs -- this runs pretty quick
working directory is /project/uma_lisa_komoroske/Tanya/scripts/EEMS
Usage: 
```
bsub -q short -W 0:30 -R rusage[mem=1000] -n 2 -R span\[hosts=1\] "bed2diffs_v1 --bfile 6SV_unfiltered_SNPs --nthreads 2" #8377918
```
bed2diffs generates two files. The matrix of average pairwise differences is written to a text file without row names, column names or comments, with extension .diffs. The size of the matrix is NxN where N is the number of samples; there are only 0s on the main diagonal. The order of the samples in the diffs matrix is the same as in the fam file, but in any case, the order is explicitly written to a text file with one sample per line, with extension order.
### Load the dissimilarity matrix into R
```
diffs<- read.table("/project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.diffs")
diffs
```
Note: you should see a matrix with 0.00 along the diagonal

### datapath.coord 
.coord is a list of sample coordinates (2 coords per sample, one sample per line)
We have generated these coords for previous mapping, and had to convert them from NAD83 zone 19 to WGS84 in R. See the "study area mapping" RMarkdown for directions on that.

### datapath.outer
.outer is a list of habitat coordinates as a sequence of vertices that outline a closed polygon. 
This was a pain in the ass to get. I tried to coerce various SpatialPolygonsDataFrames to lines and then points, which worked great. But EEMS .outer format requires the points to be listed in order, counterclockwise. So ultimately, after many hours of toiling we ended drawing points along the Canada lynx species critical habitat shapefile and formatting it according to the specifications in EEMS manual. 
The final .outer file is called /Users/tanyalama/Box/project_canada_lynx_wgs/R_canada_lynx_wgs/EEMS/IUCN_redlist_Canada_lynx_spp_distribution/EEMS_outer_arcgis_pts.csv

### Edit the parameter file params-chain1.ini
We changed the relevant fields in the first params file. Namely, the paths and the number of individuals and nSites (n SNPs). 

## Run EEMS (on the command line) 
Run EEMS (on the command line)with three different random seeds. The argument --seed is optional;it can be specified in the parameter file as well, and if not specified, the seed is randomly assigned.
We'll try this in an interactive session 
bsub -q interactive -R rusage[mem=24000] -n 1 -W 30 -Is bash
```
bsub -q long -W 8:00 -R rusage[mem=48000] -n 1 -R span\[hosts=1\] "runeems_snps --params ./src/params-chain1.ini"
./runeems_snps --params params-chain2.ini --seed 456
./runeems_snps --params params-chain3.ini --seed 789
```

# Visualize EEMS results
Finally, the EEMS results can be visualized with the function eems.plots defined in the R package rEEMSplot. The package is not on CRAN, so install it from source instead. (The code is in the directory plotting.)

## Part 1: Install rEEMSplots
## Check that the current directory contains the rEEMSplots source directory
if (file.exists("./rEEMSplots")) {
  install.packages("rEEMSplots", repos = NULL, type = "source")
} else {
  stop("Move to the directory that contains the rEEMSplots source to install the package.")
}

## Part 2: Generate graphics
library(rEEMSplots)

mcmcpath = "./data/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes200-simno1"
plotpath = "./plot/barrier-schemeX-nIndiv300-nSites3000-EEMS-nDemes200-simno1-rEEMSplots"

eems.plots(mcmcpath, plotpath, longlat = TRUE)


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `EEMS` `migration rate` `migration` `landscape genomics` `estimated effective migration surface`
