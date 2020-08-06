---
title: 'EEMS'
disqus: hackmd
---

Estimated effetive migration surface (EEMS) using genome-wide SNPs
===

[![hackmd-github-sync-badge](https://hackmd.io/lgVVPmECSuCYcML6ieOFTQ/badge)](https://hackmd.io/lgVVPmECSuCYcML6ieOFTQ)

from the EEMS manual: 
This is an implementation of the EEMS method for analyzing and visualizing spatial population structure from geo-referenced genetic samples (Canada lynx). EEMS uses effective migration to model the relationship between genetics and geography, and outputs an estimated effective migration surface (EEMS). It's basically a visual representation of population structure that can highlight fine-scale population structure, and potential regions of higher-than average and lower-than-average historic gene flow. In other words, it **highlights features of the landscape where genetic similarity decays faster than  expectations under isolation-by-distance**. 

EEMS constructs a triangular grid over available habitat. In this situation we are using a shapefile encompassing the USFWS designation of critical habitat (US) and the species distribution (Canada) as "available habitat". EEMS then assigns each geo-referenced individual to its closest vertex (deme) and estimates a migration rate for every edge and a diversity rate for every deme in the grid. 

## Objective: 
Use EEMS to visualize patterns in 1. estimated effective migration rate 2. genomic diversity across the study area

## Additional reading 
* The EEMS paper: http://www.nature.com/ng/journal/v48/n1/full/ng.3464.html
* Example implementation: https://osf.io/xbf5n/ 
* Author git Repo: https://github.com/dipetkov/eems
* My Git: https://github.com/ECOtlama/project_canada_lynx_wgs/tree/master/EEMS

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

The VCF we want to use is: /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_HighQualSites_processed.vcf.gz  which includes all individuals in our three lynx populations (NFLD, NSLR, SSLR). 

Convert the VCF to binary (bin/fam) format using PLINK.
```
plink --vcf 6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz --allow-extra-chr '0'
```

### Workaround for plink -allow-extra-chr
plinkio is having trouble with an "unrecognized chromosome code" because we're using a non-human genome with more chromosomes than expected and chromosome names that differ from the expected names e.g. those in hg19 (that start with a letter, not a number). You can expect to have this issue with downstream analyses if you're using a VGP genome or a draft genome with scaffolds rather than chromosomes.

The solution is this: 
#### 1. use plink 1.90 to convert VCF to ped/map format
```
plink --vcf  6SV_unfiltered_SNPs.vcf.gz --recode --allow-extra-chr
```
#### 2. unload plink and module load plink/1.07
#### 3. use plink 1.07 to convert ped/map format to binary (bed/bim/fam). 
```
plink 6SV_unfiltered_SNPs --make-bed
```
#### 4. proceed with bed2diffs using the 6SV_SNPs bed/bim/fam files
 

### bed2diffs -- this runs pretty quick
working directory is /project/uma_lisa_komoroske/Tanya/scripts/EEMS
Usage: 
```
bsub -q short -W 0:30 -R rusage[mem=1000] -n 2 -R span\[hosts=1\] "bed2diffs_v1 --bfile 6SV_unfiltered_SNPs --nthreads 2"
```
bed2diffs generates two files. The matrix of average pairwise differences is written to a text file without row names, column names or comments, with extension .diffs. The size of the matrix is NxN where N is the number of samples; there are only 0s on the main diagonal. The order of the samples in the diffs matrix is the same as in the fam file, but in any case, the order is explicitly written to a text file with one sample per line, with extension order.

### Load the dissimilarity matrix into R
```
diffs<- read.table("/project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_SNPs.diffs")
```
Note: you should see a matrix with 0.00 along the diagonal

### datapath.coord 
.coord is a list of sample coordinates (2 coords per sample, one sample per line)
We have generated these coords for previous mapping, and had to convert them from NAD83 zone 19 to WGS84 in R. See the "study area mapping" RMarkdown for directions on that.

### datapath.outer
.outer is a list of habitat coordinates as a sequence of vertices that outline a closed polygon. 

This was a pain in the ass. I tried to coerce various SpatialPolygonsDataFrames to lines and then points, which looked great. But EEMS .outer format requires the points to be listed in counterclockwise order. So ultimately, after many hours of toiling we ended up drawing points along our available habitat shapefile and formatting it according to the specifications in EEMS manual. 
The final .outer file is called /Users/tanyalama/Box/project_canada_lynx_wgs/R_canada_lynx_wgs/EEMS/IUCN_redlist_Canada_lynx_spp_distribution/EEMS_outer_arcgis_pts.csv

### Edit the parameter file params-chain1.ini
We changed the relevant fields in the first params file (see manual). Namely, the paths and the number of individuals and nSites (n SNPs -- this is listed in the plink conversion output, make note). 

## Run at least 3 rounds of EEMS (on the command line) 
Run EEMS with params-chain1.ini, -chain2.ini, -chain3.ini. The argument --seed is optional;it can be specified in the parameter file or in the command, and if not specified, the seed is randomly assigned. (try this in a interactive session if it's your first time!)
```
bsub -q long -W 8:00 -R rusage[mem=48000] -n 1 -R span\[hosts=1\] "runeems_snps --params ./src/params-chain1.ini"
```

# Visualize EEMS results
Finally, the EEMS results can be visualized with the function eems.plots defined in the R package rEEMSplot. The package is not on CRAN, so install it from source instead. (The plotting code is in the directory "plotting".)

Visit [my github to see the RMarkdown](https://raw.githubusercontent.com/ECOtlama/project_canada_lynx_wgs/master/EEMS/visualization_mapping_EEMS_results.Rmd)
I describe local installation of the rEEMSplot package and steps for visualization of the results. 

# Results
## Estimated effective migration rate
The between-demes component characterizes the expected genetic dissimilarity between two individuals from distinct demes and is a function of the effective migration rate. It represents dissimilarity that is due to the spatial structure of the population and is not a consequence of the local diversity in the two demes. Colors in the estimated effective migration surface correspond to local deviations from isolation by distance: in particular, effective migration is low **(orange) in geographic regions where genetic similarity decays quickly**.

![](https://i.imgur.com/8pF9SZW.png)

## Estimated effective diversity rate
EEMS accounts for local variation in genetic diversity when modeling spatial patterns in genetic dissimilarity. Genetic dissimilarity can therefore be partitioned into two components 1. within-demes (diversity) and 2. between demes (migration). The "within-demes" component, W, characterizes expected genetic dissimilarity between individuals from the same deme, and is a function of the effective diversity rate. Diversity parameters can be visualized on the log10 scale, so that **blue indicates higher than average diversity and orange indicates lower than average diversity**. 


![](https://i.imgur.com/4KpOQqT.jpg)

## Interpretation
EEMS confirms the results from our other analyses (PCA, sNMF) indicating that the St. Lawrence River is not a complete barrier to gene flow. We've detected multiple bi-directional dispersal events north-south of the St. Lawrence River and one dispersal event between Newfoundland and mainland Canada. However, genetic similarity north-south of the St. Lawrence River decays faster than expected under isolation-by-distance. 

This connectivity pattern relates directly to metrics of genetic diversity, which are well-supported by genomic diversity rates output by EEMS. We see (above) that lower-than-expected historic gene flow between island populations and mainland Canada has eroded genetic diversity in the Cape Breton Island and Newfoundland populations. Despite evidence of connectivity, genetic diversity south of the St Lawrence River (including Maine) is lower than that north of the St. Lawrence River. 

See the Status Report on lynx in Nova Scotia for further interpretation about which populations have been in decline/extirpated/protected etc in the last 50 years.

## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `EEMS` `migration rate` `migration` `landscape genomics` `estimated effective migration surface`
