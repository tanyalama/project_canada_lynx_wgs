---
title: 'LFMM'
disqus: hackmd
---

LFMM latent factor mixed models
===
LFMMs belong to a broad class of statistical models for association studies including genome-wide association studies (GWAS), epigenome-wide association studies (EWAS), and genome-environment association studies (GEAS) among others. LFMMs were used in (Frichot et al. 2013) for testing correlations between loci and environmental variables. Other implementations of LFMMs are available in the R packages lfmm, LEA (Frichot et al. 2013), sva (Leek and Storey 2007), and cate (Wang et al. 2017).

The R package lfmm implements new algorithms for parameter estimation in latent factor mixed models (LFMM). The new methods are computationally efficient, and provide statistically optimal corrections resulting in improved power and control for false discoveries (i.e. corrections for multiple testing). The package lfmm provides two main functions for estimating latent confounders (or factors): lfmm_ridge and lfmm_lasso. Those functions are based on optimal solutions of regularized least-squares problems. A short tutorial provides brief examples on how the R packages lfmm can be used for fitting latent factor mixed models and evaluating association between a response matrix (SNP genotype or methylation levels) and a variable of interest (phenotype or exposure levels) in genome-wide (GW), genome-environment (GE), epigenome-wide (EW) association studies. Corresponding software is available at the following url https://bcm-uga.github.io/lfmm/.

## Table of Contents

[TOC]

---
# Part 1: Prepare SNP data

## Input Files
LFMM requires two input files: SNP data in 012 format with missing data coded as 9, and a table with environmental variables of interest in which rows are samples and column are variables. Neither can have headers or row names (the original order will be preserved in the output). 

We need to use the LD-pruned SNP set from the LEA package. I find this interesting and didn't get a good justification for using LD-pruned SNPs for LFMM. Certainly SNPs which could be outliers were dropped in LD-pruning, and now won't be assessed for association with environmental variables (but ok...)

### SNP data
We'll use vcftools to convert our LD-pruned SNPset in snp x major mode to  012 format (homozygous, heterozygous, other homozygous). 
Working directory is /project/uma_lisa_komoroske/Tanya/scripts/LFMM
We'll try this in an interactive session 
```
bsub -q interactive -R rusage[mem=24000] -n 1 -W 30 -Is bash

vcftools --gzvcf 6SV.vcf.gz  --012 --out 6SV012
```
The output file will be filename.012

### Encode missing data as -9
Use the sed command to encode missing data (-1 in 012 format) as 9, and then use the cut command to retain everything after column 2. Save the object as snp.lfmm. This is what we'll use as the SNP input for LFMM in R.
```
sed 's/-1/9/g' 6SV012.012 | cut -f2- > snp.lfmm
```
View the result and move the file to /project/uma_lisa_komoroske/Tanya/scripts/LFMM folder. 

### WorldClim climate data (current)
We will use current and projected climate data from WorldClim. Re: data choice, I explored other, potentially more appropriate data sources (e.g. CHELSA) but because our study area expands into Canada, US-centric climate data wasn't an option and I had to opt for a lower-resolution global database like WorldClim. This has been used in a lot of published works, and Brenna says that it's not something I would likely be criticized for from reviewers.

Download climate data layers to the /climate_data folder in GeoTiff format. We will prepare the climate data for lfmm in R. 

Note that instead of using climate variables directly, we will use principal components of the climate variables as uncorrelated, synthetic climate variables. When dealing with correlated climate variables (as is typically the case), this strategy may be preferrable. This also reduces the number of multiple tests from nSNP * 5 to nSNP * 2 if we use just the first two PC's.

See RMarkdown for details on how to format raw WorldClim data for R

### LFMM
See RMarkdown for details on running LFMM

### Compile geneset of interest from LFMM results 
LFMM determined which loci are outliers as well as their association with available climate variables.

See RMarkdown for details on this.

Now that we have a small list of outlier genes, we can subset our VCF using positions identified by LFMM. We will pull the SNPs at those positions, annotate them and run some simple statistics such as allele frequency on each locus. 

### Subset the VCF at select positions
Use nano to write a list of positions (from lfmm results) that includes (without headers):
CHR tab POS 
We named this positions.txt

Then subset the 6SV VCF at the given loci: 
```
zcat /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV.vcf.gz | vcftools --gzvcf /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV.vcf.gz --positions lfmm_results_two_factor_positions.txt --recode --out lfmm_results_two_factor_snpset.vcf
```
4811/53265 SNPS meaningful associations w/ climate 

### Calculate allele frequency at each locus
```
vcftools --vcf subset_snp_list.vcf --freq --out subset_snp_list_allele_freq 
```
### Annotate our snpset with SNPSift
module load anaconda2/4.4.0
module load  java/1.8.0_171
We've locally installed SNPEff and the mLynCan4 reference genome annotation, so usage is: 
```
java -Xmx10G -jar /project/uma_lisa_komoroske/bin/snpEff/snpEff.jar mLynCan4 ./lfmm_results_snpset.vcf > snpset_ann.vcf
```
### Extract gene names and describe function 
```
java -jar /project/uma_lisa_komoroske/bin/snpEff/SnpSift.jar extractFields lfmm_results_snpset.ann.vcf "ANN[*].GENE:" | awk -F"|" '{print $4}'

java -jar /project/uma_lisa_komoroske/bin/snpEff/SnpSift.jar extractFields LD_pruned_lfmm_results_all_ann.vcf CHROM POS REF ALT "ANN[*].GENE:" | awk -F'[\t|]' '{print $1,$2,$3,$4,$6,$8}' OFS="\t" > snpsift_LD_pruned_lfmm_results_all_variant_type.txt

java -jar /project/uma_lisa_komoroske/bin/snpEff/SnpSift.jar extractFields LD_pruned_lfmm_results_all_ann.vcf "ANN[*].GENE:" | awk -F"|" '{print $4}' OFS="\t" > snpsift_LD_pruned_lfmm_results_all
```

### Our top two genes, associated with three/four climatic variables are:
[TFRC](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TFRC)
[ABCA3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ABCA3&keywords=ABCA3)

### Download projected climate data 

We're going to use [SSP2 (Shared Socioeconomic Pathway)](https://www.carbonbrief.org/explainer-how-shared-socioeconomic-pathways-explore-future-climate-change). SSP2 represents a “middle of the road” scenario historical patterns of development are continued throughout the 21st century. In other words, "business as usual". 

## Refining LFMM analysis
1. Selecting better BioClim variables for lynx
We can see the full list of BioClim variables [here](https://www.worldclim.org/data/bioclim.html). I've selected: 

"tseas" BIO4 = Temperature Seasonality (standard deviation ×100)

"tmin"  BIO6 = Min Temperature of Coldest Month

"tmean" BIO11 = Mean Temperature of Coldest Quarter

"pseas" BIO15 = Precipitation Seasonality (Coefficient of Variation)

"pcol"  BIO19 = Precipitation of Coldest Quarter

Run K=3. Download and run climate projection comparison. Calculate allele frequencies. Use the snpsift_circ_therm_genelist_6SV.txt for LFMM 

---
# Part 3
I'm going to run LFMM on a consensus set of outlier SNPs identified via different methods. The set we'll run with LFMM is a consensus set of 88 outlier SNPs identified in 54 genes associated with thermal tolerance and circadian rhythm in lynx. Let's see how many of these SNPs are associated with climate variables.

### Convert outlier_loci_circ_therm_genelist_6SV.vcf to 012 format
```
vcftools --vcf /project/uma_lisa_komoroske/Tanya/scripts/outlier/outlier_loci_circ_therm_genelist_6SV.vcf  --012 --out outlier_loci_circ_therm_genelist_6SV #61 individuals and 88 loci
```
### Encode missing data as -9
Use the sed command to encode missing data (-1 in 012 format) as 9, and then use the cut command to retain everything after column 2. 
```
sed 's/-1/9/g' 6SV012.012 | cut -f2- > snp.lfmm
```
### Climate Variables
We'll use the lynx_clim.env input, which includes five climate BioClim variables we selected (see above). 

### Run LFMM in R
LFMM needs to be run on the cluster via RStudio OnDemand. Make sure you setwd to /project/uma_lisa_komoroske/Tanya/scripts/LFMM/outlier_loci_circ_therm and use the outlier_loci_circ_therm_snp.lfmm and lynx_clim.env input files in the /outlier_loci_circ_therm folder.

### Annotate with SNPEff
Use the position data you've extracted for SNPs with a significant (<0.01, <0.05) association with climate variables to decipher which genes they are associated with. 

module module load anaconda2/4.4.0
module load  java/1.8.0_171
module load vcftools/0.1.14

bsub -q interactive -R rusage[mem=24000] -n 1 -W 30 -Is bash

java -Xmx20G -jar /project/uma_lisa_komoroske/bin/snpEff/snpEff.jar mLynCan4 ./outlier_loci_circ_therm_genelist_6SV.vcf > outlier_loci_circ_therm_genelist_6SV.ann.vcf

Post-Processing snpEff Annotations using snpSift
And then pull the two outlier SNPs out of that annotated VCF and see which genes they associate with

Then subset the 6SV VCF at the given loci: 
```
zcat ./outlier_loci_circ_therm_genelist_6SV.vcf | vcftools --vcf ./outlier_loci_circ_therm_genelist_6SV.vcf --positions positions_outlier_loci.txt --out lfmm_results_outlier_loci_circ_therm_0.05_significance.vcf

zcat ./mLynCan4v1p_filtered_reduced.vcf.gz | vcftools --gzvcf mLynCan4v1p_filtered_reduced.vcf.gz --positions future_LD_pruned_lfmm_results_all_positions --recode --out future_outlier_LD_pruned_snps.vcf
```
### Annotate with SNPEff
```
java -Xmx20G -jar /project/uma_lisa_komoroske/bin/snpEff/snpEff.jar mLynCan4 ./lfmm_results_outlier_loci_circ_therm_0.05_significance.vcf > lfmm_results_outlier_loci_circ_therm_0.05_significance_ann.vcf

java -Xmx20G -jar /project/uma_lisa_komoroske/bin/snpEff/snpEff.jar mLynCan4 ./future_outlier_LD_pruned_snps.vcf > future_outlier_LD_pruned_snps_ann.vcf
```
### Extract gene names and describe function 
```
java -jar /project/uma_lisa_komoroske/bin/snpEff/SnpSift.jar extractFields lfmm_results_outlier_loci_circ_therm_0.05_significance_ann.vcf "ANN[*].GENE:" | awk -F"|" '{print $4}'

### Calculate allele frequency for the population at each locus
vcftools --vcf lfmm_results_outlier_loci_circ_therm_0.05_significance_ann.vcf --freq --out lfmm_results_outlier_loci_circ_therm_0.05_significance_freq 
```

# Part 4: 1650 SNPs
Let's run this one more time with a larger set of "consensus" SNPs within our ~50+ genes of interest. This vcf is called 
1650_snps_circ_therm_57_genes_6SV.vcf

vcftools --vcf 1650_snps_circ_therm_57_genes_6SV.vcf  --012 --recode --out 1650_snps_circ_therm_57_genes_6SV

### Convert to 012 format with vcftools

# Results
**1970-2000 climate conditions**
9.03% of SNPS (4811/53k) were significantly <0.01 associated with climate variables
**2061-2080 climate conditions**
8.6% of SNPs (4603/54k) are significantly <0.01 associated with climate variables

Yet to say whether they are the same genes/ gene overlap

## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `LFMM` `adaptive genomics` `climate change` `WorldClim` `selection` `adaptation` `R`
