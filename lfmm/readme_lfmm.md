---
title: 'LFMM'
disqus: hackmd
---

LFMM latent factor mixed models
===


[![hackmd-github-sync-badge](https://hackmd.io/QB6ERAK0Tbe1ZeOmRUYO7g/badge)](https://hackmd.io/QB6ERAK0Tbe1ZeOmRUYO7g)
LFMMs belong to a broad class of statistical models for association studies including genome-wide association studies (GWAS), epigenome-wide association studies (EWAS), and genome-environment association studies (GEAS) among others. LFMMs were used in (Frichot et al. 2013) for testing correlations between loci and environmental variables. Other implementations of LFMMs are available in the R packages lfmm, LEA (Frichot et al. 2013), sva (Leek and Storey 2007), and cate (Wang et al. 2017).

The R package lfmm implements new algorithms for parameter estimation in latent factor mixed models (LFMM). The new methods are computationally efficient, and provide statistically optimal corrections resulting in improved power and control for false discoveries (i.e. corrections for multiple testing). The package lfmm provides two main functions for estimating latent confounders (or factors): lfmm_ridge and lfmm_lasso. Those functions are based on optimal solutions of regularized least-squares problems. A short tutorial provides brief examples on how the R packages lfmm can be used for fitting latent factor mixed models and evaluating association between a response matrix (SNP genotype or methylation levels) and a variable of interest (phenotype or exposure levels) in genome-wide (GW), genome-environment (GE), epigenome-wide (EW) association studies. Corresponding software is available at the following url https://bcm-uga.github.io/lfmm/.

## Table of Contents

[TOC]

---

## Input Data
LFMM requires two input files: SNP data in 012 format with missing data coded as 9, and a table with environmental variables of interest in which rows are samples and column are variables. Neither can have headers or row names (the original order will be preserved in the output). 

We need to use the LD-pruned SNP set from the LEA package. I find this interesting and didn't get a good justification for using LD-pruned SNPs for LFMM. Certainly SNPs which could be outliers were dropped in LD-pruning, and now won't be assessed for association with environmental variables (but ok...)

### SNP data
We'll use vcftools to convert our LD-pruned SNPset in snp x major mode to  012 format (homozygous, heterozygous, other homozygous). 
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
View the result and move the file to /LFMM folder. 

### WorldClim climate data (current)
We will use current and projected climate data from WorldClim. Re: data choice, I explored other, potentially more appropriate data sources (e.g. CHELSA) but because our study area expands into Canada, US-centric climate data wasn't an option and I had to opt for a lower-resolution global database like WorldClim. This has been used in a lot of published works, and Brenna says that it's not something I would likely be criticized for from reviewers.

Download climate data layers to the /climate_data folder in GeoTiff format. We will prepare the climate data for lfmm in R. 

Note that instead of using climate variables directly, we will use principal components of the climate variables as uncorrelated, synthetic climate variables. When dealing with correlated climate variables (as is typically the case), this strategy may be preferrable. This also reduces the number of multiple tests from nSNP * 5 to nSNP * 2 if we use just the first two PC's. I discussed this with Brenna and she didn't have a preference and said that as long as multicollinearity and multiple testing are addressed, it should be OK to do either, and methods like GDM and RDA actually just use the variables directly, rather than PCs. We've run it both ways.

See LFMM.Rmd for details on how to format raw WorldClim data for R

### WorldClim Variables
We can see the full list of BioClim variables [here](https://www.worldclim.org/data/bioclim.html). I've selected: 

"tseas" BIO4 = Temperature Seasonality (standard deviation ×100)

"tmin"  BIO6 = Min Temperature of Coldest Month

"tmean" BIO11 = Mean Temperature of Coldest Quarter

"pseas" BIO15 = Precipitation Seasonality (Coefficient of Variation)

"pcol"  BIO19 = Precipitation of Coldest Quarter

---
## LFMM
See LFMM.Rmd for details on running LFMM

---
### Compile geneset of interest from LFMM results 
LFMM determined that 9% of SNPs were significantly associated with climate and met our threshold for large effect size as well. The positions of those outliers were compiled in R, and we'll proceed with functional annotation.

See LFMM.Rmd for details on this.

### Subset the VCF at select positions
Use nano to write a positions matrix from LFMM_results that includes (without headers):
CHR tab POS 

We named this positions.txt

Then subset the 6SV.vcf at the given positions: 
```
zcat /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV.vcf.gz | vcftools --gzvcf /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV.vcf.gz --positions positions.txt --recode --out subset_snpset.vcf
```
## Allele frequency
### Calculate allele frequency at each locus
We thought it would be interesting to calculate allele frequency at each locus. In particular, I'm interested in the genes that are associated with climate and functionally relate to circadian rhythm (I have a short list of them).
```
vcftools --vcf subset_snp_list.vcf --freq --out subset_snp_list_allele_freq 
```
---
## Functional annotation with SNPEff and SNPSift
We'll now annotate the SNPs that have been identified as significantly associated with climate variables. We have access to the mLynCan4 reference genome and annotation, and that's what will serve for functional annotation here: 

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
---

## Results (current climatic conditions)
#### 4811/53265 or ~9% of SNPs have meaningful associations with at least one climate variable after multiple testing correction

#### After annotation, these SNPs were located within 1525 genes
#### Our top two genes associated with three/four climatic variables under current conditions are:
[TFRC](https://www.genecards.org/cgi-bin/carddisp.pl?gene=TFRC)
[ABCA3](https://www.genecards.org/cgi-bin/carddisp.pl?gene=ABCA3&keywords=ABCA3)

#### Seven genes from our circadian rhythm gene list were associated with at least one climatic variable under current conditions:
These genes are part of the [Circadian Entrainment Superpath](https://pathcards.genecards.org/card/circadian_entrainment): 

* RORC
* RORA
* NPAS2
* NPAS3
* CSNK1A1
* PRKCB
* PTPRC

### Allele frequency
We looked at the distribution of homref/het/homalt individuals between our three populations (North of the St. Lawrence River, South of the St. Lawrence River, and Newfoundland). We plotted the allele frequency for each population at an intronic "modifier" SNP in NPAS3: 

![](https://i.imgur.com/lFAlhgG.png)
**Figure**: white = homozygous reference; 
black = heterozygous; 
grey = homozygous alternate. Pie charts are layed over WorldClim bioclimatic variable Bio06(minimum temperature of the coldest month for 1970-2000).

We can see that the population North of the St Lawrence River has the largest number of heterozygous and homozygous alternate individuals at this locus. All individuals on the island population of Newfoundland were homozygous alternate at this SNP. This finding aligns well with our EEMS genomic diversity findings, which suggested that genomic diversity is poorest in Newfoundland and richest North of the St. Lawrence River. 

---
## LFMM under projected future (2061-2080) conditions

We used WorldClim projected bioclimatic variables under [SSP2 (Shared Socioeconomic Pathway)](https://www.carbonbrief.org/explainer-how-shared-socioeconomic-pathways-explore-future-climate-change). SSP2 represents a “middle of the road” scenario for climate projections, under which historical patterns of development are continued throughout the 21st century. In other words, this is the "business as usual" scenario.

We ran the same LFMM analysis for associations between our SNPset and future climatic conditions, as above/ seen in the LFMM.Rmd script. Here is a brief summary of the results.

## Results (future climatic conditions)
**1970-2000 climate conditions**
9.03% of SNPS (4811/53k) were significantly <0.01 associated with at least one climate variable

**2061-2080 climate conditions**
8.6% of SNPs (4603/53k) were significantly <0.01 associated with at least one climate variable
#### Two genes from our circiadian rhythm gene list *remain* associated with at least one climatic variable under future conditions.
* NPAS2
* NPAS3

---

## What's next? 
1. Maybe we'll compare the % of genes in functional groups for current and future conditions and make some more comparisons about overlap
2. Futher interpret the role/function of the circadian rhythm genes we identified that have meaningful associates with climate
3. Compelement with GDM or GF to describe "genetic offset"? We've already run GDM, but I think GF addresses genetic offset more directly and is more easily interpretable...
4. I need to track down the code and results from the time I ran this on PC components instead of individual climatic variables (from my recollection the results were not too different)
5. Add mapPies code (on desktop) to LFMM.Rmd on github
---

## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `LFMM` `adaptive genomics` `climate change` `WorldClim` `selection` `adaptation` `R`
