---
title: 'Population Structure Analyses Using SNPRelate'
disqus: hackmd
---

Population Structure Analyses Using SNPRelate
---

This script is best run in the cluster via interactive mode or the RStudio OnDemand interface. It cannot be run locally due to the volume of the data. 

This script will do the last couple rounds of filtering (LD-based SNP pruning and biallelic SNPs only) and then proceed into analyses of population structure and population admixture. 

These analyses will inform chapter 2 re: whether dispersal events connect populations of lynx North and South of the St. Lawrence River. 

Table of Contents
---
[TOC]

### Prepare our workspace and load required packages
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
setwd("/home/tl50a")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install("ape")
require(ape)
BiocManager::install("vcfR") #trying this some >no to all >non zero exit status
require(vcfR) #not working but oh well
require(SNPRelate) #working fine! and this is what we need for the pruning steps anyhow
BiocManager::install("LEA")
require(LEA)
require(dplyr)
BiocManager::install("calibrate")
require(calibrate)
require(ggplot2)
#BiocManager::install("gcc")
require(gcc)
require(vegan)
require(ggfortify)

closeAllConnections()
snpgdsClose(test1.gds) 
```
### Step1: Data Wrangling
```{r}
# Load our allsamples.vcf.gz
vcf.fn <- ("/project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz")
```

### Step2: Plot some stats about our VCF
```{r}
install.packages(c("vcfR"))
library(vcfR) #non-zero exit status
vcf <- read.vcfR( vcf.fn, verbose = FALSE )
#Then plot important statistics summed over entire VCF
chrom <- create.chromR(name='Whole Genome Sequencing', vcf=vcf) 
plot(chrom) #plot the data

#Then extract the allele depths per each sample (DP field of VCF) and plot distribution of allele depths of all sites per each sample. NB: You may inspect and visualize other fields of VCF, e.g. allele depth (AD) or genotype quality (GQ)
# read depth distribution per individual 
dp<- extract.gt(vcf, element = 'DP', as.numeric=TRUE)
#pdf("DP_RAD_data.pdf", width = 10, height=3) # boxplot #where did this go?
#par(mar=c(8,4,1,1)) 

#Plot Read Depth (DP) per individual 
boxplot(dp, las=3, col=c("#BOBCAT1", "#BOBCAT2"), ylab="Read Depth (DP)",las=2, cex=0.4, cex.axis=0.5)
#zoom to smaller values 
boxplot(dp, las=3, col=c("#C0C0C0", "#808080"), ylab="Read Depth (DP)",
las=2, cex=0.4, cex.axis=0.5, ylim=c(0,50))
abline(h=8, col="red") 
```

### Convert our vcf to gds format and select only biallelic SNPs.This will only select biallelic SNPs, which is another pruning method because we have a lot of multiallelic SNPs in our dataset
```{r}
snpgdsVCF2GDS(vcf.fn, "test2.gds", method="biallelic.only") 
    #Number of samples: 61 (ALL lynx samples, no bobcat, no hybrid...import 1823879 variants.

# Review our "test1.gds" file
snpgdsSummary("test2.gds")
  #The file name: /home/tl50a/test1.gds 
  #The total number of samples: 61 
  #The total number of SNPs: 1823879 
  #SNP genotypes are stored in SNP-major mode (Sample X SNP)

# Open the .gds file
genofile <- snpgdsOpen("/home/tl50a/test2.gds", readonly= FALSE, allow.duplicate = TRUE)
# Assign sample.id names
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
 #[1] "L155.variant4.variant2"         "LIC11.variant4.variant2"...
```
We've already filtered a lot of our SNPs by selecting biallelic SNPs only. We will reduce our SNPset further by removing monoSNPs and LDpruning. 
Questions: 1. What is meant by monosnp (only one individual has SNP?). A monomorphic site is one site in which all the individuals have the same form (genotype) 2. What is meant by ld.threshold? 
### Proceed to snpgdsLDpruning
```{r}
#Try different LD thresholds for sensitivity analysis. #Need to include autosome.only = FALSE clause to run.
set.seed(1) #should be 1000, reduced for trial runs
snpset.remove.monosnp <- snpgdsLDpruning(genofile, autosome.only = FALSE, remove.monosnp = TRUE, ld.threshold=0.1, verbose = TRUE)
    #snpset.remove.monosnp=TRUE 53,421 markers are selected in total/1,823,879 with ld=0.1
#consider different values for missing.rate and maf? => same biological result in structure etc?

# Create an object to hold our ~53,421 selected SNPs
snpset.id <- unlist(snpset.remove.monosnp)

# 1 Take our selected snps from snpset.id and generate a .ped file
snpgdsGDS2PED(genofile, "ped1.fn", sample.id=sample.id, snp.id=snpset.id, use.snp.rsid=TRUE, format=c("A/G/C/T", "A/B", "1/2"), verbose=TRUE)

# 2 Convert the ped file to geno file for snmf
ped2geno("/home/tl50a/ped1.fn.ped", output.file = ("/home/tl50a/ped1.fn.ped.geno"), force = TRUE)

	#- number of detected individuals:	61
	#- number of detected loci:		53265 #previously only 28,977. Interesting that lynx only gen. more SNPs than lynx + bobcat calls?
  #"/home/tl50a/ped1.fn.ped.geno" #this is now a genofile which we will use

# 3 save our filtered/reduced SNPset as a vcf
# We were using a command [seqGDS2VCF("ped1.fn.ped.geno", "output.vcf")] from the SeqArray package on BiocManager but that is not working. This is the solution:
# We were able to get this converted to a VCF using PLINK on the cluster. See my HackMD for instructions on how to do that (it's very simple. https://hackmd.io/o4-w-bDgSv2HRn6agxd2CA)
```

#We now have a genofile with filtered SNPs that we will use for analyses
#Let's start running our population analyses. 
#Try running sNMF analysis with our new .geno file
```{r}
#repetitions should be 1000, but we have 1 here for Rmarkdown to run efficiently
obj.snmf = snmf("/home/tl50a/ped1.fn.ped.geno", K = 1:10, project = "new", repetitions = 1, tolerance = 0.00001, entropy=TRUE, ploidy = 2)

obj.snmf = load.snmfProject("ped1.fn.ped.snmfProject")
```
#Examine the results

Number of iterations: 134

Least-square error: 220638.793373
Write individual ancestry coefficient file /home/tl50a/ped1.fn.ped.snmf/K10/run1/ped1.fn.ped_r1.10.Q:		OK.
Write ancestral allele frequency coefficient file /home/tl50a/ped1.fn.ped.snmf/K10/run1/ped1.fn.ped_r1.10.G:	OK.

[1] "*************************************"
[1] "*    cross-entropy estimation       *"
[1] "*************************************"
summary of the options:

        -n (number of individuals)         61
        -L (number of loci)                53265
        -K (number of ancestral pops)      10
        -x (genotype file)                 /home/tl50a/ped1.fn.ped.geno
        -q (individual admixture)          /home/tl50a/ped1.fn.ped.snmf/K10/run1/ped1.fn.ped_r1.10.Q
        -g (ancestral frequencies)         /home/tl50a/ped1.fn.ped.snmf/K10/run1/ped1.fn.ped_r1.10.G
        -i (with masked genotypes)         /home/tl50a/ped1.fn.ped.snmf/masked/ped1.fn.ped_I.geno
        - diploid

Cross-Entropy (all data):	 0.122603
Cross-Entropy (masked data):	 0.307019
The project is saved into :
 ped1.fn.ped.snmfProject 

To load the project, use:
 project = load.snmfProject("ped1.fn.ped.snmfProject")

To remove the project, use:
 remove.snmfProject("ped1.fn.ped.snmfProject")

```{r echo= T, include = T}
# plot cross-entropy criterion of all runs of the project
#This plot helps us identify the number of clusters to identify as the most likely for our population (looks like k=3)
plot(obj.snmf, cex = 1.2, col = "lightblue", pch = 19)

# get the cross-entropy value for each run at K = ?
ce <- cross.entropy(obj.snmf, K = 3) #we only did 1 run for R Markdown. More realistically we would perform many runs and select the best one

# select the run with the lowest cross-entropy value for K = ?
best <- which.min(ce) # the best run is #1 at k=3

#At the end of the run, the qmatrix object contains the matrix of ancestry coefficients for each individual and for K = 3 clusters. The Q-matrix has 24 rows and 3 columns, and it is traditionally displayed using a barplot representation. For this representation, we just use the barplot function of R (Figure 1).
qmatrix = Q(obj.snmf, K = 3, run = best)
qmatrix

# Name the cluster assignment for each individual
cluster<- apply(qmatrix, 1, which.max) #this corresponds with the 1:24 order 
[1] 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 1 2 2 3 3 1 3 1 1 3 1 1
[42] 1 1 3 2 1 3 1 1 3 3 3 3 1 1 1 1 1 3 1 3 #same
> 
# e.g. Individual #29 (BOBCAT) was assigned to cluster 3

#Ancestry proportions
barchart(obj.snmf, K = 3, run = best,  border = NA, space = .2, col = c("orange","violet","lightgreen"), xlab = "Individuals", ylab = "Ancestry proportions",  main = "Ancestry matrix", horiz = FALSE,   names.arg = c("a475","a494","b276","a202","a33","a772","a803","a818","b114","b124","b13","b188","b554","c165","c323","f264","f457","fha_024","fha_042","fha_043","l09_007","L155","LIC11","LIC20","LIC23","LIC24","LIC27B","LIC28","LIC31","LIC32","LIC36","LIC46","LIC47","LIC48","LIC54","LIC57","LIC60","LIC8","LIC9","LIT2","LIT5","LRK10","LRK11","LRK12","LRK13","LRK17","LRK22","a109","a182","a507","a697","a794","a857","b23","b90","c548","cb15","cb42","cb7","l09_003","l09_015"), cex.names=0.65, las = 2) -> bp 

#sample_id order
sample.id
"L155", "LIC11", "LIC20", "LIC23", "LIC24", "LIC27B", "LIC28", "LIC31", "LIC32", "LIC36", "LIC46", "LIC47", "LIC48", "LIC54", "LIC57", "LIC60", "LIC8", "LIC9", "LIT2", "LIT5", "LRK10", "LRK11", "LRK12", "LRK13", "LRK17", "LRK22", "a109", "a182", "a202", "a33", "a475", "a494", "a507", "a697", "a772", "a794", "a803", "a818", "a857", "b114", "b124", "b13", "b188", "b23", "b276","b554", "b90", "c165", "c323", "c548", "cb15", "cn42", "cb7", "f264", "f457", "fha_024", "fha_042", "fha_043", "l09_003", "l09_007", "l09_015"

#In order to correctly label the ancestry proportions for each individual, we need to extract the order of the sample id's from the obj.snmf. 
ids<- bp$order
#1] 31 32 45 29 30 35 37 38 40 41 42 43 46 48 49 54 55 56 57 58 60  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 33 34 36 39 44 47 50 [57] 51 52 53 59 61

#And re-order them in the correct order for the Admixture plot
"a475","a494","b276","a202","a33","a772","a803","a818","b114","b124","b13","b188","b554","c165","c323","f264","f457","fha_024","fha_042","fha_043","l09_007","L155","LIC11","LIC20","LIC23","LIC24","LIC27B","LIC28","LIC31","LIC32","LIC36","LIC46","LIC47","LIC48","LIC54","LIC57","LIC60","LIC8","LIC9","LIT2","LIT5","LRK10","LRK11","LRK12","LRK13","LRK17","LRK22","a109","a182","a507","a697","a794","a857","b23","b90","c548","cb15","cb42","cb7","l09_003","l09_015"

```

### Ancestry proportions
```{r echo = F, include = FALSE}
barchart(obj.snmf, K = 3, run = best,  border = NA, space = .2, col = c("orange","violet","lightgreen"), xlab = "Individuals", ylab = "Admixture coefficients", horiz = FALSE,   names.arg = c("a475","a494","b276","a202","a33","a772","a803","a818","b114","b124","b13","b188","b554","c165","c323","f264","f457","fha_024","fha_042","fha_043","l09_007","L155","LIC11","LIC20","LIC23","LIC24","LIC27B","LIC28","LIC31","LIC32","LIC36","LIC46","LIC47","LIC48","LIC54","LIC57","LIC60","LIC8","LIC9","LIT2","LIT5","LRK10","LRK11","LRK12","LRK13","LRK17","LRK22","a109","a182","a507","a697","a794","a857","b23","b90","c548","cb15","cb42","cb7","l09_003","l09_015"),cex.names=0.4, las = 2) -> bp 

```

### Estimating Individual Admixture Proportions from NGS data
This was our most successful plot. names.arg needs to be in the order of sample.id (not the bp$order)
```{r echo=T, include = T}
barplot(t(qmatrix),col=c("orange","violet","lightgreen"),border = NA, space = .2,xlab = "Individuals", ylab = "Admixture coefficients", main = "Ancestry matrix", horiz = FALSE, names.arg = c("L155", "LIC11", "LIC20", "LIC23", "LIC24", "LIC27B", "LIC28", "LIC31", "LIC32", "LIC36", "LIC46", "LIC47", "LIC48", "LIC54", "LIC57", "LIC60", "LIC8", "LIC9", "LIT2", "LIT5", "LRK10", "LRK11", "LRK12", "LRK13", "LRK17", "LRK22", "a109", "a182", "a202", "a33", "a475", "a494", "a507", "a697", "a772", "a794", "a803", "a818", "a857", "b114", "b124", "b13", "b188", "b23", "b276","b554", "b90", "c165", "c323", "c548", "cb15", "cn42", "cb7", "f264", "f457", "fha_024", "fha_042", "fha_043", "l09_003", "l09_007", "l09_015"),cex.names=0.5, las = 2)
#legend v1, n= 19;  v2, n= 3;  v3, n=2
```

### Other post-run sNMF treatments
Call the ancestral genotype frequency matrix, G, for the best run of k=3. 
```{r echo=T, include =T}
gmatrix = G(obj.snmf, K = 3, run = best)
head(gmatrix)
barplot(gmatrix,  border = NA, space = 0.2,  col = c("orange","violet","lightgreen"),  xlab = "Populations", ylab = "Ancestry proportions",  main = "Ancestry matrix") ->gp
```

### Fst  Estimation
Given two or more populations, Fst can be estimated by the method of Weir & Cockerham (1984).
```{r}
#get sample id information
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

figurelabels<- c("L155", "LIC11", "LIC20", "LIC23", "LIC24", "LIC27B", "LIC28", "LIC31", "LIC32", "LIC36", "LIC46", "LIC47", "LIC48", "LIC54", "LIC57", "LIC60", "LIC8", "LIC9", "LIT2", "LIT5", "LRK10", "LRK11", "LRK12", "LRK13", "LRK17", "LRK22", "a109", "a182", "a202", "a33", "a475", "a494", "a507", "a697", "a772", "a794", "a803", "a818", "a857", "b114", "b124", "b13", "b188", "b23", "b276","b554", "b90", "c165", "c323", "c548", "cb15", "cn42", "cb7", "f264", "f457", "fha_024", "fha_042", "fha_043", "l09_003", "l09_007", "l09_015") 

#Add population information to your gds genofile 
add.gdsn(genofile, "pop_code", c( 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 2, 2, 3, 3, 1, 3, 1, 1, 3, 1, 1, 1, 1, 3, 2, 1, 3, 1, 1, 3, 3, 3, 3, 1, 1, 1, 1, 1, 3, 1, 3)) #all identified to 1 population for starters

pop_code <- read.gdsn(index.gdsn(genofile, "pop_code"))

pop_code_factor<- as.factor(pop_code) #Fst needs the population information to be a factor, not a number

v <- snpgdsFst(genofile, population=pop_code_factor, method=c("W&C84"), sample.id=sample.id, remove.monosnp = FALSE,
    snp.id=snpset.id, autosome.only=FALSE, maf=0.05,
    missing.rate=NaN, with.id=TRUE, verbose=TRUE)
v$Fst
v$MeanFst
summary(v$FstSNP)
```
#Fst estimation on genotypes:
Excluding 47,370 SNPs (monomorphic: TRUE, MAF: 0.05, missing rate: NaN)
Working space: 61 samples, 5,993 SNPs
Method: Weir & Cockerham, 1984
#of Populations: 3
    1 (18), 2 (3), 3 (40)
$Fst
[1] 0.134448

$MeanFst
[1] 0.1184698


## Principal Component Analysis (PCA)  
The functions in SNPRelate for PCA include calculating the genetic covariance matrix from genotypes, computing the correlation coefficients between sample loadings and genotypes for each SNP, calculating SNP eigenvectors (loadings), and estimating the sample loadings of a new dataset from specified SNP eigenvectors.  
```{r PCA, echo=TRUE, warning=FALSE}
# Run PCA using the ldpruned set of 30k SNPs (snpset.id) instead of the full set. 
pca <- snpgdsPCA(genofile, sample.id=sample.id, snp.id=snpset.id, autosome.only = FALSE, remove.monosnp = TRUE, verbose= TRUE)
  #Excluding 0 SNP (monomorphic: TRUE, MAF: NaN, missing rate: NaN)
  #Working space: 65 samples, 28,977 SNPs

#The following code shows how to calculate the percent of variation that is accounted for by the top principal components. It is clear to see the first two eigenvectors hold the largest percentage of variance among the population, although the total variance accounted for is still less the one-quarter of the total.

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#get population information
id_code<- read.gdsn(index.gdsn(genofile, "sample.id"))

head(cbind(sample.id, pop_code))

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(pop_code)[match(pca$sample.id,sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
tab #view the dataframe

plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1") #three distinct groups

#Draw
plot.new()
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
textxy(tab$EV2, tab$EV1, labs=figurelabels2[1:61], cex = 0.6, pos = 4)
#legend("bottomleft", legend=c("Population 1","Population 2", "Population 3"), cex = 1,pch="o", col=c("black","red","green"))


figurelabels2<- c("", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "a475", "a494", "", "", "", "", "", "", "", "", "", "", "", "", "b276","", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "") 

#Plot the principal component pairs for the first four PCs:
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)

#Parallel coordinates plot for the top principal components:
datpop <- factor(pop_code)[match(pca$sample.id, sample.id)]
parcoord(pca$eigenvect[,1:16], col=datpop)

#To calculate the SNP correlations between eigenvactors and SNP genotypes:
# Get chromosome index
chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)
savepar <- par(mfrow=c(2,1), mai=c(0.45, 0.55, 0.1, 0.25))
for (i in 1:2)
{
    plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="", ylab=paste("PC", i), pch="+")
} #col=chr
par(savepar)
```
## IBD analysis
In addition to these summary statistics, we also want to filter on relatedness criteria. We use the SNPRelate package to perform identity-by-descent (IBD) analysis. This package requires that the data be transformed into a GDS format file. IBD analysis is performed on only a subset of SNPs that are in linkage equilibrium by iteratively removing adjacent SNPs that exceed an LD threshold in a sliding window using the snpgdsLDpruning function.

ld.thresh <- 0.1    # LD cut-off
kin.thresh <- 0.1   # Kinship cut-off

### Find IBD coefficients using Method of Moments procedure.  Include pairwise kinship.
```{r echo=TRUE, warning=FALSE}
# Estimate IBD coefficients
ibd <- snpgdsIBDMoM(genofile, sample.id=sample.id, snp.id=snpset.id, autosome.only = FALSE, kinship = TRUE, num.thread=2, verbose= TRUE)
#maf (minor allele frequency) threshold is the frequency at which the second most common allele occurs in a given population. 

# Make a dataframe
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)
ibd.coeff 

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
    xlab="k0", ylab="k1") #main="Lynx Bobcat and Hybrid Sample Relatedness (MoM)"
lines(c(0,1), c(1,0), col="red", lty=2)
#textxy(ibd.coeff$k0, ibd.coeff$k1, labs = figurelabels2[1:65], cex = 0.45, pos = 3)

# Check if there are any candidates for relatedness
#Set thresholds 
ld.thresh <- 0.1    
kin.thresh <- 0.2599999 #what should the threshold be?

# Check if there are any candidates for relatedness
ibdcoeff <- ibd.coeff[ ibd.coeff$kinship >= kin.thresh, ]

# iteratively remove samples with high kinship starting with the sample with the most pairings
#related.samples <- NULL while ( nrow(ibdcoeff) > 0 ) {

    # count the number of occurrences of each and take the top one
    sample.counts <- arrange(count(c(ibdcoeff$ID1, ibdcoeff$ID2)), -freq)
    rm.sample <- sample.counts[1, 'x']
    cat("Removing sample", as.character(rm.sample), 'too closely related to', sample.counts[1, 'freq'],'other samples.\n')

    # remove from ibdcoeff and add to list
    ibdcoeff <- ibdcoeff[ibdcoeff$ID1 != rm.sample & ibdcoeff$ID2 != rm.sample,]
    related.samples <- c(as.character(rm.sample), related.samples)
}

# filter genotype and clinical to include only unrelated samples
genotype <- genotype[ !(rownames(genotype) %in% related.samples), ]
clinical <- clinical[ !(clinical$FamID %in% related.samples), ]

geno.sample.ids <- rownames(genotype)

cat(length(related.samples), "similar samples removed due to correlation coefficient >=", kin.thresh,"\n") 

    ```
    
Estimating IBD: Maximum Likelihood Estimation  
#Identical by descent (IBD) describes a matching segment of DNA shared by two or more individuals that has been inherited from a common ancestor without any intervening recombination. Alleles that are IBD have been interested by a recent common ancestor (and indicate recent gene flow).
```{r message=TRUE, warning=FALSE}
# Estimate IBD coefficients
set.seed(1) #should be 1000
ibd <- snpgdsIBDMLE(genofile, sample.id = NULL, snp.id=snpset.id,autosome.only = FALSE, remove.monosnp= TRUE, maf=0.05, missing.rate=0.05, kinship = TRUE, num.thread=2, verbose = TRUE)

# Make a dataframe
ibd.coeff <- snpgdsIBDSelection(ibd)

plot.new()
plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1), col=c("black"),
    xlab="k0", ylab="k1") #as.factor(ibd.coeff$ID2) #main="Canada Lynx Relatedness (MLE)"
lines(c(0,1), c(1,0), col="red", lty=2)
#textxy(ibd.coeff$k0, ibd.coeff$k1, labs = figurelabels[1:65], cex = 0.45, pos = 3)
#legend("left", legend=as.factor(ibd.coeff$ID2), cex = 0.4,pch="o", col=as.factor(ibd.coeff$ID2)) #legend optional
```

## Identity-By-State Analysis
For n study individuals, snpgdsIBS() can be used to create a n×n matrix of genome-wide average IBS pairwise identities:
```{r IBS, echo=TRUE, warning=FALSE}
ibs <- snpgdsIBS(genofile, sample.id=sample.id, snp.id=snpset.id, num.thread=2,autosome.only = FALSE)
pop.idx <- order(pop_code)
#The heat map is shown:
# individulas in the same population are clustered together
plot.new()
image(ibs$ibs[pop.idx, pop.idx], col=terrain.colors(16))
textxy(ibs$ibs, ibs$pop.idx, labs = ibs$snp.id, cex = 0.45, pos = 3)

#To perform multidimensional scaling analysis on the n×n matrix of genome-wide IBS pairwise distances:
loc <- cmdscale(1 - ibs$ibs, k = 2)
x <- loc[, 1]; y <- loc[, 2]
race <- as.factor(pop_code) #All showing =1 because we set them all as a single population

plot(x, y, col=race, xlab = "", ylab = "") #main = "Multidimensional Scaling Analysis (IBS)"
legend("topleft", legend=levels(race), pch="o", text.col=1:nlevels(race))

#To perform cluster analysis on the n×n matrix of genome-wide IBS pairwise distances, and determine the groups by a permutation score:

set.seed(1)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile, num.thread=2, autosome.only = FALSE))
# Determine groups of individuals automatically
rv <- snpgdsCutTree(ibs.hc)
table(rv$samp.group)
#Here is the population information we have known:
# Determine groups of individuals by population information
rv2 <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code), label.H=TRUE, label.Z=TRUE)
plot(rv2$dendrogram) #, main="HapMap Phase II"
#textxy(ibd.coeff$k0, ibd.coeff$k1, labs = figurelabels[1:65], cex = 0.45, pos = 3)
#legend("topright", legend=levels(race), col=1:nlevels(race), pch=19, ncol=4)
```
## Sample Level Filtering
We need to filter our samples. The second state of data pre-processing is removing individuals due to missing data, sample contamination, correlation (for population genomics) and gender ambiguity or discordance. In our study, we address these issues by filtering on call rate, heterozygosity, cryptic relatedness and duplicates using identity-by-descent, and we visually assess ancestry

Getting started
Sample level quality control for missing data and heterozygosity is achieved using the row.summary function from snpStats.
BiocManager::install("snpStats")
library(snpStats)
library(plyr)

geno<- read.pedfile("/home/tl50a/ped.fn.ped")
## Obtain the SnpMatrix object (genotypes) table from geno list
#Note: Phenotypes and covariates will be read from the clinical data file, below
genotype <- geno$genotype
print(genotype)
## Create SNP summary statistics (MAF, call rate, etc.)
snpsum.col <- col.summary(genotype)
print(head(snpsum.col))

Using these summary statistics, we keep the subset of SNPs that meet our criteria for minimum call rate and minor allele frequency.
#Setting thresholds
call <- 0.95
minor <- 0.01

#Filter on MAF and call rate
use <- with(snpsum.col, (!is.na(MAF) & MAF > minor) & Call.rate >= call)
use[is.na(use)] <- FALSE                # Remove NA's as well

cat(ncol(genotype)-sum(use),"SNPs will be removed due to low MAF or call rate.\n") 
#28,687 SNPs will be removed due to low MAF or call rate.

## Subset genotype and SNP summary data for SNPs that pass call rate and MAF criteria
genotype <- genotype[,use]
snpsum.col <- snpsum.col[use,]

print(genotype)
A SnpMatrix with  69 rows and  1828 columns # 1828 SNPs remain
Row names:  BOBCAT1.variant2.variant2 ... l09_015.variant2.variant.variant 
Col names:  locus.8 ... locus.30515 
```

# Create sample statistics (Call rate, Heterozygosity)
snpsum.row <- row.summary(genotype)

# Add the F stat (inbreeding coefficient) to snpsum.row
MAF <- snpsum.col$MAF
callmatrix <- !is.na(genotype)
hetExp <- callmatrix %*% (2*MAF*(1-MAF))
hetObs <- with(snpsum.row, Heterozygosity*(ncol(genotype))*Call.rate)
snpsum.row$hetF <- 1-(hetObs/hetExp)

head(snpsum.row)
Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
```


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `R` `SNPRelate` `sNMF` `LDpruning`
