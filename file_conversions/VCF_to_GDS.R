#the goal here is to take a filtered VCF and convert it to genofile, GDS and ped/map formats usable by various genomics packages in R (LEA, SNPRelate, vcfR)
#Prepare our workspace and load required packages
setwd("/home/tl50a")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

BiocManager::install("ape")
require(ape)
BiocManager::install("vcfR") #non zero exit status if you are using 3.6.0
require(vcfR) #not working with 3.6.0
require(SNPRelate) #This is what we need for the LD pruning step anyhow
BiocManager::install("LEA") #need this for sNMF structure-like analyses
require(LEA)
require(dplyr)
BiocManager::install("calibrate")
require(calibrate)
require(ggplot2) #for plotting things
#BiocManager::install("gcc")
require(gcc)
require(vegan)
require(ggfortify)

closeAllConnections()#if you're returning to a GDS object, you'll need to open/close them every session or you'll have errors below
snpgdsClose(test2.gds) 
```
#Step1: Data Wrangling
```{r}
# Load our VCF
vcf.fn <- ("/project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6_SV_rmClusterSNP_BiSNP_SV_HardFilter_SV_all_mLynCan4_v1.p_HighQualSites_processed.vcf.gz")
```

#Step2: Plot some stats about our VCF. If vcfR did not load properly for you above, this will not work. You'll have to update R to 4.0.0
```{r}
install.packages(c("vcfR"))
library(vcfR) #non-zero exit status, so we weren't able to run this on the cluster. Had to do it locally on a subset of the data
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

# Convert our vcf to GDS format and select only biallelic SNPs. This will only select biallelic SNPs, which is another pruning method because we have a lot of multiallelic SNPs in our dataset
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

We've already filtered a lot of our SNPs out by selecting biallelic SNPs only. We will reduce our SNPset further by removing monoSNPs and LDpruning. 
Questions: 1. What is meant by monosnp. A monomorphic site is one site in which all the individuals have the same form (genotype) 

# Proceed to snpgdsLDpruning
```{r}
#Try different LD thresholds for sensitivity analysis. #Need to include autosome.only = FALSE clause to run.
set.seed(1) #should be 1000, reduced for trial runs
snpset.remove.monosnp <- snpgdsLDpruning(genofile, autosome.only = FALSE, remove.monosnp = TRUE, ld.threshold=0.1, verbose = TRUE)
    #snpset.remove.monosnp=TRUE 53,421 markers are selected in total/1,823,879 with ld=0.1
#consider different values for missing.rate and maf? => same biological result in structure etc?

# Create an object to hold our ~53,421 selected SNPs
snpset.id <- unlist(snpset.remove.monosnp)

# 1 Take our snpset.id snpset and generate a .ped file
snpgdsGDS2PED(genofile, "ped1.fn", sample.id=sample.id, snp.id=snpset.id, use.snp.rsid=TRUE, format=c("A/G/C/T", "A/B", "1/2"), verbose=TRUE)

# 2 Convert the ped file to geno file for snmf. The genofile is the input for much of the code below. 
ped2geno("/home/tl50a/ped1.fn.ped", output.file = ("/home/tl50a/ped1.fn.ped.geno"), force = TRUE)

	#- number of detected individuals:	61
	#- number of detected loci:		53265 #previously only 28,977. Interesting that lynx only gen. more SNPs than lynx + bobcat calls?
  #"/home/tl50a/ped1.fn.ped.geno" #this is now a genofile which we will use

# 3 You might want to save our filtered/reduced SNPset as a vcf for other analyses. 
# We were using a command [seqGDS2VCF("ped1.fn.ped.geno", "output.vcf")] from the SeqArray package on BiocManager but that is not working on R 3.6.0. This is the solution:
# We were able to get this converted to a VCF using PLINK on the cluster. See my HackMD for instructions on how to do that (it's very simple. https://hackmd.io/o4-w-bDgSv2HRn6agxd2CA)
```

#We now have a genofile with filtered SNPs that we will use for analyses in R.
#Let's start running our population analyses in the project_canada_lynx_wgs/SNPRelate folder
