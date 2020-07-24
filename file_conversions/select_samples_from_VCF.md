---
title: 'Subsetting a VCF'
disqus: hackmd
---

Subsetting a VCF
===

We've resolved from much trial and error that the best way of subsetting a VCF is using **vcftools**

This script can be used to drop/select individuals from a given VCF into subsets as you'd like

Each population-level VCF and associated ped/map file will be located in the Rgenomics/mLynCan4_v1.p folder along with the "mother" 6SV VCF


## Table of Contents

[TOC]

**Our output files will include:** 
6SV_southSLR.vcf.gz .map/.ped
6SV_northSLR.vcf.gz .map/.ped
6SV_NFLD.vcf.gz .map/.ped
6SV_bobcat.vcf.gz .map/.ped

northsouth.vcf.gz
northNFLD.vcf.gz
southNFLD.vcf.gz

module load plink/1.90b6.9



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
## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `Templates` `Documentation`
