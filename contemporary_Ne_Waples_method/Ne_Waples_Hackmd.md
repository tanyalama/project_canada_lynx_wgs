---
title: 'Estimation of contemporary effective population size (Ne) in Canada lynx'
disqus: hackmd
---

Estimation of contemporary effective population size (Ne) in Canada lynx
===
Author: Tanya Lama tanya.m.lama@gmail.com

## Objective
Relationships between LD, Ne, and Nc

## Background reading (Nielsen and Slatkin)
Wright-Fischer Model: p. 22-27
Effective population size: p. 43-46
Linkage Disequilibrium: p. 108-112, including boxes 6.1-6.3

## Rough estimates of the census population sizes (Nc).
(from MDIFW) I think the range is 300-1500

## Table of Contents
[TOC]
module load plink/1.90b6.9

## Start an interactive session: 
```
bsub -q interactive -R rusage[mem=24000] -n 1 -W 300 -Is bash
```

## Calculate LD with Plink with a wrapper script
We've written a wrapper script to calculate LD chr by chr for each populations (e.g. /southSLR). This script uses plink to calculate r2 LD and generate a bed file required for our next step, calculating Ne using an Rscript. 

Usage is: 
```

```

**This is a full list of chromosomes**
NC_044303.1 NC_044304.1 NC_044305.1 NC_044306.1 NC_044307.1 NC_044308.1 NC_044309.1 NC_044310.1 NC_044311.1 NC_044312.1 NC_044313.1 NC_044314.1 NC_044315.1 NC_044316.1 NC_044317.1 NC_044318.1 NC_044319.1 NC_044320.1 NC_044321.1 NW_022059692.1 NW_022059693.1 NW_022059694.1 NW_022059695.1 NW_022059696.1 NW_022059697.1 NW_022059698.1 NW_022059699.1 NW_022059700.1 NW_022059701.1 NW_022059702.1 NW_022059703.1 NW_022059704.1 NW_022059705.1 NW_022059706.1 NW_022059707.1 NW_022059708.1 NW_022059709.1 NW_022059710.1 NW_022059711.1 NW_022059712.1 NW_022059713.1 NW_022059714.1 NW_022059715.1 NW_022059716.1 NW_022059717.1 NW_022059718.1 NW_022059719.1 NW_022059720.1 NW_022059721.1 NW_022059722.1 NW_022059723.1 NW_022059724.1 NW_022059725.1 NW_022059726.1 NW_022059727.1 NW_022059728.1 NW_022059729.1 NW_022059730.1 NW_022059731.1 NW_022059732.1 NW_022059733.1 NW_022059734.1 NW_022059735.1 NW_022059736.1 NW_022059737.1 NW_022059738.1


```
bsub -q short -R rusage[mem=130000] -n 1 -W 4:00 -R span\[hosts=1\] "plink --vcf  /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/southSLR/6SV_southSLR.vcf.gz --r2 square 'yes-really' --out /project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_LD/southSLR/southSLR_allchr --allow-extra-chr"
```

2. try just one chr? --chr NC_044303.1 #8272221
bsub -q short -R rusage[mem=130000] -n 1 -W 4:00 -R span\[hosts=1\] "plink --vcf  /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/southSLR/6SV_southSLR.vcf.gz --r2 square 'yes-really' --out ./work/southSLR_NC0443031 --chr NC_044303.1 --allow-extra-chr" #this is running pretty swiftly

bsub -q short -R rusage[mem=240000] -n 1 -W 4:00 -R span\[hosts=1\] "plink --vcf  /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/NFLD/6SV_NFLD.vcf.gz --r2 square 'yes-really' --out ./work/NFLD_NC0443041 --chr NC_044304.1 --allow-extra-chr" #this is running pretty swiftly 8273569

southSLR_NC0443031 --chr NC_044303.1 is done!
We can now move forward with the R script!
This worked but is now taking forever to run on R. Let's try with a smaller chromosome

## Make .bed .bim .fam
In addition to .ld, we'll use plink to make these additional file types
plink --vcf  /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/southSLR/6SV_southSLR.vcf.gz --make-bed --out ./work/southSLR_NC0443031 --chr NC_044303.1 --allow-extra-chr

plink --vcf  /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/NFLD/6SV_NFLD.vcf.gz --make-bed --out ./work/NFLD_NC0443041 --chr NC_044304.1 --allow-extra-chr

## Estimate Ne in R
Working directory is /project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/scripts

We created an R script that sources a more complex Ne calculation from source('./scripts/R_functions.r')

1. Load dependencies
module load R/3.6.1
module load gcc/8.1.0

2. Usage is:
```
bsub -q long -R rusage[mem=60000] -n 4 -W 24:00 -R span\[hosts=1\] "Rscript ./NFLD_NC0443041_estimate_Ne.r" #8279516
#actually used 80000.00 MB so now running on 120

bsub -q long -R rusage[mem=160000] -n 1 -W 24:00 -R span\[hosts=1\] "Rscript ./southSLR_NC0443031_estimate_Ne.r" #8279515

#We're having a hard time getting any of our Ne estimates to run and have had the R script exit up to 160G


```
3. Output with .Ne suffix is located in /project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work

## Results


| nFLD | southSLR | northSLR |
| -------- | -------- | -------- |
| n     | n     | n     |


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `Ne``Effective Population Size`
