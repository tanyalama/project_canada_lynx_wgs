#! /bin/bash

# This script is for the autosomes (chr 1 to 66)
individuals_array=(NC_044317.1 NC_044318.1 NC_044319.1 NC_044320.1 NC_044321.1)

for i in "${individuals_array[@]}"
do

plink --vcf  /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/mainelynx/6SV_mainelynx.vcf.gz --r2 square 'yes-really' --out /project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_LD/mainelynx/mainelynx_${i} --chr ${i} --make-bed --allow-extra-chr

done
