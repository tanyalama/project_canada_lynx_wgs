#!/bin/bash
#Uncomment the line below and fill in the name of the individuals.
individuals_array=(NC_044303.1
NC_044304.1
NC_044305.1
NC_044306.1
NC_044307.1
NC_044308.1
NC_044309.1
NC_044310.1
NC_044311.1
NC_044312.1
NC_044313.1
NC_044314.1
NC_044315.1
NC_044316.1
NC_044317.1
NC_044318.1
NC_044319.1
NC_044320.1
NC_044321.1
NW_022059692.1
NW_022059693.1
NW_022059695.1
NW_022059696.1
NW_022059697.1
NW_022059699.1
NW_022059700.1
NW_022059701.1
NW_022059706.1
NW_022059707.1
NW_022059710.1
NW_022059714.1
NW_022059715.1
NW_022059716.1
NW_022059717.1
NW_022059718.1
NW_022059722.1
NW_022059723.1
NW_022059724.1
NW_022059725.1
NW_022059726.1
NW_022059727.1
NW_022059728.1
NW_022059729.1
NW_022059730.1
NW_022059731.1
NW_022059733.1
NW_022059734.1
NW_022059735.1
NW_022059736.1
NW_022059738.1)

for i in "${individuals_array[@]}"
do
vcftools --gzvcf snpset_LD0.2_MAF0.05_missing_0.05.vcf.gz  --chr ${i} --recode  --out ${i}

done
