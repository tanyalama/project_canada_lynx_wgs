#! /bin/bash

# This script is for the autosomes (chr 1 to 66)
individuals_array=(NC_044304.1 NC_044305.1 NC_044306.1 NC_044307.1 NC_044308.1 NC_044309.1 NC_044310.1 NC_044311.1 NC_044312.1 NC_044313.1 NC_044314.1 NC_044315.1 NC_044316.1 NC_044317.1 NC_044318.1 NC_044319.1 NC_044320.1 NC_044321.1 NW_022059692.1 )

for i in "${individuals_array[@]}"
do

singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/mainelynx/6SV_mainelynx.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/mainelynx/6SV_mainelynx.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/mainelynx/mainelynx_${i} -chr ${i} -samplesize 26 -itemsTH n -threads 8

done
