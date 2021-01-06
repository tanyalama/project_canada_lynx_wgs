#!/bin/bash

# Uncomment the line below and fill in the name of the folders to check hashes
individuals_array=(6SV_LDpruned_biallelic_maine_lynx  6SV_unfiltered_maine_lynx 6SV_unfiltered_NFLDnorth_lynx 6SV_unfiltered_northSLR_lynx 6SV_unfiltered_southSLR_lynx 6SV_unfiltered_lynx 6SV_unfiltered_NFLD_lynx 6SV_unfiltered_NFLDsouth_lynx 6SV_unfiltered_northsouth_lynx)

for i in "${individuals_array[@]}"
do

md5sum ./${i}/* > ./${i}/hashes

done

#bsub -q short -W 1:00 -R rusage[mem=1000] -n 2 -R span\[hosts=1\] "./md5sum_folders.sh"

md5sum_folders.sh (END) 
