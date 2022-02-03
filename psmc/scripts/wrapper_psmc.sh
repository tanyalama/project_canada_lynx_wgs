#!/bin/bash

# Uncomment the line below and fill in the name of the individuals.
individuals_array=(LIC60 LIC8 LIC9 LIT2 LIT5 LRK10 LRK11 LRK12 LRK13 LRK17 LRK22)

for i in "${individuals_array[@]}"
do

psmc -p "4+25*2+4+6" -o ${i}.psmc ${i}.psmcfa
done

#bsub -q long -W 48:00 -R rusage[mem=8000] -n 1 "./wrapper_psmc.sh"
#2502 18h for 9 samples
