#! /bin/bash

samtools="/share/pkg/samtools/1.2" # Please give the path to GATK here

REFERENCE=/project/uma_lisa_komoroske/Tanya/download/refs/mLynCan4_v1.p/GCF_007474595.1_mLynCan4_v1.p_genomic.fa
Chrom=chr$SGE_TASK_ID
Individual=$1

bsub -q long -W 48:00 -R rusage[mem=1000] -n 1 "samtools mpileup -Q 30 -q 30 -u -v \
-f /project/uma_lisa_komoroske/Tanya/download/refs/mLynCan4_v1.p/GCF_007474595.1_mLynCan4_v1.p_genomic.fa /project/uma_lisa_komoroske/Tanya/analyses/bams_mLynCan4_v1/LRK22_RemoveBadReads.bam |  
bcftools call -c |  
vcfutils.pl vcf2fq -d 3 -D 500 -Q 30 > LRK22.fq"
