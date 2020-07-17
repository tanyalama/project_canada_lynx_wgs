source('/home/tl50a/Ne_Waples/R_functions.r')
Ne<- get_Ne("/project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_LD/NFLD/NFLD_NC0443041")
write.table(Ne$Ne_est, "/project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_Ne_est/NFLD/NFLD_NC0443041.Ne", sep = '\t')

#bsub -q long -R rusage[mem=60000] -n 4 -W 12:00 -R span\[hosts=1\] "Rscript ./NFLD_estimate_Ne.r"
