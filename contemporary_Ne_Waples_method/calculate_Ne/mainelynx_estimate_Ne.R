source('/home/tl50a/Ne_Waples/R_functions.r')
Ne<- get_Ne("/project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_LD/mainelynx/mainelynx_NC_044318.1")
write.table(Ne$Ne_est, "/project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_Ne_est/mainelynx/mainelynx_NC_044318.1.Ne", sep = '\t')

#bsub -q long -R rusage[mem=60000] -n 4 -W 6:00 -R span\[hosts=1\] "Rscript ./mainelynx_estimate_Ne.r"
