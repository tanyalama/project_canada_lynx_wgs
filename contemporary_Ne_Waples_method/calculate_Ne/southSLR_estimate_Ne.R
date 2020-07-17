source('/home/tl50a/Ne_Waples/R_functions.r')
Ne<- get_Ne("/project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_LD/southSLR/southSLR_NC_044317.1")
write.table(Ne$Ne_est, "/project/uma_lisa_komoroske/Tanya/scripts/Ne_Waples/work/calculate_Ne_est/southSLR/southSLR_NC_044317.1.Ne", sep = '\t')
