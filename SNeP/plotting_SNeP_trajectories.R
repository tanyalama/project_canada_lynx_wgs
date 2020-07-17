---
title: "SNeP: trends in recent effective population size trajectories using genome-wide SNP data"
author: "Tanya Lama"
date: "7/10/2020"
output: html_document
editor_options: 
  chunk_output_type: console
---
## Objective
We're interested in looking at trends in "contemporary" demographic history (most recent ~1000 generations) that would be complementary to PSMC and likely more of interest to our state agency folks than PSMC. SNeP and IBDNe are two methods to do this, but SNeP has proven easier to use because IBDNe requires a full pedigree and linkage map which we don't have available at this time.

You can see our full implementation of SNeP here: https://hackmd.io/@tlama/SNeP 

This script is dedicated to visualizations and sensitivity analysis. 

1. We initially ran SNeP on just two chromosomes to explore the resultant output files. SNeP ran successfully on ped/map files from genome-wide SNP data for 54 Canada lynx from three populations. The results are presented here: 
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
lynx_2chrom_test<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/54lynx_2chrom_test/54lynx_2chrom_test.NeAll", header=TRUE)

NW_022059695.1<- read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/allsamples_NW_022059695.1.NeAll", header=TRUE) 
NW_022059697.1<- read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/allsamples_NW_022059697.1.NeAll", header=TRUE) 
allsamples_22chrs<- read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/allsamples_22chrs.NeAll", header=TRUE)
NW_022059735.1<- read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/allsamples_NW_022059735.1.NeAll", header=TRUE) 
NFLD_NW_022059695<- read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/NFLD/NFLD_NW_022059695.1.NeAll", header=TRUE)
NFLD_NC_044303.1<- read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/NFLD/NFLD_NC_044303.1.NeAll", header=TRUE)




#Plot the Ne Trend
# Load Dependencies
library(ggplot2)
library(dplyr)

# Plot
lynx_2chrom_test %>%
  tail(10) %>%
  ggplot( aes(x=GenAgo, y=Ne)) +
    geom_line() +
    geom_point()

ggplot(lynx_20chrom_test, aes(x=GenAgo, y=Ne), ADD=TRUE) +
    geom_line() +
    geom_point()

plot( lynx_20chrom_test$Ne~lynx_20chrom_test$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=rgb(0.2,0.4,0.1,0.7) , lwd=3 , pch=1 , ylim=c(0,100211), )
lines(lynx_2chrom_test$Ne~lynx_2chrom_test$GenAgo , col=rgb(0.8,0.4,0.1,0.7) , lwd=3 , pch=19 , type="b")

plot(NW_022059695.1$Ne~NW_022059695.1$GenAgo,xlim=c(0,1000))


#this is probably the trend line we actually want
plot(NW_022059695.1$Ne~NW_022059695.1$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=c("black") , lwd=1 , pch='o', xlim=c(0,1000), ylim=c(0,100000), axes=FALSE, cex =0.7)
axis(2, at=NW_022059695.1$Ne, col.axis="red", las=1, cex.axis=0.7) #yes!
axis(1, at=NW_022059695.1$GenAgo, col.axis="red", las=2, cex.axis=0.7) #yes
lines(NFLD_NC_044303.1$Ne~NFLD_NC_044303.1$GenAgo, col=c("blue") , lwd=1 , pch=19 , type="c")
legend("topleft", c("southSLR", "Newfoundland"), col = c("black", "blue"),
       text.col = "green4", lty = c(1, 2), cex=0.65, bty = "n")


lines(NW_022059695.1$Ne~NW_022059695.1$GenAgo, col=rgb(0.9,0.4,0.5,0.7) , lwd=3 , pch=19 , type="b")
lines(NW_022059697.1$Ne~NW_022059697.1$GenAgo, col=rgb(0.9,0.1,0.6,0.7) , lwd=3 , pch=19 , type="b")
lines(allsamples_22chrs$Ne~allsamples_22chrs$GenAgo, col=rgb(0.1,0.9,0.9,0.7) , lwd=3 , pch=19 , type="b")
lines(NW_022059735.1$Ne~NW_022059735.1$GenAgo, col=c("black") , lwd=3 , pch=19 , type="c") #makes no sense

lines(NFLD_NW_022059695$Ne~NFLD_NW_022059695$GenAgo, col=c("red") , lwd=1 , pch=19 , type="c") #NFLD a flat line
```
#Plotting sensitivity analysis 
```{r}
#Reading in tables from the sensitivity analysis
sens7<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/sensitivity_analysis/4_sensitivity_itemsTH1.NeAll", header=TRUE)
sens4<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/sensitivity_analysis/7_sensitivity_samplesize_itemsTH_mindist_maxsnp.NeAll", header=TRUE)

plot(sens7$Ne~sens7$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=c("blue") , lwd=1 , pch='o', axes=FALSE, cex =0.7)
axis(2, at=sens7$Ne, col.axis="red", las=1, cex.axis=0.7) #yes!
axis(1, at=sens7$GenAgo, col.axis="red", las=2, cex.axis=0.7) #yes
legend("topleft", c("sens7"), col = c("blue"), text.col = "green4", lty = c(1, 2), cex=0.65, bty = "n")
#Adding 4
lines(sens4$Ne~sens4$GenAgo, col=c("red") , lwd=1 , pch=19 , type="c") 

# I would just drop the obsercation at 332, and then the lines would be pretty comparable!
sens7<-sens7[which(sens7$Ne< 417986),]

#most recent 200 generations actually look pretty good.

```
#Plotting NewFoundland
```{r}
plot(NFLD_NW_022059695reduced$Ne~NFLD_NW_022059695reduced$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=c("blue") , lwd=1 , pch='o', axes=FALSE, cex =0.7)
axis(2, at=NFLD_NW_022059695reduced$Ne, col.axis="red", las=1, cex.axis=0.7) #yes!
axis(1, at=NFLD_NW_022059695reduced$GenAgo, col.axis="red", las=2, cex.axis=0.7) #yes
legend("topleft", c("Newfoundland"), col = c("blue"), text.col = "green4", lty = c(1, 2), cex=0.65, bty = "n")
```
#Plotting a populaiton comparison NFLD southSLR bobcat (missing northSLR)
#Reading in tables
```{r}
bobcat<-read.table("//project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/bobcats_NC_044303.1.NeAll", header=TRUE)
southSLR<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/southSLR_NC_044303.1.NeAll", header=TRUE)
northSLR303b<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/2northSLR_NC_044303.1.NeAll", header=TRUE) 
northSLR304<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/northSLR_NC_044304.1.NeAll", header=TRUE)  
north304full<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/northSLR_full_NC_044303.1.NeAll", header=TRUE)
north303full<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/northSLR_full_NC_044303.1.NeAll", header=TRUE)
northSLR_to_NW_022059693.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/population_comparison/northSLR_to_NW_022059693.1.NeAll", header=TRUE)

#Reading in date for all of southSLR
southSLR_NC_044304.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044304.1.NeAll", header=TRUE)
southSLR_NC_044308.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044308.1.NeAll", header=TRUE) 
southSLR_NC_044312.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044312.1.NeAll", header=TRUE)  
southSLR_NC_044316.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044316.1.NeAll", header=TRUE)  
southSLR_NC_044320.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044320.1.NeAll", header=TRUE)
southSLR_NC_044305.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044305.1.NeAll", header=TRUE)  
southSLR_NC_044309.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044309.1.NeAll", header=TRUE)  
southSLR_NC_044313.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044313.1.NeAll", header=TRUE)  
southSLR_NC_044317.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044317.1.NeAll", header=TRUE)  
southSLR_NC_044321.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044321.1.NeAll", header=TRUE)
southSLR_NC_044306.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044306.1.NeAll", header=TRUE)  
southSLR_NC_044310.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044310.1.NeAll", header=TRUE)  
southSLR_NC_044314.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044314.1.NeAll", header=TRUE)  
southSLR_NC_044318.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044318.1.NeAll", header=TRUE)  
southSLR_NW_022059692.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NW_022059692.1.NeAll", header=TRUE)
southSLR_NC_044307.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044307.1.NeAll", header=TRUE)  
southSLR_NC_044311.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044311.1.NeAll", header=TRUE)  
southSLR_NC_044315.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044315.1.NeAll", header=TRUE)  
southSLR_NC_044319.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/southSLR/southSLR_NC_044319.1.NeAll", header=TRUE)

#plotting southSLR lines
lines(southSLR_NC_044304.1.NeAll$Ne~southSLR_NC_044304.1.NeAll$GenAgo, col=c("red") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044308.1.NeAll$Ne~southSLR_NC_044308.1.NeAll$GenAgo, col=c("blue") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044312.1.NeAll$Ne~southSLR_NC_044312.1.NeAll$GenAgo, col=c("green") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044316.1.NeAll$Ne~southSLR_NC_044316.1.NeAll$GenAgo, col=c("orange") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044320.1.NeAll$Ne~southSLR_NC_044320.1.NeAll$GenAgo, col=c("yellow") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044305.1.NeAll$Ne~southSLR_NC_044305.1.NeAll$GenAgo, col=c("pink") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044309.1.NeAll$Ne~southSLR_NC_044309.1.NeAll$GenAgo, col=c("purple") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044313.1.NeAll$Ne~southSLR_NC_044313.1.NeAll$GenAgo, col=c("violet") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044317.1.NeAll$Ne~southSLR_NC_044317.1.NeAll$GenAgo, col=c("green") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044321.1.NeAll$Ne~southSLR_NC_044321.1.NeAll$GenAgo, col=c("magenta") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044306.1.NeAll$Ne~southSLR_NC_044306.1.NeAll$GenAgo, col=c("darkgreen") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044310.1.NeAll$Ne~southSLR_NC_044310.1.NeAll$GenAgo, col=c("black") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044314.1.NeAll$Ne~southSLR_NC_044314.1.NeAll$GenAgo, col=c("darkblue") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044318.1.NeAll$Ne~southSLR_NC_044318.1.NeAll$GenAgo, col=c("darkred") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NW_022059692.1.NeAll$Ne~southSLR_NW_022059692.1.NeAll$GenAgo, col=c("darkorange") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044307.1.NeAll$Ne~southSLR_NC_044307.1.NeAll$GenAgo, col=c("black") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044311.1.NeAll$Ne~southSLR_NC_044311.1.NeAll$GenAgo, col=c("grey") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044315.1.NeAll$Ne~southSLR_NC_044315.1.NeAll$GenAgo, col=c("gray") , lwd=1 , pch=19 , type="c") 
lines(southSLR_NC_044319.1.NeAll$Ne~southSLR_NC_044319.1.NeAll$GenAgo, col=c("hotpink") , lwd=1 , pch=19 , type="c") 
  
#Plotting northSLR trend lines
northSLR_NC_044305.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044305.1.NeAll", header=TRUE)
northSLR_NC_044309.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044309.1.NeAll", header=TRUE)
northSLR_NC_044306.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044306.1.NeAll", header=TRUE)
northSLR_NC_044310.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044310.1.NeAll", header=TRUE)
northSLR_NC_044307.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044307.1.NeAll", header=TRUE)
northSLR_NC_044311.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044311.1.NeAll", header=TRUE)
northSLR_NC_044308.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044308.1.NeAll", header=TRUE)
northSLR_NC_044312.1.NeAll<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/northSLR/northSLR_NC_044312.1.NeAll", header=TRUE)

lines(northSLR_NC_044305.1.NeAll$Ne~northSLR_NC_044305.1.NeAll$GenAgo, col=c("red") , lwd=1 , pch=19 , type="b") 
lines(northSLR_NC_044309.1.NeAll$Ne~northSLR_NC_044309.1.NeAll$GenAgo, col=c("green") , lwd=1 , pch=19 , type="b") 
lines(northSLR_NC_044306.1.NeAll$Ne~northSLR_NC_044306.1.NeAll$GenAgo, col=c("yellow") , lwd=1 , pch=19 , type="b") 
lines(northSLR_NC_044310.1.NeAll$Ne~northSLR_NC_044310.1.NeAll$GenAgo, col=c("purple") , lwd=1 , pch=19 , type="b") 
lines(northSLR_NC_044307.1.NeAll$Ne~northSLR_NC_044307.1.NeAll$GenAgo, col=c("pink") , lwd=1 , pch=19 , type="b") 
lines(northSLR_NC_044311.1.NeAll$Ne~northSLR_NC_044311.1.NeAll$GenAgo, col=c("magenta") , lwd=1 , pch=19 , type="b") 
lines(northSLR_NC_044308.1.NeAll$Ne~northSLR_NC_044308.1.NeAll$GenAgo, col=c("lightblue") , lwd=1 , pch=19 , type="b") 
lines(northSLR_NC_044312.1.NeAll$Ne~northSLR_NC_044312.1.NeAll$GenAgo, col=c("black") , lwd=1 , pch=19 , type="c") 

  

plot(bobcat$Ne~bobcat$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=c("blue") , lwd=1 , pch='o', axes=FALSE, cex =0.7)
axis(2, at=bobcat$Ne, col.axis="red", las=1, cex.axis=0.7) #yes!
axis(1, at=bobcat$GenAgo, col.axis="red", las=2, cex.axis=0.7) #yes
legend("topleft", c("bobcat"), col = c("blue"), text.col = "green4", lty = c(1, 2), cex=0.65, bty = "n")

plot(southSLR$Ne~southSLR$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=c("blue") , lwd=1 , pch='o', axes=FALSE, cex =0.7)
axis(2, at=southSLR$Ne, col.axis="black", las=1, cex.axis=0.7) #yes!
axis(1, at=southSLR$GenAgo, col.axis="black", las=2, cex.axis=0.7) #yes
legend("topleft", c("southSLR", "NFLD", "bobcat"), col = c("blue","orange", "red"), text.col = "green4", lty = c(1, 2), cex=0.65, bty = "n") 

#dropping this peak 
southSLR<-southSLR[which(southSLR$Ne< 60542),]
southSLR<-southSLR[which(southSLR$Ne< 10858 ),]
southSLR<-southSLR[which(southSLR$Ne< 8681 ),]

#adding lines
lines(bobcat$Ne~bobcat$GenAgo, col=c("red") , lwd=1 , pch=19 , type="c") 
lines(NFLD_NC_044303.1$Ne~NFLD_NC_044303.1$GenAgo, col=c("orange") , lwd=1 , pch=19 , type="c")
lines(northSLR$Ne~northSLR$GenAgo, col=c("magenta") , lwd=1 , pch=19 , type="c")
lines(northSLR304$Ne~northSLR304$GenAgo, col=c("green") , lwd=1 , pch=19 , type="c") #this one mirrors the exact trajectory of southSLR. Running the full CHR tonight
lines(northSLR303b$Ne~northSLR303b$GenAgo, col=c("red") , lwd=1 , pch=19 , type="b")
lines(north303full$Ne~north303full$GenAgo, col=c("purple") , lwd=1 , pch=19 , type="b")
lines(southSLR$Ne~southSLR$GenAgo, col=c("green") , lwd=1 , pch=19 , type="c")
lines(northSLR_to_NW_022059693.1.NeAll$Ne~northSLR_to_NW_022059693.1.NeAll$GenAgo, col=c("green") , lwd=1 , pch=19 , type="b")

#we can try this in reverse, using northSLR as the main plot and adding southSLR as a line
plot(north303full$Ne~north303full$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=c("blue") , lwd=3 , pch='o', axes=FALSE, cex =0.7, xlim=c(0,100))
axis(2, at=north303full$Ne, col.axis="black", las=1, cex.axis=0.7) #yes!
axis(1, at=north303full$GenAgo, col.axis="black", las=2, cex.axis=0.7) #yes
legend("topleft", c("northSLR", "southSLR", "NFLD", "bobcat"), col = c("blue","orange", "red"), text.col = "green4", lty = c(1, 2), cex=0.65, bty = "n")

#I think this figure would make more sense without bobcat, and with northSLR

#Zooming in on the first 100 gen
plot(southSLR$Ne~southSLR$GenAgo, type="b" , bty="l" , xlab="Generations Ago" , ylab="Effective Population Size (Ne)" , col=c("blue") , lwd=1 , pch='o', axes=FALSE, cex =0.7, xlim=c(0,100), ylim=c(0,4000))
axis(2, at=southSLR$Ne, col.axis="black", las=1, cex.axis=0.7) #yes!
axis(1, at=southSLR$GenAgo, col.axis="black", las=2, cex.axis=0.7) #yes
legend("topleft", c("southSLR", "NFLD"), col = c("blue","orange", "red"), text.col = "green4", lty = c(1, 2), cex=0.65, bty = "n")

I think our SLR plot was looking good. northSLR isn't looking great. There's very little detail there, but that was just a quick run. Now it's time to run the full suite, with as many chromosomes as we can.






p <- ggplot2::ggplot(NFLD_NW_022059695reduced)
        p <- p + ggplot2::geom_point(data = NFLD_NW_022059695reduced, aes(x = GenAgo, 
            y = Ne, color = c("black"))
            
        p <- p + ggplot2::xlim(0, max(teilsatz$to)) + ggplot2::ggtitle(paste("Chromosome", 
            unique(teilsatz$chrom)))
        p <- p + ggplot2::xlab("Chromosome position (Mbps)")
        p <- p + theme(plot.title = element_text(hjust = 0.5)) +theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), aspect.ratio = .5, axis.text.y=element_blank(), legend.title = element_blank(), legend.position = "none", axis.title.y = element_blank())
        print(p)
grid.arrange(p,q, nrow = 2)




#NFLD a flat line



smooth.spline(NW_022059695.1, y = Ne, w = NULL, df, spar = NULL, lambda = NULL, cv = FALSE,
              all.knots = FALSE, nknots = .nknots.smspl,
              keep.data = TRUE, df.offset = 0, penalty = 1,
              control.spar = list(), tol = 1e-6 * IQR(x), keep.stuff = FALSE)

plot(NFLD_NW_022059695$Ne~NFLD_NW_022059695$GenAgo, main = "Stopping Distance versus Speed")

NFLD_NW_022059695
```

#Sensitivity Analysis
2. Now we have to ask ourselves if SNeP is sensitive to the number of chromosomes (and thus SNPs) included in the analysis, and the number of samples/population. 
Let's try with 20 chrom
```{r}
lynx_20chrom_test<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/54lynx_NC_044303.1_NW_022059692.1.NeAll", header=TRUE)


3. Let's feed SNeP the full suite of chromosomes, and evaluate the output. 

We are having a really hard time with one chromosome in particular, it may be the SNP density (or lack thereof?) of SNPs in the 21st chromosome, that SNeP is getting stuck on. When we remove that 21st chromosome, it carries on pretty seamlessly. 
```{r}
lynx_65chrom<-read.table("/project/uma_lisa_komoroske/Tanya/scripts/SNeP/)

bsub -q short -R rusage[mem=1000] -n 1 -W 1:00 -R span\[hosts=1\] "singularity exec /share/pkg/singularity/images/xenial_16.04LTS.sif /project/uma_lisa_komoroske/bin/SNeP1.1 -ped /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.ped -map /project/uma_lisa_komoroske/Tanya/Rgenomics/mLynCan4_v1.p_lynxonly/6SV_unfiltered_SNPs.map -out /project/uma_lisa_komoroske/Tanya/scripts/SNeP/allsamples_NW_022059694.1 -chr NW_022059694.1 -samplesize 54 -threads 8"
