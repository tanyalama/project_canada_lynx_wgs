---
title: 'generalized dissimilarity modeling'
disqus: hackmd
---
By: Tanya Lama

These are some additional details re: generalized dissimilarity modeling technique. GDM is part of our gene x environment association analysis pipeline, preceded by LFMM and ultimately used to quantify adaptive potential and measure genetic offset across the landscape. 

Generalized Dissimilarity Modeling GDM
===

## Table of Contents

[TOC]

## Predicting biological distances between sites
The gdm predict function requires a site-pair table in the same format as that used to fit the model. For this we will run the same gdm.input lines, but with clim.points.future (our projected climate data for 2061-2080). Luckily we have the same set of environmental predictors available at our locations for a future time point, so this is fairly straight forward. 

## #1. SNP data
Remains the same, snp from snp.forR above

## #2. Future climate data (2061-2080)
```{r eval=FALSE, include=FALSE}
clim.points.future <- read.csv("/project/uma_lisa_komoroske/Tanya/scripts/GDM/LFMM_subset_snpset_0.1/clim.points.future.csv", header = T, sep=",")
rownames(clim.points.future)<- rownames(clim.points) #looks good

clim.points.future$ID<-c("L155.variant4.variant2",
"LIC11.variant4.variant2",
"LIC20.variant4.variant2",
"LIC23.variant4.variant2",
"LIC24.variant4.variant2",
"LIC27B.variant5.variant2",
"LIC28.variant5.variant2",
"LIC31.variant5.variant2",
"LIC32.variant5.variant2",
"LIC36.variant5.variant2",
"LIC46.variant5.variant2",
"LIC47.variant5.variant2",
"LIC48.variant6.variant2",
"LIC54.variant6.variant2",
"LIC57.variant6.variant2",
"LIC60.variant6.variant2",
"LIC8.variant4.variant2",
"LIC9.variant4.variant2",
"LIT2.variant6.variant2",
"LIT5.variant6.variant2",
"LRK10.variant6.variant2",
"LRK11.variant7.variant2",
"LRK12.variant7.variant2",
"LRK13.variant7.variant2",
"LRK17.variant7.variant2",
"LRK22.variant7.variant2",
"a109.variant.variant.variant2",
"a182.variant.variant.variant2",
"a202.variant.variant.variant2",
"a33.variant.variant.variant2",
"a475.variant.variant.variant2",
"a494.variant.variant.variant2",
"a507.variant.variant.variant2",
"a697.variant","a772.variant",
"a794.variant",
"a803.variant",
"a818.variant",
"a857.variant",
"b114.variant",
"b124.variant2.variant.variant2",
"b13.variant2.variant2",
"b188.variant2.variant.variant2",
"b23.variant2.variant.variant2",
"b276.variant2.variant.variant2",
"b554.variant2.variant.variant2",
"b90.variant2.variant.variant2",
"c165.variant2.variant.variant2",
"c323.variant2.variant2",
"c548.variant2.variant2",
"cb15.variant2.variant2",
"cb42.variant2.variant2",
"cb7.variant2.variant2",
"f264.variant2.variant2",
"f457.variant3.variant2",
"fha_024.variant3.variant2",
"fha_042.variant3.variant2",
"fha_043.variant3.variant2",
"l09_003.variant3.variant2",
"l09_007.variant3.variant2",
"l09_015.variant3.variant2")

#We had to adjust the colnames and order of variables in clim.points.future because they didn't match clim.points (current). The naming and order of the future variables was slightly different. We did this in excel because I made a few errors reordering things in R the first time around (sorry). 
```
## #3. SNP distance matrix
snp.dist.1 remains the same

## #4. Prepare the gdm.input.future
```{r eval=FALSE, include=FALSE}
gdm.input.future <- formatsitepair(bioData=snp.dist.1, bioFormat=3, predData=clim.points.future, siteColumn="ID", XColumn="Longitude", YColumn="Latitude") 
```
## #5. Predict biological distance between 1970-2000 current and 2061-2080 future
```{r eval=FALSE, include=FALSE}
gdm.1.pred<- predict(gdm, gdm.input.future) #gdm is from our first run with current climate data

#Plot predicted vs. observed compositional dissimilarity
par(mfrow=c(1,1))
plot(gdm.input.future$distance, gdm.1.pred, xlab="Observed dissimilarity",
ylab="Predicted dissimilarity", xlim=c(0,1), ylim=c(0,1), pch=20, col=rgb(0,0,1,0.5))
lines(c(-1,2), c(-1,2))

#What does it mean that the points are horizontal to the predicted line...
```
## Variable Importance
Let's use the gdm.varImp function again to assess variable importance. 
```{r eval=FALSE, include=FALSE}
gdm.importance.future <- gdm.varImp(gdm.input.future, geo=T, splines=NULL, nPerm=1, parallel=F, cores=1) #Percent explained is 17%

#pdf("future_GDM_VariableImportance.pdf")
barplot(sort(gdm.importance.future[[2]][,1], decreasing=T), horiz=F, las=2, cex.axis=0.8, cex.names=0.5, ylab="Percent Deviance Explained")
#dev.off()

#This is a different order from our gdm on current climate data..
```
## Turnover Functions
Similar to GF, we can get "turnover functions" showing the relationship of genetic composition with geographic and climate gradients.

These are the "turnover functions" described in Fitzpatrick et al. showing how allelic composition changes along the spatial or environmental gradients. The shapes are nonlinear and large jumps show steep genetic changes along certain portions of the environmental gradient. The height that the function achieved on the right side of the plot is the total importance and should match the barplot. First, organize the variables by importance and then plot:
```{r eval=FALSE, include=FALSE}
plot(gdm.1.pred, plot.layout = c(3, 4))
```
What do these plots suggest? How do the variable importance and turnver functions compare to results from GF?

## GDM_TurnoverFunctions_Uncertainty
You may also wonder how conifdent we are in the modeled relationships. With gdm, we can generate confidence intervals using the plotUncertainty function:
```{r eval=FALSE, include=FALSE]}
pdf("GDM_TurnoverFunctions_Uncertainty.pdf")
plotUncertainty(gdm.input, sampleSites=0.70, bsIters=100, geo=T, plot.layout=c(3,4)) #2 pages because there's too many
dev.off()
```

## Predicting biological change through time
The gdm predict funtion can make predictions across time, for example under climate change scenarios to estimate the magnitude of expected change in biological composition in response to environmental change. In this case, we have to provide rasters for the two time periods of interest. 

## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::

###### tags: `chapter3_adaptation`
