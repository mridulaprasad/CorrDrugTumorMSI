# CorrDrugTumorMSI
A R package to correlate drug distribution with tumor tissue types in mass spectrometry imaging (MSI) data. 

## Description
CorreDrugTumorMSI provides a complete analysis pipeline to analyze untargeted drug mass spectrometry imaging (MSI) data. The main functionality of the package includes: preprocessing of MSI data, identified unsupervised clusters, create a drug spatial map, quantitative analysis of identified clusters and drug spatial map and selection of spatially relevant molecular ions from the identified clusters. For further details see the package vignette

## Installation

* Install HomogenMSI: In the shell type

```r
library(devtools)
install_github("pietrofranceschi/HomogenMSI", dependencies = TRUE) 
```

* To additionally build the vignette the previous command should include build_vignettes = TRUE

```r
library(devtools)
install_github("pietrofranceschi/HomogenMSI", dependencies = TRUE, build_vignettes = TRUE) 

```
* Load the package: In shell type
```r
library(HomogenMSI)
```
