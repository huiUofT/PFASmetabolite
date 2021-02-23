
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PFASMetabolite

<!-- badges: start -->

<!-- badges: end -->

The goal of PFASMetabolite is to detect peak features of PFAS
metabolites from zebrafish

## Installation

You can install the the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("huiUofT/PFASmetabolite")
```

## Example

This is a basic example which shows you how to solve a common problem:

Step 1: loading packages

``` r
  rm(list=ls())
  library(xcms)
  library(devtools)
  library(PFASMetabolite)
  
  #'set up the path
  setwd("C:/Rprogram/EHP_PFAS/PFASmetabolite")
```

Step 2: extracting peak features from raw mass spectrometry data

Step 3: finding out peaks features significantly higher in treatment
groups than controls

``` r
  #'Defining the keywords to extract data from control or treatment
  control<-'S_6'
  treat<-'S18_6'
  
  #'10 is the fold change between control and treatment groups, 0.05 is the p value cutoff
  peaks<-DiffPeaks(peak, 10, 0.05,control,treat)
```

Step 4: identifying primary isotopic peaks, and remove redudant isotopic
peaks

``` r
  #2 ppm is the mass tolerance to find isotopic peaks
  #'0.80 is the correlation coefficient to extract isotopic peaks
  #'3 is the intensity ratio cutoff between primary and other isotopic peaks
  peaks.iso<-FindIsotope(peaks,2,0.80,2)
#> [1] "MS file..." "1"         
#> [1] "MS file..." "2"         
#> [1] "MS file..." "3"         
#> [1] "MS file..." "4"         
#> [1] "MS file..." "5"         
#> [1] "MS file..." "6"         
#> [1] "MS file..." "7"         
#> [1] "MS file..." "8"
```
