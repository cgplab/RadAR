# RadAR

Radiomics Analysis with R (RadAR) is a package for R to perform comprehensive analysis of high-dimensional radiomic datasets.

## Introduction

The quantitative analysis of biomedical images, referred to as Radiomics, is emerging as a promising approach to facilitate clinical decisions and improve patients’ stratification. A typical radiomic workflow includes image acquisition, segmentation, feature extraction and statistical analysis of high-dimensional datasets (Gillies et al., 2016). While procedures for primary radiomic analyses have been established during the last years, processing of resulting radiomic datasets remains challenging due to the lack of dedicated tools.
Here, we present RadAR (Radiomics Analysis with R), a new software to perform comprehensive analysis of radiomic features. RadAR allows the users to carry out the entire processing of radiomic datasets, from data import to feature processing and visualization and implements statistical functions tailored to these kinds of data. We applied RadAR to two radiomic datasets and showed that it was able to recapitulate expected results, demonstrating its reliability and proving that RadAR may represent a valuable tool for the radiomic community. 

## Installation
RadAR  is freely available under GPL-3 license at  <https://github.com/cgplab/RadAR>. 
First, install biocViews to facilitate the installation of the package dependecies:

```{r, eval=FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("biocViews")
```

Then, you can install RadAR from github with:

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("cgplab/RadAR")
```

To include vignette use:

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("cgplab/RadAR", build_vignettes = T)
```

## Usage 

To show how to use RadAR for the analysis of radiomic datasets, a step-by-step tutorial is included in the package.

## Citation

<<<<<<< HEAD
Benelli M, Barucci A, Zoppetti N, Calusi S, Redapi L, Della Gala G, Piffer S, Bernardi L, Fusi F, Pallotta S. Comprehensive analysis of radiomic datasets by RadAR. *Submitted*.
