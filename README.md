# RadAR

Radiomics Analysis with R (RadAR) is a package for R to perform comprehensive analysis of high-dimensional radiomic datasets.

## Introduction

The quantitative analysis of biomedical images, referred to as Radiomics, is emerging as a promising approach to facilitate clinical decisions and improve patients’ stratification. The typical radiomic workflow includes image acquisition, segmentation, feature extraction and analysis of high-dimensional datasets. While procedures for primary radiomic analyses have been established during the last years, the processing of the resulting radiomic datasets remains challenging due to the lack of specific tools. 

Here, we present RadAR (Radiomics Analysis with R), a new software to perform comprehensive analysis of radiomic features. RadAR allows the users to carry out the entire processing of radiomic datasets, from data import to feature processing and visualization and implements multiple statistical methods for the analysis of these data. We used RadAR to analyse the radiomic profiles of more than 500 cancer patients from publicly available datasets and showed that it was able to recapitulate expected results, demonstrating its reliability and proving that RadAR may represent a valuable tool for the radiomic community.

## Installation
RadAR  is freely available under MIT license at  <https://github.com/cgplab/RadAR>. 
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

Benelli M, Barucci A, Zoppetti N, Calusi S, Redapi L, Della Gala G, Piffer S, Bernardi L, Fusi F, Pallotta S. Comprehensive analysis of radiomic datasets by RadAR. Cancer Research; in press. https://cancerres.aacrjournals.org/content/early/2020/06/13/0008-5472.CAN-20-0332
