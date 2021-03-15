# RadAR

Radiomics Analysis with R (RadAR) is a package for R to perform comprehensive analysis of high-dimensional radiomic datasets.

## Introduction

Quantitative analysis of biomedical images, referred to as radiomics, 
is emerging as a promising approach to facilitate clinical decisions 
and improve patient stratification. The typical radiomic workflow 
includes image acquisition, segmentation, feature extraction, 
and analysis of high-dimensional datasets. While procedures for primary radiomic 
analyses have been established in recent years, processing the resulting radiomic datasets 
remains a challenge due to the lack of specific tools for doing so. 

Here we present RadAR (Radiomics Analysis with R), a new software to perform 
comprehensive analysis of radiomic features. RadAR allows users to process radiomic 
datasets in their entirety, from data import to feature processing and visualization, 
and implements multiple statistical methods for analysis of these data. 
We used RadAR to analyse the radiomic profiles of more than 850 cancer patients 
from publicly available datasets and showed that it was able to recapitulate expected results. 
These results demonstrate RadAR as a reliable and valuable tool for the radiomics community.

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

To include vignettes use:

```{r, eval = FALSE}
# install.packages("devtools")
devtools::install_github("cgplab/RadAR", build_vignettes = T)
```

## Usage 

To show how to use RadAR for the analysis of radiomic datasets, a step-by-step tutorial is included in the package.

## Citation

Benelli M, Barucci A, Zoppetti N, Calusi S, Redapi L, Della Gala G, Piffer S, Bernardi L, Fusi F, Pallotta S. Comprehensive analysis of radiomic datasets by RadAR. Cancer Research 80 (15), 3170-3174. https://cancerres.aacrjournals.org/content/80/15/3170
