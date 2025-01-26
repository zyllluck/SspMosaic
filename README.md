# SSpMosaic
## Robust integration and annotation of single-cell and spatial multi-omics data using interpretable gene programs 

<p align="center">
<img src="https://github.com/zyllluck/SSpMosaic/blob/main/Workflow.png" width="800" />
</p>
The integration and annotation of single-cell and spatial multi-omics data have revolutionized our understanding of cellular heterogeneity and tissue microenvironments. Gene programs, representing co-regulated or co-expressed gene sets reflective of cellular functions, offer a biologically interpretable means to bridge diverse data types. Here, we present SSpMosaic, a novel framework that leverages interpretable metaprograms — higher-order representations of gene programs learned across datasets — to achieve robust data integration and cell-type annotation. By aligning gene programs across diverse single-cell data, SSpMosaic ensures consistent and accurate integration across batches, modalities, and species. Through learning metaprograms from reference datasets, SSpMosaic accurately annotates cell types within query datasets, enabling discovery and annotation of novel cell types through gene program analysis. Furthermore, SSpMosaic extends to spatial transcriptomics, allowing cell-type deconvolution, spatial domain detection, and identification of spatially resolved gene programs to reveal tissue organization and cellular interactions. Finally, SSpMosaic can be adapted for reference-free spatial organization characterization, facilitating the discovery of novel tissue structures and cellular niches without relying on reference single-cell datasets.

Installation
------------
You can install the released version of SSpMosaic from Github with the following code.

## Dependencies 
* R version >= 4.0.0.
* Dependent R packages: nnls, Seurat, SeuratObject, WGCNA, methods, Matrix, RcppAnnoy, hdWGCNA, igraph, qlcMatrix, stats, ggplot2, grDevices, uwot, utils, stringr, patchwork, dplyr

``` r
# install devtools if necessary
install.packages('devtools')

# install the SSpMosaic package
devtools::install_github('zyllluck/SSpMosaic')

# load package
library(SSpMosaic)

```


How to cite `SSpMosaic`
-------------------

How to use `SSpMosaic`
-------------------
SSpMosaic tutorials are as follows.

## SSpMosaic program

SSpMosaic program generation tutorial is in [generate_module](https://zyllluck.github.io/SSpMosaic/program_generation.html)

SSpMosaic module network-propagation tutorial is in [network-propagation](https://zyllluck.github.io/SSpMosaic/network_propagation.html)

## SSpMosaic downstream functions

SSpMosaic integration tutorial is in [integration](https://zyllluck.github.io/SSpMosaic/integration_tutorial.html)

SSpMosaic single-cell annotation tutorial is in [sc_annotation](https://zyllluck.github.io/SSpMosaic/sc_annotation_tutorial.html)

SSpMosaic deconvolution tutorial is in [deconvolution](https://zyllluck.github.io/SSpMosaic/spatial_deconvolution_tutorial.html)

Issues
------------
All feedback, bug reports and suggestions are warmly welcomed! Please make sure to raise issues with a detailed and reproducible example and also please provide the output of your sessionInfo() in R! 
