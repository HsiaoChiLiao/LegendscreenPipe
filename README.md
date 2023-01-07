# Analysing single-cell protein intensities with LegendscreenPipe (Beta-version)
### Author: Hsiao-Chi Liao
### Date: 7 Jan 2022
<br>

# Introduction

## For data from the LegendScreen assay:

This end-to-end toolbox carefully pre-processes the raw data in FCS
format (one file per well), and further imputes the ‘missing’
exploratory markers (Legend) in the wells without measurement. The
pipeline starts by performing background correction on raw intensities
to remove the noise from electronic baseline restoration and
fluorescence compensation by adapting the normal-exponential convolution
model. Secondly, unwanted technical variation such as well effects is
removed by applying the log-normal model with plate, column, and row
factors. Thirdly, imputation is done by using the informative backbone
markers as predictors. Lastly, cluster analysis is performed on both
normalised backbone measurements and the completed dataset.

## For data from the fluorescence flow cytometry (FFC) assay:

This end-to-end toolbox carefully pre-processes the raw data in FCS
format (one file per well or batch). The pipeline starts by performing
background correction on raw intensities to remove the noise from
electronic baseline restoration and fluorescence compensation by
adapting the normal-exponential convolution model. Secondly, unwanted
technical variation such as batch/well effects is removed by applying
the log-normal model with batch information or plate, column, and row
factors. Lastly, cluster analysis is performed on the normalised
backbone measurements.

## Notes:

\*For the plate-based experiments, the file names of the FCS files MUST
contain the following specified format (if set file_meta == "auto"): <br> Plate information: Plate1, Plate2,
…, Plate9 <br> Well information: A1, A2, …, A12, B1, …, H1, …, H12

\*For `file_meta == "usr"`, prepare `filename_meta.csv` in the following format and save the CSV file in `FCSpath` 
|  Filenam  |  Plate  |  Well  |  Column  |  Row  |  Well.lab  |
|:---------:|:-------:|:------:|:--------:|:-----:|:----------:|
|p1_a12.fcs |  Plate1 |   A12  |  Col.12  | Row.01|   P1_A12   |
|p2_d08.fcs |  Plate2 |   D08  |  Col.08  | Row.04|   P2_D08   |
|p3_g1.fcs  |  Plate3 |   G01  |  Col.01  | Row.07|   P3_G01   |
<br>
Note that the "Filenam" column refers to the file names of the FCS files in the `FCSpath`.

\*The function `pipBackbone` is designed for normalising data from the
fluorescence flow cytometry (FFC) assay. The main difference between the
functions `pipLegendscreen` and `pipBackbone` is that `pipBackbone` does
not contain the ‘exploratory’ `Legend` markers, and thus does not
require imputation. Users may use `pipBackbone` for normalising data from the LegendScreen assay without doing imputation.
<br>

**COMING SOON**: <br> Experiments (FFC) from different batches… <br> 
Prepare `filename_meta.csv` in the following format and save the CSV file in `FCSpath`
|  Filenam  |  Batch  |
|:---------:|:-------:|
|090122.fcs |  Batch1 |
|070122.fcs |  Batch2 |
|010122.fcs |  Batch3 |
<br>

# Four sections in this Vignette

Section 0 (S0). Required packages <br> Section 1 (S1). Installing `LegendscreenPipe` and downloading the sample dataset <br> Section 2 (S2). Analysing the single-cell
protein intensites with `LegendscreenPipe` <br> Section 3 (S3).
Description of the output

## S0. Loading required packages and checking if they’re installed

Along with the `LegendscreenPipe` package, we will also load the
following packages required for the whole analysis.

``` r
library(flowCore)
library(Biobase)
library(stringr)
library(uwot)
library(Rphenograph)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(Rfast)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(glmnetUtils)
library(e1071)
library(xgboost)
library(foreach)
library(doParallel)
library(parallel)
library(pbapply)
library(reshape)
library(gtools)
```

## S1. Installing `LegendscreenPipe` and downloading the sample dataset

### Please install the package from GitHub.

``` r
if (!requireNamespace("devtools", quietly = TRUE)){
    install.packages("devtools")
}
devtools::install_github("HsiaoChiLiao/LegendscreenPipe", build_vignettes=FALSE)
```
### The sample dataset is stored in the folder 'data-raw' on GitHub
https://github.com/HsiaoChiLiao/LegendscreenPipe.git
<br>
This is a subset (1000 cells/well) of data downloaded from Becht et al., 2021 (DOI: 10.1126/sciadv.abg0505).
<br>

## S2. Analysing the single-cell protein intensites with `LegendscreenPipe`

``` r
library(LegendscreenPipe)
```

### `pipLegendscreen`: for data from the LegendScreen assay

``` r
# running with build-in FCS files:
pipLegendscreen(
FCSpath="/PathToThePackage/LegendscreenPipe/data-raw/",  #must contain 'LegendscreenPipe/data-raw'
Outpath="/OutputPath/Output/",  #the folder for saving output should be named 'Output'
file_meta="auto",
bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
"F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"), 
chans = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
"F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"), 
yvar="Legend", control.wells = c("P1_A01", "P2_A02", "P3_G02"), 
bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01, 
pe.mean.sd=3, pe.lower.quntile=0.1, pe.min.quntile=0.01,
trans.dat="cent.lgc", 
models.use = c("XGBoost"),
extra_args_regression_params = list(list(nrounds = 1500, eta = 0.03)),
cores=2L)

# check the usage
help(pipLegendscreen, package = "LegendscreenPipe")
```

### `pipBackbone`: for data from the fluorescence flow cytometry (FFC) assay

``` r
# running with build-in FCS files:
pipBackbone(
FCSpath="/PathToThePackage/LegendscreenPipe/data-raw/",  #must contain 'LegendscreenPipe/data-raw'
Outpath="/OutputPath/Output/",  #the folder for saving output should be named 'Output'
file_meta="auto",
plate_based = TRUE,
bkb.v = c("FSC-H", "FSC-W", "SSC-H", "SSC-W", "CD69-CD301b", "MHCII", "CD4", "CD44", "CD8", "CD11c", "CD11b",
"F480", "Ly6C", "Lineage", "CD45a488", "CD24", "CD103"), 
bkb.upper.quntile=0.9, bkb.lower.quntile=0.1, bkb.min.quntile=0.01, 
trans.dat="cent.lgc")

# check the usage
help(pipBackbone, package = "LegendscreenPipe")
```

All the output will be stored in the `Outpath`.

## S3. Description of the output

Three folders will be automatically generated in the `Output` folder.
<br> <br> 1. `intermediary`: Intermediary results will be saved in the
.rds and .RData formats and will be stored here. <br> <br> 2.
`downstream`: Final results will be saved in the .rds format and will be
stored here. The results include normalised backbone measurements (on
both linear and log scale), the completed dataset with imputed Legend
markers, UMAP coordinates derived from both normalised backbones and the
completed dataset, and metadata for cells including cluster labels
derived from both normalised backbones and the completed dataset. <br>
<br> 3. `graph`: Figures will be stored here, including **scatter
plots** for comparing background corrected and raw intensities for each
protein marker, **heatmaps** for presenting the biological and unwanted
effects in the data before and after removal of unwanted variation,
**boxplots** of R-sq for visualising the accuracy of imputed Legend
markers, and **UMAP plots** for showing the cluster structure. <br> <br>
<br>

# THANK YOU

Thanks for carrying out analyses with our `LegendscreenPipe` package.
Please feel free to raise any issues and/or questions and send them to
[Hsiao-Chi’s email](hsiaochi.liao@student.unimelb.edu.au).

## Information about the R session when this vignette was built

``` r
sessionInfo()
#> R version 4.2.0 (2022-04-22)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur/Monterey 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.2.0  magrittr_2.0.3  fastmap_1.1.0   cli_3.4.1      
#>  [5] tools_4.2.0     htmltools_0.5.3 rstudioapi_0.14 yaml_2.3.6     
#>  [9] stringi_1.7.8   rmarkdown_2.18  knitr_1.41      stringr_1.4.1  
#> [13] xfun_0.35       digest_0.6.30   rlang_1.0.6     evaluate_0.18
```
