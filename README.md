# alevin2bioc

![](https://github.com/mikelove/alevin2bioc/workflows/.github/workflows/basic_checks.yaml/badge.svg)

Build URL: https://mikelove.github.io/alevin2bioc

Docker image: https://hub.docker.com/repository/docker/mikelove/alevin2bioc

## Introduction

This package is a workflow for 
[BioC2020](https://bioc2020.bioconductor.org/)
and provides an online vignette describing how to import 
[alevin](https://salmon.readthedocs.io/en/latest/alevin.html)
scRNA-seq quantifications into R/Bioconductor. The developers of the
workflow are listed in the sidebar. 

## Packages

The following software are used in this workflow:

<img width="150" alt="alevin" src="https://i.imgur.com/Y9VPCsR.png"/> <img width="150" alt="tximeta" src="https://github.com/Bioconductor/BiocStickers/blob/master/tximeta/tximeta.png?raw=true"/> <img width="150" alt="SingleCellExperiment" src="https://github.com/Bioconductor/BiocStickers/blob/master/SingleCellExperiment/SingleCellExperiment.png?raw=true"/> <img width="150" alt="fishpond" src="https://github.com/Bioconductor/BiocStickers/blob/master/fishpond/fishpond.png?raw=true"/> <img width="150" alt="scran" src="https://github.com/Bioconductor/BiocStickers/blob/master/scran/scran.png?raw=true"/> <img width="150" alt="Seurat" src="https://i.imgur.com/FEFIXCc.jpeg"/>



## Vignette

The workflow vignette can be accessed by clicking `Get started` in the
top navigation bar.

## Running the workflow on your machine

The vignette is designed for R `4.0` and the `devel` branch of
Bioconductor. The easiest way to try out the workflow is to use the
Docker image (see section below), which has all the software
pre-configured to the correct versions.  The workflow should also work
with the `release` branch of Bioconductor (v3.11) except for one
package, *fishpond*, which has additional functionality demonstrated
here that only is supported in the `devel` branch. If you are using
the `release` branch on your laptop, and want to try this workflow, it
should be sufficient to install `release` branch of all packages, and
then only updating *fishpond* to `devel` with:

```
devtools::install_github("mikelove/fishpond", dependencies=FALSE)
```

## Docker image

A Docker image can be run for this workflow:

```
docker run -e PASSWORD=abc -p 8787:8787 mikelove/alevin2bioc
```

Once running, navigate to <http://localhost:8787/> and then login with
`rstudio:yourchosenpassword`.

## Miniconda message

If you load *Seurat* during the workflow, it loads *reticulate* which
triggers the following message:

```
No non-system installation of Python could be found. 
Would you like to download and install Miniconda?
```

Just type `n` in response.
