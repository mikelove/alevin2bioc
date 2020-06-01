# alevin2bioc

## Introduction

This package is a workflow for 
[BioC2020](https://bioc2020.bioconductor.org/)
and provides an online vignette describing how to import 
[alevin](https://salmon.readthedocs.io/en/latest/alevin.html)
scRNA-seq quantifications into R/Bioconductor. The developers of the
workflow are listed in the sidebar. 

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
