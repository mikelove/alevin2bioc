---
output:
  rmarkdown::html_document:
    highlight: pygments
    toc: true
    toc_depth: 3
    fig_width: 5
bibliography: "/Library/Frameworks/R.framework/Versions/4.0/Resources/library/alevin2bioc/vignettes/library.bib"
vignette: >
  %\VignetteIndexEntry{alevin2bioc}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding[utf8]{inputenc}
---

# Tutorial for importing alevin scRNA-seq quantifications into R/Bioconductor

## Instructor(s) name(s) and contact information

[Michael Love](https://mikelove.github.io),
[Avi Srivastava](https://k3yavi.github.io)

## Introduction

*alevin* is a method for ... [@alevin]. It extends the methods in the
*Salmon* software [@salmon].

The data we will use is ... (citation).

## Running alevin

In order to run *alevin*, we must first ...

Instructions on downloading *alevin* can be found here ...

Instructions on indexing a set of reference transcripts can be found
here ... For this experiment, we used the 
[GENCODE](https://www.gencodegenes.org/) 
human reference transcripts [@gencode].

## Importing alevin data with tximeta

We will use *tximeta* ... [@tximeta].

First we specify the path where the quantification data is stored. In
this tutorial, the data is stored in an R package, so we need to use
the `system.file` command. For typical use, you would not use
`system.file`, but would just specify the path to the directory of the
output from *alevin*. 


```r
# normally you just would use:
# dir <- "/path/to/alevin/output"
extdata <- system.file("extdata", package="alevin2bioc")
dir <- file.path(extdata, "pbmc_1k")
```


```r
files <- file.path(dir, "alevin", "quants_mat.gz")
file.exists(files)
```

```
## [1] TRUE
```



We can import the *alevin* quantification using the following call to
`tximeta`. The extra argument to *alevin* will filter out cells based
on *alevin*'s post-quantification quality control methods (see paper
for details).


```r
library(tximeta)
se <- tximeta(files, type="alevin", alevinArgs=list(filterBarcodes=TRUE))
```

```
## importing quantifications
```

```
## reading in alevin gene-level counts across cells with fishpond
```

```
## filtering down to 1142 cell barcodes
```

```
## reading in alevin inferential variance with fishpond
```

```
## found matching transcriptome:
## [ GENCODE - Homo sapiens - release 33 ]
```

```
## loading existing TxDb created: 2020-05-13 16:17:46
```

```
## generating gene ranges
```

```
## loading existing gene ranges created: 2020-05-13 16:17:51
```

```
## fetching genome info for GENCODE
```

`tximeta` returns a *SummarizedExperiment* [@Lawrence2013]. We can
easily convert this object into a *SingleCellExperiment*
[@Amezquita2020] which has specific slots designed for single-cell
experiment data.


```r
suppressPackageStartupMessages(library(SingleCellExperiment))
sce <- as(se, "SingleCellExperiment")
```

## Benefits of tximeta

We can automatically add IDs, because *tximeta* knows the type of
identifiers on the rows of the `sce` object:


```r
library(org.Hs.eg.db)
sce <- addIds(sce, "SYMBOL")
```

```
## mapping to new IDs using 'org.Hs.eg.db' data package
## if all matching IDs are desired, and '1:many mappings' are reported,
## set multiVals='list' to obtain all the matching IDs
```

```
## it appears the rows are gene IDs, setting 'gene' to TRUE
```

```
## 'select()' returned 1:many mapping between keys and columns
```


```r
mcols(sce)
```

```
## DataFrame with 60233 rows and 2 columns
##                             gene_id      SYMBOL
##                         <character> <character>
## ENSG00000223972.5 ENSG00000223972.5     DDX11L1
## ENSG00000227232.5 ENSG00000227232.5          NA
## ENSG00000278267.1 ENSG00000278267.1   MIR6859-1
## ENSG00000243485.5 ENSG00000243485.5          NA
## ENSG00000284332.1 ENSG00000284332.1   MIR1302-2
## ...                             ...         ...
## ENSG00000198695.2 ENSG00000198695.2         ND6
## ENSG00000210194.1 ENSG00000210194.1          NA
## ENSG00000198727.2 ENSG00000198727.2        CYTB
## ENSG00000210195.2 ENSG00000210195.2          NA
## ENSG00000210196.2 ENSG00000210196.2          NA
```

Also, because the provenance was detected, we also have the ranges of
the genes in their proper genomic context. So it is easy to find, for
example, genes near a particular position in the genome, in this case
4 genes that overlap the range `chr1:10,000,000-10,100,000`.


```r
x <- GRanges("chr1", IRanges(10e6,10.1e6))
sce[sce %over% x,]
```

```
## class: SingleCellExperiment 
## dim: 4 1142 
## metadata(6): tximetaInfo quantInfo ... txomeInfo txdbInfo
## assays(3): counts variance mean
## rownames(4): ENSG00000162444.12 ENSG00000130939.20 ENSG00000224340.1 ENSG00000233623.2
## rowData names(2): gene_id SYMBOL
## colnames(1142): AAGACAACAGATCACT CAGGTATGTCAGTCTA ... CTCTGGTAGTAGGATT AGGAGGTAGGTCACAG
## colData names(0):
## reducedDimNames(0):
## altExpNames(0):
```

## Add cell annotations

Cell annotations were generated using *Seurat* [@seurat]. The script is
saved in this package in `inst/scripts/seurat.R`.


```r
ids <- readRDS(file.path(extdata, "idents.rds"))
top10 <- read.csv(file.path(extdata, "top10.csv"))
```


```r
idx <- colnames(sce) %in% names(ids)
table(idx)
```

```
## idx
## FALSE  TRUE 
##    74  1068
```

```r
sce <- sce[,idx]
sce$cluster <- ids[colnames(sce)]
```


```r
table(sce$cluster)
```

```
## 
##           CD4 T CD14+ Monocytes               B           CD8 T              NK 
##             346             346             189             133              54
```

Note that the different clusters have different total counts, for
example:


```r
cs <- colSums(assays(sce)[["counts"]])
# cells with higher number of UMI
more.umi <- cs > 10000
(tab <- table(more.umi, sce$cluster))
```

```
##         
## more.umi CD4 T CD14+ Monocytes   B CD8 T  NK
##    FALSE   322             159 165   130  54
##    TRUE     24             187  24     3   0
```

```r
100 * round(prop.table(tab,2),2) # percent
```

```
##         
## more.umi CD4 T CD14+ Monocytes   B CD8 T  NK
##    FALSE    93              46  87    98 100
##    TRUE      7              54  13     2   0
```

```r
# cell with lower number of UMI
fewer.umi <- cs < 5000
(tab <- table(fewer.umi, sce$cluster))
```

```
##          
## fewer.umi CD4 T CD14+ Monocytes   B CD8 T  NK
##     FALSE   316             300 149   107  39
##     TRUE     30              46  40    26  15
```

```r
100 * round(prop.table(tab,2),2) # percent
```

```
##          
## fewer.umi CD4 T CD14+ Monocytes  B CD8 T NK
##     FALSE    91              87 79    80 72
##     TRUE      9              13 21    20 28
```


```r
head(top10)
```

```
##   X         p_val avg_logFC pct.1 pct.2     p_val_adj cluster               gene symbol
## 1 1 3.399988e-135 1.2134065 0.997 0.602 6.075778e-131       0 ENSG00000111716.14   LDHB
## 2 2 2.755855e-116 1.0112296 0.879 0.177 4.924713e-112       0 ENSG00000081059.20   TCF7
## 3 3 4.214635e-114 0.9936307 0.968 0.194 7.531553e-110       0  ENSG00000167286.9   CD3D
## 4 4 1.644166e-111 0.9595734 0.855 0.154 2.938125e-107       0 ENSG00000127152.18 BCL11B
## 5 5 9.344999e-110 1.0665597 0.980 0.266 1.669951e-105       0  ENSG00000277734.8   <NA>
## 6 6 2.602057e-105 1.1441765 0.916 0.175 4.649875e-101       0 ENSG00000168685.15   IL7R
##   cluster.name
## 1        CD4 T
## 2        CD4 T
## 3        CD4 T
## 4        CD4 T
## 5        CD4 T
## 6        CD4 T
```

## Plotting counts with uncertainty

For a later demonstration of scaling, we will sort the cells by the
total count (this is not something you would necessarily do in a
typical analysis).


```r
o <- order(colSums(assays(sce)[["counts"]]), decreasing=TRUE)
sce <- sce[,o]
```


```r
library(fishpond)
plotInfReps(sce, idx="ENSG00000167286.9",
            x="cluster", mainCol="SYMBOL",
            legend=TRUE)
```

<img src="alevin2bioc_files/figure-html/plot-basic-1.png" width="480" />


```r
plotInfReps(sce[,sample(ncol(sce),200)],
            idx="ENSG00000167286.9",
            x="cluster", mainCol="SYMBOL",
            legend=TRUE)
```

<img src="alevin2bioc_files/figure-html/plot-medium-1.png" width="480" />


```r
plotInfReps(sce[,sample(ncol(sce),100)],
            idx="ENSG00000167286.9",
            x="cluster", mainCol="SYMBOL",
            legend=TRUE)
```

<img src="alevin2bioc_files/figure-html/plot-small-1.png" width="480" />


```r
plotInfReps(sce, idx="ENSG00000167286.9",
            x="cluster", mainCol="SYMBOL",
            reorder=FALSE)
```

<img src="alevin2bioc_files/figure-html/plot-no-order-1.png" width="480" />

## Scaling with size factors

We use `computeSumFactors` [@Lun2016] from the *scran* package
[@scran] ...


```r
library(scran)
sce <- computeSumFactors(sce, clusters=sce$cluster)
plot(cs, sizeFactors(sce), xlab="column sum", ylab="sum factor")
```

<img src="alevin2bioc_files/figure-html/size-factors-1.png" width="480" />


```r
par(mfrow=c(2,1), mar=c(2.5,4.5,1,1))
plotInfReps(sce, idx="ENSG00000167286.9",
            x="cluster", main="",
            reorder=FALSE)
plotInfReps(sce, idx="ENSG00000167286.9",
            x="cluster", main="",
            applySF=TRUE, reorder=FALSE)
```

<img src="alevin2bioc_files/figure-html/scaling-1.png" width="480" />

## Inferential variance (uncertainty)

Inferential variance is a focus of the *Swish* nonparametric
statistical method [@swish] which we do not demonstrate here ...


```r
var <- as.vector(assays(sce)[["variance"]])
mu <- as.vector(assays(sce)[["mean"]])
idx <- mu > 3
df <- data.frame(log10mean=log10(mu[idx]),
                 log10var=log10(var[idx]))
```


```r
library(ggplot2)
ggplot(df, aes(log10mean, log10var)) +
  geom_hex(bins=100) + 
  geom_abline(intercept=0, slope=1, col="red")
```

<img src="alevin2bioc_files/figure-html/var-mean-1.png" width="480" />


```r
rratio <- rowMeans(assays(sce)[["variance"]] /
                   (assays(sce)[["mean"]] + 1))
rmu <- rowMeans(assays(sce)[["mean"]])
idx <- rmu > 3
df <- data.frame(log10mean=log10(rmu[idx]),
                 rratio=rratio[idx],
                 gene=mcols(sce)$SYMBOL[idx])
```


```r
with(df, plot(log10mean, rratio))
high.uncert <- c(187,188)
with(df[high.uncert,],
     points(log10mean, rratio, pch=20, col="red"))
with(df[high.uncert,],
     text(log10mean, rratio, gene, pos=4))
```

<img src="alevin2bioc_files/figure-html/high-var-mean-ratio-1.png" width="480" />

```r
mcols(sce)$SYMBOL[idx][c(668,669)]
```

```
## <NA> <NA> 
##   NA   NA
```

Read-through of nearby gene `ENSG00000265681.7` (RPL17):


```r
plotInfReps(sce[,1:100], idx="ENSG00000215472.10",
            x="cluster", mainCol="SYMBOL")
```

<img src="alevin2bioc_files/figure-html/plot-high-uncert-rd-1.png" width="480" />


```r
plotInfReps(sce[,1:100], idx="ENSG00000265681.7",
            x="cluster", mainCol="SYMBOL")
```

<img src="alevin2bioc_files/figure-html/plot-high-uncert-gene-1.png" width="480" />

## Downstream analysis with Seurat

We now load *Seurat* [@seurat] ...


```r
library(Seurat)
cts <- assays(sce)[["counts"]]
pbmc <- CreateSeuratObject(cts)
```


```r
mt.genes <- rownames(sce)[as.logical(seqnames(sce) == "chrM")]
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features=mt.genes)
feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
VlnPlot(pbmc, features=feats, ncol=3)
```

<img src="alevin2bioc_files/figure-html/seurat-violin-1.png" width="480" />

## Links for further reading

## References
