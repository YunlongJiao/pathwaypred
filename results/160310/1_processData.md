---
title: "Process ICGC BRCA data"
author: "Yunlong Jiao"
date: "10 March 2016"
output: html_document
---

This script processes data and generates matrices for naive models later to make predictions.


```r
knitr::opts_chunk$set(error = FALSE)
set.seed(35875954)
source("../../src/func.R")
datapath <- '../../data/BRCA_saved_data/'
```

# Features

Feature matrices are provided by Marta stored at `crom01:/fsclinic/common/analisis/data/BRCA_saved_data`. Note that all feature matrices have to be initially provided by features in rows and patients in cols (to facilitate processing in a universal style). These matrix objects should have a name ending with `[.]vals$`.


```r
# gene expression
load(paste0(datapath, '110_genes_vals.RData'))

# gene expression wrt those present in pathways
load(paste0(datapath, '110_mini_genes_vals.RData'))

# pathway-wise features
load(paste0(datapath, '110_hipathia_matrices.RData'))

xlist <- ls(pattern = '[.]vals$')
# show all feature matrices
xlist
```

```
## [1] "eff.vals"        "fun.vals"        "genes.vals"      "go.vals"        
## [5] "mini.genes.vals" "path.vals"
```

```r
for (xname in xlist) { # process each feature matrices
  # check for no duplicated row- or col- names
  stopifnot(!any(duplicated(rownames(get(xname)))))
  stopifnot(!any(duplicated(colnames(get(xname)))))
  stopifnot(!any(is.na(get(xname))))
  assign(xname, t(get(xname))) # make patients in rows and features in cols
  cat('-------------> ', xname, ' <------------- \n')
  print(str(get(xname)))
  cat('\n')
}
```

```
## ------------->  eff.vals  <------------- 
##  num [1:881, 1:1038] 0.00407 0.0044 0.00438 0.00426 0.00426 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:1038] "hsa04014__42" "hsa04014__43" "hsa04014__44" "hsa04014__33" ...
## NULL
## 
## ------------->  fun.vals  <------------- 
##  num [1:881, 1:81] 0.09 0.0865 0.0904 0.0936 0.0929 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:81] "Lipid degradation" "Lipid metabolism" "Transcription regulation" "Apoptosis" ...
## NULL
## 
## ------------->  genes.vals  <------------- 
##  num [1:881, 1:18708] 0.417 0.418 0.418 0.426 0.401 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:18708] "1" "29974" "87769" "2" ...
## NULL
## 
## ------------->  go.vals  <------------- 
##  num [1:881, 1:370] 0.0247 0.0262 0.0259 0.0259 0.0265 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:370] "glycerophospholipid catabolic process" "phospholipid metabolic process" "multicellular organismal lipid catabolic process" "positive regulation of phospholipase activity" ...
## NULL
## 
## ------------->  mini.genes.vals  <------------- 
##  num [1:881, 1:2212] 0.513 0.529 0.525 0.505 0.532 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:2212] "5594" "5595" "5604" "5605" ...
## NULL
## 
## ------------->  path.vals  <------------- 
##  num [1:881, 1:6101] 0.00209 0.0022 0.00205 0.00186 0.00218 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:6101] "hsa04014__14 - 42" "hsa04014__14 - 43" "hsa04014__14 - 44" "hsa04014__14 - 33" ...
## NULL
```

```r
# set samples
samplelist <- unique(lapply(xlist, function(xname) rownames(get(xname))))
stopifnot(length(samplelist) == 1)
samplelist <- unlist(samplelist)
nsample <- length(samplelist)
```

# Groups

Binary classes are created from clinical info for classification. There vector objects are named ending with `[.]grps$`.


```r
#### TODO
rbinom.grps <- factor(rbinom(nsample, 1, 0.5), labels = c('neg', 'pos'))
names(rbinom.grps) <- samplelist

ylist <- ls(pattern = '[.]grps$')
# show all label vectors
ylist
```

```
## [1] "rbinom.grps"
```

# Predictors

See `../../src/func.R` for all predictors implemented. These function objects should all have a name starting with `^predictor`.


```r
prlist <- ls(pattern = '^predictor')
prlist
```

```
##  [1] "predictorDT"         "predictorGBM"        "predictorKNN"       
##  [4] "predictorLDA"        "predictorLogit"      "predictorLogitLasso"
##  [7] "predictorNNet"       "predictorRF"         "predictorSparseSVM" 
## [10] "predictorSVM"
```

# Save up !!


```r
# write out combinations of datasets for running on cluster
param <- expand.grid(xlist, ylist, prlist, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
write.table(param, file = '2_cluster_param.txt', quote = FALSE, row.names = TRUE, col.names = FALSE)

# save entire image to be loaded later
if (!dir.exists('Robj')) dir.create('Robj')
save(list = c(xlist, ylist), file = 'Robj/dat.RData')
```

# Session info


```r
sessionInfo()
```

```
## R version 3.2.1 (2015-06-18)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] knitr_1.12.3
## 
## loaded via a namespace (and not attached):
## [1] magrittr_1.5   formatR_1.3    tools_3.2.1    stringi_1.0-1 
## [5] methods_3.2.1  stringr_1.0.0  evaluate_0.8.3
```
