# Process ICGC BRCA data
# Yunlong Jiao, 21 April 2016

This script processes data and generates matrices for naive models later to make predictions.


```r
knitr::opts_chunk$set(error = FALSE)
set.seed(35875954)
source("../../src/func.R")
datapath <- '../../data/BRCA_saved_data/'
```

## Features

Feature matrices are stored at ../../data/BRCA_saved_data/. Note that all feature matrices have to be initially provided by features in rows and patients in cols (to facilitate processing in a universal style). These matrix objects should have a name ending with `[.]vals$`.


```r
# gene expression
load(paste0(datapath, '110_genes_vals.RData'))

# gene expression wrt those present in pathways
load(paste0(datapath, '110_mini_genes_vals.RData'))

# pathway-wise features
load(paste0(datapath, '110_hipathia_matrices.RData'))

# complement of mini.genes from all.genes
stopifnot(rownames(mini.genes.vals) %in% rownames(genes.vals))
other.genes.vals <- genes.vals[setdiff(rownames(genes.vals), rownames(mini.genes.vals)), ]

# get all feature matrices
xlist <- ls(pattern = '[.]vals$')
for (xname in xlist) {
  # process each feature matrices
  x <- get(xname)
  # check no NA values
  stopifnot(!any(is.na(x)))
  # make patients in rows and features in cols
  x <- t(x)
  # reform feature names and check for duplicated names
  colnames(x) <- paste0('X_', gsub('[^[:alnum:]_]', '_', colnames(x)))
  stopifnot(!any(duplicated(rownames(x))))
  stopifnot(!any(duplicated(colnames(x))))
  # assign back
  assign(xname, x)
}
# set samples
samplelist <- unique(lapply(xlist, function(xname) rownames(get(xname))))
stopifnot(length(samplelist) == 1)
samplelist <- unlist(samplelist)
nsample <- length(samplelist)

# create mixed feature matrices
xlist2combn <- list("eff.and.other.genes.vals" = c("eff.vals", "other.genes.vals"), 
                    "path.and.other.genes.vals" = c("path.vals", "other.genes.vals"), 
                    "path.and.genes.vals" = c("path.vals", "genes.vals"))
for (xname in names(xlist2combn)) {
  x <- mget(xlist2combn[[xname]])
  x <- do.call("cbind", x)
  assign(xname, x)
}
xlist <- c(xlist, names(xlist2combn))

# show all feature matrices
xlist
```

```
##  [1] "eff.vals"                  "fun.vals"                 
##  [3] "genes.vals"                "go.vals"                  
##  [5] "mini.genes.vals"           "other.genes.vals"         
##  [7] "path.vals"                 "eff.and.other.genes.vals" 
##  [9] "path.and.other.genes.vals" "path.and.genes.vals"
```

```r
cat('\n-------------> \t preview \t <-------------\n')
```

```
## 
## -------------> 	 preview 	 <-------------
```

```r
for (xname in xlist) {
  cat('\n-------------> \t ', xname, ' \t <-------------\n')
  print(str(get(xname)))
  cat('\n')
}
```

```
## 
## -------------> 	  eff.vals  	 <-------------
##  num [1:881, 1:1038] 0.00407 0.0044 0.00438 0.00426 0.00426 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:1038] "X_hsa04014__42" "X_hsa04014__43" "X_hsa04014__44" "X_hsa04014__33" ...
## NULL
## 
## 
## -------------> 	  fun.vals  	 <-------------
##  num [1:881, 1:81] 0.09 0.0865 0.0904 0.0936 0.0929 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:81] "X_Lipid_degradation" "X_Lipid_metabolism" "X_Transcription_regulation" "X_Apoptosis" ...
## NULL
## 
## 
## -------------> 	  genes.vals  	 <-------------
##  num [1:881, 1:18708] 0.417 0.418 0.418 0.426 0.401 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:18708] "X_1" "X_29974" "X_87769" "X_2" ...
## NULL
## 
## 
## -------------> 	  go.vals  	 <-------------
##  num [1:881, 1:370] 0.0247 0.0262 0.0259 0.0259 0.0265 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:370] "X_glycerophospholipid_catabolic_process" "X_phospholipid_metabolic_process" "X_multicellular_organismal_lipid_catabolic_process" "X_positive_regulation_of_phospholipase_activity" ...
## NULL
## 
## 
## -------------> 	  mini.genes.vals  	 <-------------
##  num [1:881, 1:2212] 0.513 0.529 0.525 0.505 0.532 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:2212] "X_5594" "X_5595" "X_5604" "X_5605" ...
## NULL
## 
## 
## -------------> 	  other.genes.vals  	 <-------------
##  num [1:881, 1:16496] 0.417 0.418 0.418 0.426 0.401 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:16496] "X_1" "X_29974" "X_87769" "X_2" ...
## NULL
## 
## 
## -------------> 	  path.vals  	 <-------------
##  num [1:881, 1:6101] 0.00209 0.0022 0.00205 0.00186 0.00218 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:6101] "X_hsa04014__14___42" "X_hsa04014__14___43" "X_hsa04014__14___44" "X_hsa04014__14___33" ...
## NULL
## 
## 
## -------------> 	  eff.and.other.genes.vals  	 <-------------
##  num [1:881, 1:17534] 0.00407 0.0044 0.00438 0.00426 0.00426 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:17534] "X_hsa04014__42" "X_hsa04014__43" "X_hsa04014__44" "X_hsa04014__33" ...
## NULL
## 
## 
## -------------> 	  path.and.other.genes.vals  	 <-------------
##  num [1:881, 1:22597] 0.00209 0.0022 0.00205 0.00186 0.00218 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:22597] "X_hsa04014__14___42" "X_hsa04014__14___43" "X_hsa04014__14___44" "X_hsa04014__14___33" ...
## NULL
## 
## 
## -------------> 	  path.and.genes.vals  	 <-------------
##  num [1:881, 1:24809] 0.00209 0.0022 0.00205 0.00186 0.00218 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:24809] "X_hsa04014__14___42" "X_hsa04014__14___43" "X_hsa04014__14___44" "X_hsa04014__14___33" ...
## NULL
```

## Groups (multi-class only)

Multiple classes are created from clinical info for classification. There vector objects are named ending with `[.]grps$`.

**NOTE always put THE positive group to the last level in the ordered factor as the evaluation measure such as TPR, FPR will be specific to that level.**

#### Luminal vs Her2 vs Basal (subtypes)


```r
y <- read.table(paste0(datapath, 'sampleCancerType.txt'), header = FALSE)
head(y$V2)
```

```
## [1] "Basal" "Basal" "Basal" "Basal" "Basal" "Basal"
```

```r
unique(y$V2)
```

```
## [1] "Basal"   "Luminal" "Unknown" "Her2"
```

```r
subtype.grps <- ordered(gsub('[^[:alnum:]_]', '_', y$V2), levels = c("Luminal", "Her2", "Basal"), labels = c("T1_Luminal", "T2_Her2", "T3_Basal"))
id.na <- is.na(subtype.grps)
subtype.grps <- subtype.grps[!id.na]
names(subtype.grps) <- y$V1[!id.na]
head(subtype.grps)
```

```
## TCGA.A2.A0T2.01A.11R.A084.07 TCGA.A1.A0SK.01A.12R.A084.07 
##                     T3_Basal                     T3_Basal 
## TCGA.A2.A0CM.01A.31R.A034.07 TCGA.AR.A1AR.01A.31R.A137.07 
##                     T3_Basal                     T3_Basal 
## TCGA.BH.A18V.11A.52R.A12D.07 TCGA.BH.A18V.01A.11R.A12D.07 
##                     T3_Basal                     T3_Basal 
## Levels: T1_Luminal < T2_Her2 < T3_Basal
```

```r
# sample size (how many labels are available in the sample cohort)
length(sl <- intersect(samplelist, names(subtype.grps)))
```

```
## [1] 495
```

```r
# contrasting class size
table(subtype.grps[sl])
```

```
## 
## T1_Luminal    T2_Her2   T3_Basal 
##        315         86         94
```

#### Summary


```r
# get all groups
ylist <- ls(pattern = '[.]grps$')
# show all groups and their classes
for (yname in ylist) {
  cat('\n-------------> \t summary \t <-------------\n')
  cat('\n-------------> \t ', yname, ' \t <-------------\n')
  print(levels(get(yname)))
  cat('\n')
}
```

```
## 
## -------------> 	 summary 	 <-------------
## 
## -------------> 	  subtype.grps  	 <-------------
## [1] "T1_Luminal" "T2_Her2"    "T3_Basal"
```

## Predictors

See `../../src/func.R` for all predictors implemented. These function objects should all have a name starting with `^predictor`.


```r
prlist <- ls(pattern = '^predictor.')
prlist
```

```
##  [1] "predictorConstant"   "predictorGBM"        "predictorKendallSVM"
##  [4] "predictorKNN"        "predictorLDA"        "predictorLinearSVM" 
##  [7] "predictorLogitLasso" "predictorNB"         "predictorPAM"       
## [10] "predictorRadialSVM"  "predictorRF"         "predictorSparseSVM"
```

## Kernel matrices

This section computes kernel matrix to accelerate implementation of `kernlab::ksvm`. These kernel matrices are save at ../../data/BRCA_saved_data/ and named ending with `[.]kmat$`. Specifically we have the following kernel
- Kendall kernel (Jiao and Vert, 2016) for predictorKendallSVM

**NOTE before we compute kernel here **


```r
# kendall kernel
prname <- "predictorKendallSVM"
kernel <- pcaPP::cor.fk
for (xname in xlist) {
  kmatname <- paste(prname, xname, "kmat", sep = ".")
  message(kmatname)
  kmatpath <- paste0(datapath, kmatname, ".RData")
  if (file.exists(kmatpath)) {
    assign(kmatname, get(load(kmatpath)))
  } else {
    x <- get(xname)
    x <- removeConst(x)
    assign(kmatname, computeKernelMatrix(x = x, kernel = kernel))
    save(list = kmatname, file = kmatpath)
  }
}
kmatlist <- ls(pattern = '[.]kmat$')
kmatlist
```

```
##  [1] "predictorKendallSVM.eff.and.other.genes.vals.kmat" 
##  [2] "predictorKendallSVM.eff.vals.kmat"                 
##  [3] "predictorKendallSVM.fun.vals.kmat"                 
##  [4] "predictorKendallSVM.genes.vals.kmat"               
##  [5] "predictorKendallSVM.go.vals.kmat"                  
##  [6] "predictorKendallSVM.mini.genes.vals.kmat"          
##  [7] "predictorKendallSVM.other.genes.vals.kmat"         
##  [8] "predictorKendallSVM.path.and.genes.vals.kmat"      
##  [9] "predictorKendallSVM.path.and.other.genes.vals.kmat"
## [10] "predictorKendallSVM.path.vals.kmat"
```

## (Nested) cross validation parameters


```r
# CV folds for evaluation
nfolds <- 5
nrepeats <- 10
# inner CV folds for tuning predictor
nfolds.inn <- 5
nrepeats.inn <- 1
```

## Save up !!


```r
# for cluster jobs, write out full combinations of double nested CV parameters
# - for outter loop, 1:(nfolds * nrepeats) denotes index for outter CV
# - for inner loop, 0 indicates running outter loop evaluation, 1:(nfolds.inn * nrepeats.inn) indicates index for inner CV
param <- expand.grid(as.character(xlist), # feature
                     as.character(ylist), # group
                     as.character(prlist), # predictor
                     as.integer(1:(nfolds * nrepeats)), # outter CV folds index for evaluation
                     as.integer(nfolds), 
                     as.integer(nrepeats), 
                     as.integer(0:(nfolds.inn * nrepeats.inn)), # inner CV folds for tuning predictor
                     as.integer(nfolds.inn), 
                     as.integer(nrepeats.inn), 
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
write.table(param, file = '2_runPredict.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ' ')
# preview
str(param)
```

```
## 'data.frame':	36000 obs. of  9 variables:
##  $ Var1: chr  "eff.vals" "fun.vals" "genes.vals" "go.vals" ...
##  $ Var2: chr  "subtype.grps" "subtype.grps" "subtype.grps" "subtype.grps" ...
##  $ Var3: chr  "predictorConstant" "predictorConstant" "predictorConstant" "predictorConstant" ...
##  $ Var4: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ Var5: int  5 5 5 5 5 5 5 5 5 5 ...
##  $ Var6: int  10 10 10 10 10 10 10 10 10 10 ...
##  $ Var7: int  0 0 0 0 0 0 0 0 0 0 ...
##  $ Var8: int  5 5 5 5 5 5 5 5 5 5 ...
##  $ Var9: int  1 1 1 1 1 1 1 1 1 1 ...
```

```r
# save entire image to be loaded later
save(list = c(xlist, ylist, kmatlist), file = 'dat.RData')
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
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## loaded via a namespace (and not attached):
## [1] magrittr_1.5   formatR_1.3    tools_3.2.1    mvtnorm_1.0-3 
## [5] stringi_1.0-1  pcaPP_1.9-60   knitr_1.12.3   stringr_1.0.0 
## [9] evaluate_0.8.3
```
