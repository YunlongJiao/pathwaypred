# Process ICGC BRCA data
# Yunlong Jiao, 10 March 2016

This script processes data and generates matrices for naive models later to make predictions.


```r
knitr::opts_chunk$set(error = FALSE)
set.seed(35875954)
source("../../src/func.R")
datapath <- '../../data/BRCA_saved_data/'
```

# Features

Feature matrices are stored at `crom01:/fsclinic/common/analisis/data/BRCA_saved_data/`. Note that all feature matrices have to be initially provided by features in rows and patients in cols (to facilitate processing in a universal style). These matrix objects should have a name ending with `[.]vals$`.


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

# show all feature matrices
xlist <- ls(pattern = '[.]vals$')
xlist
```

```
## [1] "eff.vals"         "fun.vals"         "genes.vals"      
## [4] "go.vals"          "mini.genes.vals"  "other.genes.vals"
## [7] "path.vals"
```

```r
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
  # preview
  cat('\n-------------> ', xname, ' <-------------\n')
  print(str(get(xname)))
  cat('\n')
}
```

```
## 
## ------------->  eff.vals  <-------------
##  num [1:881, 1:1038] 0.00407 0.0044 0.00438 0.00426 0.00426 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:1038] "X_hsa04014__42" "X_hsa04014__43" "X_hsa04014__44" "X_hsa04014__33" ...
## NULL
## 
## 
## ------------->  fun.vals  <-------------
##  num [1:881, 1:81] 0.09 0.0865 0.0904 0.0936 0.0929 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:81] "X_Lipid_degradation" "X_Lipid_metabolism" "X_Transcription_regulation" "X_Apoptosis" ...
## NULL
## 
## 
## ------------->  genes.vals  <-------------
##  num [1:881, 1:18708] 0.417 0.418 0.418 0.426 0.401 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:18708] "X_1" "X_29974" "X_87769" "X_2" ...
## NULL
## 
## 
## ------------->  go.vals  <-------------
##  num [1:881, 1:370] 0.0247 0.0262 0.0259 0.0259 0.0265 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:370] "X_glycerophospholipid_catabolic_process" "X_phospholipid_metabolic_process" "X_multicellular_organismal_lipid_catabolic_process" "X_positive_regulation_of_phospholipase_activity" ...
## NULL
## 
## 
## ------------->  mini.genes.vals  <-------------
##  num [1:881, 1:2212] 0.513 0.529 0.525 0.505 0.532 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:2212] "X_5594" "X_5595" "X_5604" "X_5605" ...
## NULL
## 
## 
## ------------->  other.genes.vals  <-------------
##  num [1:881, 1:16496] 0.417 0.418 0.418 0.426 0.401 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:16496] "X_1" "X_29974" "X_87769" "X_2" ...
## NULL
## 
## 
## ------------->  path.vals  <-------------
##  num [1:881, 1:6101] 0.00209 0.0022 0.00205 0.00186 0.00218 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:6101] "X_hsa04014__14___42" "X_hsa04014__14___43" "X_hsa04014__14___44" "X_hsa04014__14___33" ...
## NULL
```

```r
# set samples
samplelist <- unique(lapply(xlist, function(xname) rownames(get(xname))))
stopifnot(length(samplelist) == 1)
samplelist <- unlist(samplelist)
nsample <- length(samplelist)
```

# Groups (binary)

Binary classes are created from clinical info for classification. There vector objects are named ending with `[.]grps$`.


```r
# basal type
y <- read.table(paste0(datapath, 'sampleBasalType.txt'), header = FALSE)
head(y$V2)
```

```
## [1] "Basal" "Basal" "Basal" "Basal" "Basal" "Basal"
```

```r
basal.grps <- ordered(gsub('[^[:alnum:]_]', '_', y$V2), levels = c("Not_Basal", "Basal"), labels = c("neg_NonBasal", "pos_Basal"))
stopifnot(!any(is.na(basal.grps)))
names(basal.grps) <- y$V1
head(basal.grps)
```

```
## TCGA.A2.A0T2.01A.11R.A084.07 TCGA.A1.A0SK.01A.12R.A084.07 
##                    pos_Basal                    pos_Basal 
## TCGA.A2.A0CM.01A.31R.A034.07 TCGA.AR.A1AR.01A.31R.A137.07 
##                    pos_Basal                    pos_Basal 
## TCGA.BH.A18V.11A.52R.A12D.07 TCGA.BH.A18V.01A.11R.A12D.07 
##                    pos_Basal                    pos_Basal 
## Levels: neg_NonBasal < pos_Basal
```

```r
# sample size (how many labels are available in the sample cohort)
length(sl <- intersect(samplelist, names(basal.grps)))
```

```
## [1] 495
```

```r
table(basal.grps[sl])
```

```
## 
## neg_NonBasal    pos_Basal 
##          401           94
```

```r
# phenotype tumor
y <- read.table(paste0(datapath, 'TCGA_phenotype.txt'), header = TRUE)
head(y$cancer)
```

```
## [1] "Tumor" "Tumor" "Tumor" "Tumor" "Tumor" "Tumor"
```

```r
tumor.grps <- ordered(gsub('[^[:alnum:]_]', '_', y$cancer), levels = c("Normal", "Tumor"), labels = c("neg_Normal", "pos_Tumor"))
stopifnot(!any(is.na(tumor.grps)))
names(tumor.grps) <- rownames(y)
head(tumor.grps)
```

```
## TCGA.C4.A0F0.01A.12R.A10U.07 TCGA.C4.A0F6.01A.11R.A10U.07 
##                    pos_Tumor                    pos_Tumor 
## TCGA.BL.A0C8.01A.11R.A10U.07 TCGA.BL.A13J.01A.11R.A10U.07 
##                    pos_Tumor                    pos_Tumor 
## TCGA.BT.A20J.01A.11R.A14Y.07 TCGA.BT.A20N.01A.11R.A14Y.07 
##                    pos_Tumor                    pos_Tumor 
## Levels: neg_Normal < pos_Tumor
```

```r
# sample size (how many labels are available in the sample cohort)
length(sl <- intersect(samplelist, names(tumor.grps)))
```

```
## [1] 881
```

```r
table(tumor.grps[sl])
```

```
## 
## neg_Normal  pos_Tumor 
##        102        779
```

```r
# survival
alignmat <- read.table(paste0(datapath, 'clinical_info_sample.txt'), sep = "\t", quote = "\"", header = TRUE)
y <- read.table(paste0(datapath, 'clinical_info_donor.txt'), sep = "\t", quote = "\"", header = TRUE)
alignid <- data.frame(sample = samplelist, 
                      donor = sapply(samplelist, function(id){
                        message(".", appendLF = FALSE)
                        i <- grep(id, alignmat$submitted_sample_id)
                        if (length(i) == 1) {
                          return(alignmat$icgc_donor_id[i])
                        } else if (length(i) == 0) {
                          return(NA)
                        } else {
                          stop(id, " has more than 1 alignment!")
                        }
                      }))
# number of samples before pruning
nrow(alignid)
```

```
## [1] 881
```

```r
alignid <- subset(alignid, !is.na(alignid$donor))
id <- match(alignid$donor, y$icgc_donor_id)
alignid <- alignid[!is.na(id), ]
id <- id[!is.na(id)]
y <- y[id, ]
# number of samples after pruning
nrow(alignid)
```

```
## [1] 879
```

```r
head(y$donor_vital_status)
```

```
## [1] "alive"    "alive"    "alive"    "alive"    "deceased" "deceased"
```

```r
surv.grps <- ordered(gsub('[^[:alnum:]_]', '_', y$donor_vital_status), levels = c("alive", "deceased"), labels = c("neg_alive", "pos_deceased"))
stopifnot(!any(is.na(surv.grps)))
names(surv.grps) <- alignid$sample
head(surv.grps)
```

```
## TCGA.BH.A0W3.01A.11R.A109.07 TCGA.BH.A0W4.01A.11R.A109.07 
##                    neg_alive                    neg_alive 
## TCGA.BH.A0DX.01A.11R.A115.07 TCGA.BH.A0W7.01A.11R.A115.07 
##                    neg_alive                    neg_alive 
## TCGA.BH.A18U.11A.23R.A12D.07 TCGA.BH.A18P.01A.11R.A12D.07 
##                 pos_deceased                 pos_deceased 
## Levels: neg_alive < pos_deceased
```

```r
# sample size (how many labels are available in the sample cohort)
length(sl <- intersect(samplelist, names(surv.grps)))
```

```
## [1] 879
```

```r
table(surv.grps[sl])
```

```
## 
##    neg_alive pos_deceased 
##          755          124
```

```r
# show all label vectors
ylist <- ls(pattern = '[.]grps$')
ylist
```

```
## [1] "basal.grps" "surv.grps"  "tumor.grps"
```

```r
# check for now to deal with binary classification only
for (yname in ylist) {
  stopifnot(length(levels(get(yname))) == 2)
}
```

# Predictors

See `../../src/func.R` for all predictors implemented. These function objects should all have a name starting with `^predictor`.


```r
prlist <- ls(pattern = '^predictor.')
prlist
```

```
## [1] "predictorGBM"        "predictorKNN"        "predictorLDA"       
## [4] "predictorLinearSVM"  "predictorLogitLasso" "predictorNB"        
## [7] "predictorRF"         "predictorRadialSVM"  "predictorSparseSVM"
```

# Cross validation parameters


```r
# CV folds for evaluation
nfolds <- 5
nrepeats <- 10
# inner CV folds for tuning predictor
nfolds.inn <- 5
nrepeats.inn <- 1
```

# Save up !!


```r
# write out combinations of datasets for running on cluster
param <- expand.grid(xlist, ylist, prlist, # feature, label, predictor
                     seq(nfolds * nrepeats), nfolds, nrepeats, # CV folds for evaluation
                     seq(nfolds.inn * nrepeats.inn), nfolds.inn, nrepeats.inn, # inner CV folds for tuning predictor
                     KEEP.OUT.ATTRS = FALSE)
write.table(param, file = 'cluster_param.txt', quote = FALSE, row.names = TRUE, col.names = FALSE, sep = ' ')

# save entire image to be loaded later
save(list = c(xlist, ylist), file = 'dat.RData')
```

# Session info


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-apple-darwin13.4.0 (64-bit)
## Running under: OS X 10.11.4 (El Capitan)
## 
## locale:
## [1] C/UTF-8/C/C/C/C
## 
## attached base packages:
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] knitr_1.12.3
## 
## loaded via a namespace (and not attached):
## [1] magrittr_1.5  formatR_1.2.1 tools_3.2.3   stringi_1.0-1 stringr_1.0.0
## [6] evaluate_0.8
```
