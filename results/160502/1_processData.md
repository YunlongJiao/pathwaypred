# Process ICGC BRCA data
# Yunlong Jiao, 02 May 2016

This script processes data and generates matrices for naive models later to make predictions.


```r
knitr::opts_chunk$set(error = FALSE)
set.seed(35875954)
source("../../src/func.R")
datapath <- '../../data/BRCA_saved_data/'
```

## Features

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
  # preview
  cat('\n-------------> \t preview \t <-------------\n')
  cat('\n-------------> \t ', xname, ' \t <-------------\n')
  print(str(get(xname)))
  cat('\n')
}
```

```
## 
## -------------> 	 preview 	 <-------------
## 
## -------------> 	  eff.vals  	 <-------------
##  num [1:881, 1:1038] 0.00407 0.0044 0.00438 0.00426 0.00426 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:1038] "X_hsa04014__42" "X_hsa04014__43" "X_hsa04014__44" "X_hsa04014__33" ...
## NULL
## 
## 
## -------------> 	 preview 	 <-------------
## 
## -------------> 	  fun.vals  	 <-------------
##  num [1:881, 1:81] 0.09 0.0865 0.0904 0.0936 0.0929 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:81] "X_Lipid_degradation" "X_Lipid_metabolism" "X_Transcription_regulation" "X_Apoptosis" ...
## NULL
## 
## 
## -------------> 	 preview 	 <-------------
## 
## -------------> 	  genes.vals  	 <-------------
##  num [1:881, 1:18708] 0.417 0.418 0.418 0.426 0.401 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:18708] "X_1" "X_29974" "X_87769" "X_2" ...
## NULL
## 
## 
## -------------> 	 preview 	 <-------------
## 
## -------------> 	  go.vals  	 <-------------
##  num [1:881, 1:370] 0.0247 0.0262 0.0259 0.0259 0.0265 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:370] "X_glycerophospholipid_catabolic_process" "X_phospholipid_metabolic_process" "X_multicellular_organismal_lipid_catabolic_process" "X_positive_regulation_of_phospholipase_activity" ...
## NULL
## 
## 
## -------------> 	 preview 	 <-------------
## 
## -------------> 	  mini.genes.vals  	 <-------------
##  num [1:881, 1:2212] 0.513 0.529 0.525 0.505 0.532 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:2212] "X_5594" "X_5595" "X_5604" "X_5605" ...
## NULL
## 
## 
## -------------> 	 preview 	 <-------------
## 
## -------------> 	  other.genes.vals  	 <-------------
##  num [1:881, 1:16496] 0.417 0.418 0.418 0.426 0.401 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:881] "TCGA.BH.A0W3.01A.11R.A109.07" "TCGA.BH.A0W4.01A.11R.A109.07" "TCGA.BH.A0DX.01A.11R.A115.07" "TCGA.BH.A0W7.01A.11R.A115.07" ...
##   ..$ : chr [1:16496] "X_1" "X_29974" "X_87769" "X_2" ...
## NULL
## 
## 
## -------------> 	 preview 	 <-------------
## 
## -------------> 	  path.vals  	 <-------------
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
# show all feature matrices
xlist
```

```
## [1] "eff.vals"         "fun.vals"         "genes.vals"      
## [4] "go.vals"          "mini.genes.vals"  "other.genes.vals"
## [7] "path.vals"
```

Now we pull different feature types and combine them to a big feature matrix. To avoid high degree of co-linearity, we choose finally only (decomposed) pathway values and gene values together as the other feature matrices are designed simply by combining pathway values linearly.


```r
xlist2combn <- c("path.vals", "mini.genes.vals", "other.genes.vals")
x <- mget(xlist2combn)
for (xname in xlist2combn) {
  colnames(x[[xname]]) <- paste(xname, colnames(x[[xname]]), sep = "_")
}
x <- do.call('cbind', x)
assign(paste0(xlist2combn, collapse = "_"), x)
xlist <- c(xlist2combn, "genes.vals", 
           paste0(xlist2combn, collapse = "_"))
# show all feature matrices NEWLY created by combining different feature types
xlist
```

```
## [1] "path.vals"                                 
## [2] "mini.genes.vals"                           
## [3] "other.genes.vals"                          
## [4] "genes.vals"                                
## [5] "path.vals_mini.genes.vals_other.genes.vals"
```

## Groups (binary and multi-class)

Binary classes are created from clinical info for classification. There vector objects are named ending with `[.]grps$`.

#### Basal vs Non-basal


```r
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
# contrasting class size
table(basal.grps[sl])
```

```
## 
## neg_NonBasal    pos_Basal 
##          401           94
```

#### Tumor vs Normal


```r
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
# contrasting class size
table(tumor.grps[sl])
```

```
## 
## neg_Normal  pos_Tumor 
##        102        779
```

#### Alive vs Deceased (overall survival)


```r
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
## -------------> 	  basal.grps  	 <-------------
## [1] "neg_NonBasal" "pos_Basal"   
## 
## 
## -------------> 	 summary 	 <-------------
## 
## -------------> 	  subtype.grps  	 <-------------
## [1] "T1_Luminal" "T2_Her2"    "T3_Basal"  
## 
## 
## -------------> 	 summary 	 <-------------
## 
## -------------> 	  surv.grps  	 <-------------
## [1] "neg_alive"    "pos_deceased"
## 
## 
## -------------> 	 summary 	 <-------------
## 
## -------------> 	  tumor.grps  	 <-------------
## [1] "neg_Normal" "pos_Tumor"
```

## Predictors

See `../../src/func.R` for all predictors implemented. These function objects should all have a name starting with `^predictor`.

**NOTE we are only interested in predictors that are capable of FS in this study.**


```r
prlist <- ls(pattern = '^predictor.')
prlist <- grep("(lasso|pam|rf)", prlist, ignore.case = TRUE, value = TRUE)
prlist
```

```
## [1] "predictorLogitLasso" "predictorPAM"        "predictorRF"
```

## (Nested) cross validation parameters

**NOTE nested CV for tuning predictors is not necessary in this study but we keep it for ease of unified syntax.**


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
                     as.integer(0), # inner CV folds for tuning predictor NOT NECESSARY
                     as.integer(nfolds.inn), 
                     as.integer(nrepeats.inn), 
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
write.table(param, file = '2_runPredict.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = ' ')
# preview
str(param)
```

```
## 'data.frame':	3000 obs. of  9 variables:
##  $ Var1: chr  "path.vals" "mini.genes.vals" "other.genes.vals" "genes.vals" ...
##  $ Var2: chr  "basal.grps" "basal.grps" "basal.grps" "basal.grps" ...
##  $ Var3: chr  "predictorLogitLasso" "predictorLogitLasso" "predictorLogitLasso" "predictorLogitLasso" ...
##  $ Var4: int  1 1 1 1 1 1 1 1 1 1 ...
##  $ Var5: int  5 5 5 5 5 5 5 5 5 5 ...
##  $ Var6: int  10 10 10 10 10 10 10 10 10 10 ...
##  $ Var7: int  0 0 0 0 0 0 0 0 0 0 ...
##  $ Var8: int  5 5 5 5 5 5 5 5 5 5 ...
##  $ Var9: int  1 1 1 1 1 1 1 1 1 1 ...
```

```r
# save entire image to be loaded later
save(list = c(xlist, ylist), file = 'dat.RData')
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
## [1] magrittr_1.5   formatR_1.3    tools_3.2.1    stringi_1.0-1 
## [5] knitr_1.12.3   stringr_1.0.0  evaluate_0.8.3
```
