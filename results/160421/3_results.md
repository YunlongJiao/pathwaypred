# Results
# Yunlong Jiao, 25 Apr 2016

This script collects results from earlier run and illustrates with tables and plots.


```r
knitr::opts_chunk$set(error = FALSE, fig.width = 12, fig.height = 8, dev = "pdf", fig.keep = "high", fig.path = "result_figure/", cache.path = "result_cache/")
set.seed(70236562)
source("../../src/func.R")
library(ggplot2)
```

First read in parameters!


```r
# read in parameters
param <- read.table("2_runPredict.txt", header = FALSE, row.names = NULL, col.names = c("xname", "yname", "prname", "i.fold", "nfolds", "nrepeats", "i.fold.inn", "nfolds.inn", "nrepeats.inn"))
# features
xlist <- unique(param$xname)
xlist
```

```
## [1] "eff.vals"         "fun.vals"         "genes.vals"      
## [4] "go.vals"          "mini.genes.vals"  "other.genes.vals"
## [7] "path.vals"
```

```r
# feature types
xlist.type <- c("func-wise", "func-wise", 
                "path-wise", "path-wise", 
                "gene-wise", "gene-wise", "gene-wise")
names(xlist.type) <- c("fun.vals", "go.vals", # functionality features
                       "eff.vals", "path.vals", # pathway features
                       "mini.genes.vals", "other.genes.vals", "genes.vals") # gene features
xlist.vline <- c(2.5, 4.5) # cut out types
stopifnot(length(setdiff(xlist, names(xlist.type))) == 0)
xlist.type
```

```
##         fun.vals          go.vals         eff.vals        path.vals 
##      "func-wise"      "func-wise"      "path-wise"      "path-wise" 
##  mini.genes.vals other.genes.vals       genes.vals 
##      "gene-wise"      "gene-wise"      "gene-wise"
```

```r
# groups
ylist <- unique(param$yname)
ylist
```

```
## [1] "subtype.grps"
```

```r
# predictors
prlist <- unique(param$prname)
prlist
```

```
##  [1] "predictorConstant"   "predictorGBM"        "predictorKendallSVM"
##  [4] "predictorKNN"        "predictorLDA"        "predictorLinearSVM" 
##  [7] "predictorLogitLasso" "predictorNB"         "predictorPAM"       
## [10] "predictorRadialSVM"  "predictorRF"         "predictorSparseSVM"
```

```r
# (outter) `nfolds`-fold CV repeated `nrepeats` times for evaluation
nfolds <- unique(param$nfolds)
stopifnot(length(nfolds) == 1)
nfolds
```

```
## [1] 5
```

```r
nrepeats <- unique(param$nrepeats)
stopifnot(length(nrepeats) == 1)
nrepeats
```

```
## [1] 10
```

```r
# inner `nfolds.inn`-fold CV repeated `nrepeats.inn` times for tuning predictors
nfolds.inn <- unique(param$nfolds.inn)
stopifnot(length(nfolds.inn) == 1)
nfolds.inn
```

```
## [1] 5
```

```r
nrepeats.inn <- unique(param$nrepeats.inn)
stopifnot(length(nrepeats.inn) == 1)
nrepeats.inn
```

```
## [1] 1
```

```r
# evaluation measures
slist <- c("acc","fpr","tpr","ppv","fval", "auroc")
slist
```

```
## [1] "acc"   "fpr"   "tpr"   "ppv"   "fval"  "auroc"
```

```r
slist.prefer.large.score <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE)
names(slist.prefer.large.score) <- slist
slist.prefer.large.score
```

```
##   acc   fpr   tpr   ppv  fval auroc 
##  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
```

## Tuning predictor

In order to avoid overfitting with the choice of prediction algorithm used, across each feature `xname` X each label group `yname` X each evaluation CV (outter) fold repeat `i.fold`, the predictor `prname` is selected by nested CV (inner) runs and the best prediction performance is reported by each evaluation measure.


```r
# gather scores
scores <- list()
for (i.fold in seq(nfolds*nrepeats)) {
  message("\n", i.fold, "-th fold out of ", nfolds*nrepeats, " folds ", appendLF = FALSE)
  for (yname in ylist) {
    for (xname in xlist) {
      message(".", appendLF = FALSE)
      # read nested cv res
      cvres <- list()
      for (prname in prlist) {
        res.files <- list.files(path = 'Robj', 
                                pattern = paste('^cvres', xname, yname, prname, 
                                                i.fold, nfolds, nrepeats, 
                                                '[[:digit:]]+', nfolds.inn, nrepeats.inn, 
                                                sep = '_'), 
                                full.names = TRUE)
        res <- lapply(res.files, function(f) get(load(f)))
        cvres[[prname]] <- crossValidationCombineResults(res)
      }
      
      # generate iv res where predictor has been selected by nested cv
      ivres <- list()
      for (sname in slist) {
        tt <- sapply(cvres, '[[', "system_time")
        # order scores by increasing system time
        ss <- sapply(cvres, '[[', sname)[order(tt, decreasing = FALSE)]
        ss <- ss[!is.na(ss)]
        # fast algorithm is preferred among those returning equal score values
        best.prname <- names(ss)[order(ss, decreasing = slist.prefer.large.score[sname])]
        # get iv res NOTE while loop is to guarantee no iv res is found for such predictor
        while (length(best.prname) > 0) {
          best.files <- list.files(path = 'Robj', 
                                   pattern = paste('^ivres', xname, yname, best.prname[1], 
                                                   i.fold, nfolds, nrepeats, 
                                                   '0+', nfolds.inn, nrepeats.inn, 
                                                   sep = '_'), 
                                   full.names = TRUE)
          if (length(best.files) == 0) {
            best.prname <- best.prname[-1]
            next
          } else if (length(best.files) == 1) {
            best.prname <- best.prname[1]
            ivres[[sname]] <- get(load(best.files))
            break
          } else {
            stop("multiple ivres files found for ", best.prname[1])
          }
        }
        # in case no ivres files found for any predictors
        stopifnot(length(best.prname) == 1)
      }
      
      # assign obj
      objname <- paste('res', xname, yname, i.fold, nfolds, nrepeats, sep = '_')
      assign(objname, ivres)
      
      # record score values
      scores[[objname]] <- data.frame(y = yname, 
                                      x = xname, 
                                      type = xlist.type[xname], 
                                      predictor = sapply(slist, function(sname) ivres[[sname]][["predictor"]]), 
                                      rep = i.fold, 
                                      score = slist, 
                                      value = sapply(slist, function(sname) ivres[[sname]][[sname]]), 
                                      row.names = NULL)
    }
  }
}
scores <- do.call('rbind', scores)
rownames(scores) <- seq(nrow(scores))
# write out scores
write.table(scores, file = "scores.txt", row.names = TRUE, col.names = TRUE, sep = '\t')
# preview
head(scores)
```

```
##              y        x      type           predictor rep score     value
## 1 subtype.grps eff.vals path-wise  predictorLinearSVM   1   acc 0.8888889
## 2 subtype.grps eff.vals path-wise   predictorConstant   1   fpr 0.0000000
## 3 subtype.grps eff.vals path-wise        predictorPAM   1   tpr 0.8500000
## 4 subtype.grps eff.vals path-wise  predictorLinearSVM   1   ppv 0.9444444
## 5 subtype.grps eff.vals path-wise predictorKendallSVM   1  fval 0.8947368
## 6 subtype.grps eff.vals path-wise        predictorKNN   1 auroc 0.9615385
```

## Performance overview

In this section, each plot corresponds to the prediction performance for one specific group of phenotypic response. In each plot, along x-axis we have different feature matrices, along y-axis we have different evaluation measures, the barplot shows the CV scores evaluated against randomly splitted CV folds.

Notes that some evaluation measures (`acc`, `ppv`, etc) do associate with a specific `cutoff` threshold for predicted probability, which is always set to constant 0.5. It might be interesting to tune the threshold in case that the class sizes in training data are unbalanced.

The evaluation measure are the following: eff.vals, fun.vals, genes.vals, go.vals, mini.genes.vals, other.genes.vals, path.vals. Note that in order to remove confusion, we equivalently rename `TPR` by `Sensitivity`, and plot `Specificity = 1 - FPR` instead of `FPR` values. As a result, for any evaluation measure present in the plot, higher values indicate better performance unanimously.


```r
# plot each grps in a separate figure
for (yname in ylist) {
  d <- subset(scores, scores$y == yname)
  # order feature by types
  d$x <- factor(d$x, levels = names(xlist.type), ordered = TRUE)
  # change score values by "fpr" to that by "spec"
  d$value[d$score == "fpr"] <- 1 - d$value[d$score == "fpr"]
  d$score[d$score == "fpr"] <- "spec"
  # change score name "tpr" to "sens"
  d$score[d$score == "tpr"] <- "sens"
  p1 <- ggplot(d, aes(x = x, y = value)) + 
    geom_boxplot(aes(fill = x), alpha = 0.8) + 
    geom_vline(xintercept = xlist.vline, color = "grey", size = 2) + 
    scale_y_continuous(limits = c(0.75, 1)) + 
    facet_wrap(~score, scales = "fixed") + 
    ggtitle(paste0("predicting ", yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
}
```

```
## Warning: Removed 49 rows containing non-finite values (stat_boxplot).
```

![plot of chunk overview](result_figure/overview-1.pdf)

## Predictor table w.r.t. ACC

We focus on multi-class classification accuracy (ACC) to evaluate performance.

Now we turn to look at which `predictor` suits well in making prediction for a specific group of phenotypic response `yname` using each feature `xname`. We report the frequency of such predictor is selected across 50 repeats of evaluation CV runs.


```r
# focus only auroc
sname <- "acc"
# show each grp separately
cat('\n---> \t Frequency (%) of predictor being selected \t <---\n')
for (yname in ylist) {
  cat('\n---> \t predicting for ', yname, ' using \t <---\n')
  d <- subset(scores, scores$score == sname & scores$y == yname)
  d <- split(d, d$x)
  freq.best.prname <- lapply(d, function(u){
    s <- sort(table(u$predictor)/length(u$predictor)*100, decreasing = TRUE)
    formatC(s, format = 'f', digits = 2)
  })
  print(freq.best.prname)
}
```

```
## 
## ---> 	 Frequency (%) of predictor being selected 	 <---
## 
## ---> 	 predicting for  subtype.grps  using 	 <---
## $eff.vals
## 
##        predictorGBM predictorLogitLasso  predictorLinearSVM 
##             "74.00"             "20.00"              "4.00" 
##  predictorSparseSVM 
##              "2.00" 
## 
## $fun.vals
## 
##  predictorLinearSVM  predictorRadialSVM        predictorLDA 
##             "66.00"             "28.00"              "4.00" 
## predictorKendallSVM 
##              "2.00" 
## 
## $genes.vals
## 
##        predictorGBM predictorLogitLasso  predictorSparseSVM 
##             "72.00"             "16.00"             "12.00" 
## 
## $go.vals
## 
## predictorKendallSVM predictorLogitLasso  predictorRadialSVM 
##             "44.00"             "26.00"             "18.00" 
##  predictorLinearSVM        predictorGBM 
##              "8.00"              "4.00" 
## 
## $mini.genes.vals
## 
##        predictorGBM predictorLogitLasso 
##             "94.00"              "6.00" 
## 
## $other.genes.vals
## 
##        predictorGBM predictorLogitLasso  predictorSparseSVM 
##             "62.00"             "22.00"             "14.00" 
##         predictorRF 
##              "2.00" 
## 
## $path.vals
## 
##        predictorGBM predictorLogitLasso         predictorRF 
##             "66.00"             "28.00"              "6.00"
```

## Significance test w.r.t. ACC

Now we look at the significance of the superiority using one feature over the other in terms of performance. For each group of phenotypic response `yname`, we compute a matrix where each entry indicates the p-value of a one-sided `t.test` testing if using the feature matrix in the row is indeed superior to using the feature matrix in the col.


```r
# significance threshold
thres <- 0.05
# focus only auroc
sname <- "acc"
cat('\n---> \t p-value of t.test showing row superior to col \t <---\n')
# show each grp separately
for (yname in ylist) {
  cat('\n------> \t predicting for ', yname, ' \t <------\n')
  # gather perf
  d <- subset(scores, scores$score == sname & scores$y == yname)
  d <- lapply(split(d, d$x), "[[", "value")
  # order xlist
  xlist.ordered <- as.character(sort(factor(xlist, levels = names(xlist.type), ordered = TRUE)))
  pmatrix <- matrix(-100, nrow = length(xlist), ncol = length(xlist), 
                    dimnames = list(xlist.ordered, xlist.ordered))
  # test now
  for (i in xlist.ordered) {
    for (j in xlist.ordered) {
      tt <- t.test(x = d[[i]], y = d[[j]], alternative = 'greater', mu = 0, paired = TRUE)
      pmatrix[i,j] <- round(tt$p.value, 4)
    }
  }
  # correct p-value for multiple testing with Benjamini-Hochberg
  pmatrix <- p.adjust(pmatrix, "BH")
  print(pmatrix)
  cat('Simplify by thresholding at p-value <', thres, '\n')
  print(pmatrix < thres)
}
```

```
## 
## ---> 	 p-value of t.test showing row superior to col 	 <---
## 
## ------> 	 predicting for  subtype.grps  	 <------
##  [1]          NaN 0.2801400000 0.0000000000 0.0000000000 0.0000000000
##  [6] 0.0000000000 0.0000000000 1.0000000000          NaN 0.0000000000
## [11] 0.0000000000 0.0000000000 0.0000000000 0.0000000000 1.0000000000
## [16] 1.0000000000          NaN 0.0000000000 0.0000000000 0.0000000000
## [21] 0.0000000000 1.0000000000 1.0000000000 1.0000000000          NaN
## [26] 0.0005250000 0.8524000000 0.0112000000 1.0000000000 1.0000000000
## [31] 1.0000000000 1.0000000000          NaN 1.0000000000 1.0000000000
## [36] 1.0000000000 1.0000000000 1.0000000000 1.0000000000 0.0000000000
## [41]          NaN 0.0007411765 1.0000000000 1.0000000000 1.0000000000
## [46] 1.0000000000 0.2801400000 1.0000000000          NaN
## Simplify by thresholding at p-value < 0.05 
##  [1]    NA FALSE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE    NA  TRUE  TRUE
## [12]  TRUE  TRUE  TRUE FALSE FALSE    NA  TRUE  TRUE  TRUE  TRUE FALSE
## [23] FALSE FALSE    NA  TRUE FALSE  TRUE FALSE FALSE FALSE FALSE    NA
## [34] FALSE FALSE FALSE FALSE FALSE FALSE  TRUE    NA  TRUE FALSE FALSE
## [45] FALSE FALSE FALSE FALSE    NA
```

## Session info


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
## other attached packages:
## [1] ggplot2_2.1.0
## 
## loaded via a namespace (and not attached):
##  [1] labeling_0.3     colorspace_1.2-6 scales_0.3.0     plyr_1.8.3      
##  [5] magrittr_1.5     formatR_1.3      tools_3.2.1      gtable_0.1.2    
##  [9] Rcpp_0.12.1      stringi_1.0-1    grid_3.2.1       knitr_1.12.3    
## [13] digest_0.6.8     stringr_1.0.0    munsell_0.4.2    evaluate_0.8.3
```
