# Results
# Yunlong Jiao, 15 Apr 2016

This script collects results from earlier run and illustrates with tables and plots.


```r
knitr::opts_chunk$set(error = FALSE, fig.width = 12, fig.height = 8, dev = c("png","pdf"), fig.keep = "high", fig.path = "3_results_figure/", cache.path = "3_results_cache/")
set.seed(70236562)
source("../../src/func.R")
library(ggplot2)
library(reshape2)
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
##  [1] "eff.vals"                  "fun.vals"                 
##  [3] "genes.vals"                "go.vals"                  
##  [5] "mini.genes.vals"           "other.genes.vals"         
##  [7] "path.vals"                 "eff.and.other.genes.vals" 
##  [9] "path.and.other.genes.vals" "path.and.genes.vals"
```

```r
# feature types
xlist.type <- c("func-wise", "func-wise", 
                "path-wise", "path-wise", 
                "gene-wise", "gene-wise", "gene-wise", 
                "mix-wise", "mix-wise", "mix-wise")
names(xlist.type) <- c("fun.vals", "go.vals", # functionality features
                       "eff.vals", "path.vals", # pathway features
                       "mini.genes.vals", "other.genes.vals", "genes.vals", # gene features
                       "eff.and.other.genes.vals", "path.and.other.genes.vals", "path.and.genes.vals") # mixed type
xlist.vline <- c(2.5, 4.5, 7.5) # cut out types
stopifnot(length(setdiff(xlist, names(xlist.type))) == 0)
xlist.type
```

```
##                  fun.vals                   go.vals 
##               "func-wise"               "func-wise" 
##                  eff.vals                 path.vals 
##               "path-wise"               "path-wise" 
##           mini.genes.vals          other.genes.vals 
##               "gene-wise"               "gene-wise" 
##                genes.vals  eff.and.other.genes.vals 
##               "gene-wise"                "mix-wise" 
## path.and.other.genes.vals       path.and.genes.vals 
##                "mix-wise"                "mix-wise"
```

```r
# groups
ylist <- unique(param$yname)
ylist
```

```
## [1] "surv.grps"
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
slist <- c("acc","fpr","tpr","ppv","fval","auroc")
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
      # get obj
      objname <- paste('res', xname, yname, i.fold, nfolds, nrepeats, sep = '_')
      objpath <- paste0('Robj/',objname,'.RData')
      if (file.exists(objpath)) {
        ivres <- get(load(objpath))
      } else {
        # read nested cv res
        cvres <- list()
        for (prname in prlist) {
          res.files <- list.files(path = 'Robj', 
                                  pattern = paste('^cvres', xname, yname, prname, 
                                                  i.fold, nfolds, nrepeats, 
                                                  '[[:digit:]]+', nfolds.inn, nrepeats.inn, 
                                                  sep = '_'), 
                                  full.names = TRUE)
          res <- lapply(res.files, function(f) try(get(load(f))))
          res <- res[sapply(res, function(x) !inherits(x, "try-error"))]
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
        
        # save up
        assign(objname, ivres)
        save(list = objname, file = objpath)
        rm(list = c(objname,"cvres"))
        gc()
      }
      
      # record score values
      scores[[objname]] <- data.frame(y = yname, 
                                      x = xname, 
                                      type = xlist.type[xname], 
                                      predictor = sapply(slist, function(sname) ivres[[sname]][["predictor"]]), 
                                      rep = i.fold, 
                                      score = slist, 
                                      value = sapply(slist, function(sname) ivres[[sname]][[sname]]), 
                                      row.names = NULL)
      rm(ivres)
      gc()
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
##           y        x      type           predictor rep score     value
## 1 surv.grps eff.vals path-wise  predictorRadialSVM   1   acc 0.8750000
## 2 surv.grps eff.vals path-wise   predictorConstant   1   fpr 0.0000000
## 3 surv.grps eff.vals path-wise        predictorLDA   1   tpr 0.8095238
## 4 surv.grps eff.vals path-wise predictorKendallSVM   1   ppv       NaN
## 5 surv.grps eff.vals path-wise  predictorSparseSVM   1  fval 0.2641509
## 6 surv.grps eff.vals path-wise  predictorLinearSVM   1 auroc 0.7462366
```

## Performance overview

In this section, each plot corresponds to the prediction performance for one specific group of phenotypic response. In each plot, along x-axis we have different feature matrices, along y-axis we have different evaluation measures, the barplot shows the CV scores evaluated against randomly splitted CV folds.

Notes that some evaluation measures (`acc`, `ppv`, etc) do associate with a specific `cutoff` threshold for predicted probability, which is always set to constant 0.5. It might be interesting to tune the threshold in case that the class sizes in training data are unbalanced.

The evaluation measure are the following: acc, fpr, tpr, ppv, fval, auroc. Note that in order to remove confusion, we equivalently rename `TPR` by `Sensitivity`, and plot `Specificity = 1 - FPR` instead of `FPR` values. As a result, for any evaluation measure present in the plot, higher values indicate better performance unanimously.


```r
scores <- read.table("scores.txt", header = TRUE, row.names = 1, col.names = c("y", "x", "type", "predictor", "rep", "score", "value"), sep = '\t')
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
    facet_wrap(~score, scales = "free") + 
    ggtitle(paste0("predicting ", yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
}
```

```
## Warning: Removed 178 rows containing non-finite values (stat_boxplot).
```

![plot of chunk overview](3_results_figure/overview-1.png)

## Predictor table w.r.t. AUROC

As we see that class sizes are unbalanced so that we focus on AUROC to evaluate performance. This way we also get rid of the effect from tuning cutoff parameter. 

Now we turn to look at which `predictor` suits well in making prediction for a specific group of phenotypic response `yname` using each feature `xname`. We report the frequency of such predictor is selected across 50 repeats of evaluation CV runs.


```r
# focus only auroc
sname <- "auroc"
# show each grp separately
for (yname in ylist) {
  cat('\n---> \t Frequency (%) of predictor being selected \t <---\n')
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
## ---> 	 predicting for  surv.grps  using 	 <---
## $eff.and.other.genes.vals
## 
## predictorLinearSVM predictorRadialSVM       predictorLDA 
##            "68.00"            "24.00"             "4.00" 
##       predictorGBM        predictorRF 
##             "2.00"             "2.00" 
## 
## $eff.vals
## 
## predictorRadialSVM predictorLinearSVM predictorSparseSVM 
##            "86.00"            "10.00"             "4.00" 
## 
## $fun.vals
## 
##  predictorRadialSVM        predictorGBM predictorKendallSVM 
##             "72.00"             "20.00"              "4.00" 
##        predictorLDA  predictorSparseSVM 
##              "2.00"              "2.00" 
## 
## $genes.vals
## 
## predictorLinearSVM predictorRadialSVM       predictorLDA 
##            "60.00"            "32.00"             "4.00" 
##        predictorRF 
##             "4.00" 
## 
## $go.vals
## 
##  predictorRadialSVM  predictorLinearSVM  predictorSparseSVM 
##             "50.00"             "38.00"              "8.00" 
## predictorKendallSVM        predictorLDA 
##              "2.00"              "2.00" 
## 
## $mini.genes.vals
## 
##  predictorRadialSVM  predictorLinearSVM        predictorGBM 
##             "52.00"             "42.00"              "2.00" 
## predictorKendallSVM         predictorRF 
##              "2.00"              "2.00" 
## 
## $other.genes.vals
## 
##  predictorLinearSVM  predictorRadialSVM        predictorLDA 
##             "60.00"             "30.00"              "4.00" 
##         predictorRF predictorKendallSVM 
##              "4.00"              "2.00" 
## 
## $path.and.genes.vals
## 
## predictorLinearSVM predictorRadialSVM       predictorLDA 
##            "50.00"            "34.00"            "14.00" 
##       predictorGBM 
##             "2.00" 
## 
## $path.and.other.genes.vals
## 
## predictorLinearSVM predictorRadialSVM 
##            "64.00"            "36.00" 
## 
## $path.vals
## 
## predictorRadialSVM predictorLinearSVM predictorSparseSVM 
##            "62.00"            "32.00"             "4.00" 
##       predictorLDA 
##             "2.00"
```

## Significance test w.r.t. AUROC

Now we look at the significance of the superiority using one feature over the other in terms of performance. For each group of phenotypic response `yname`, we compute a matrix where each entry indicates the p-value of a one-sided `t.test` testing if using the feature matrix in the row is indeed superior to using the feature matrix in the col.


```r
# significance threshold
thres <- 0.05
# focus only auroc
sname <- "auroc"
# show each grp separately
for (yname in ylist) {
  cat('\n---> \t p-value of t.test showing row superior to col \t <---\n')
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
  pmatrix.adj <- p.adjust(pmatrix, "BH")
  attributes(pmatrix.adj) <- attributes(pmatrix)
  print(pmatrix.adj)
  cat('\nSimplify by thresholding at p-value <', thres, '\n')
  print(pmatrix.adj < thres)
}
```

```
## 
## ---> 	 p-value of t.test showing row superior to col 	 <---
## 
## ------> 	 predicting for  surv.grps  	 <------
##                             fun.vals     go.vals   eff.vals path.vals
## fun.vals                         NaN 0.157178571 1.00000000   1.00000
## go.vals                   1.00000000         NaN 1.00000000   1.00000
## eff.vals                  0.13716000 0.003000000        NaN   1.00000
## path.vals                 0.01200000 0.000000000 0.08475000       NaN
## mini.genes.vals           1.00000000 0.167896552 1.00000000   1.00000
## other.genes.vals          0.20758065 0.012375000 1.00000000   1.00000
## genes.vals                0.20758065 0.012000000 1.00000000   1.00000
## eff.and.other.genes.vals  0.08475000 0.004909091 0.64215000   1.00000
## path.and.other.genes.vals 0.00180000 0.000000000 0.02435294   0.64215
## path.and.genes.vals       0.06586364 0.001800000 0.54781579   1.00000
##                           mini.genes.vals other.genes.vals  genes.vals
## fun.vals                      0.831731707        1.0000000 1.000000000
## go.vals                       1.000000000        1.0000000 1.000000000
## eff.vals                      0.027000000        0.8895349 0.889534884
## path.vals                     0.004909091        0.1500000 0.147461538
## mini.genes.vals                       NaN        1.0000000 1.000000000
## other.genes.vals              0.027000000              NaN 0.948200000
## genes.vals                    0.027000000        1.0000000         NaN
## eff.and.other.genes.vals      0.004909091        0.2380909 0.218812500
## path.and.other.genes.vals     0.000000000        0.0045000 0.003857143
## path.and.genes.vals           0.008250000        0.2928857 0.288264706
##                           eff.and.other.genes.vals
## fun.vals                                1.00000000
## go.vals                                 1.00000000
## eff.vals                                1.00000000
## path.vals                               0.37508108
## mini.genes.vals                         1.00000000
## other.genes.vals                        1.00000000
## genes.vals                              1.00000000
## eff.and.other.genes.vals                       NaN
## path.and.other.genes.vals               0.05057143
## path.and.genes.vals                     0.90879545
##                           path.and.other.genes.vals path.and.genes.vals
## fun.vals                                          1               1.000
## go.vals                                           1               1.000
## eff.vals                                          1               1.000
## path.vals                                         1               0.352
## mini.genes.vals                                   1               1.000
## other.genes.vals                                  1               1.000
## genes.vals                                        1               1.000
## eff.and.other.genes.vals                          1               1.000
## path.and.other.genes.vals                       NaN               0.012
## path.and.genes.vals                               1                 NaN
## 
## Simplify by thresholding at p-value < 0.05 
##                           fun.vals go.vals eff.vals path.vals
## fun.vals                        NA   FALSE    FALSE     FALSE
## go.vals                      FALSE      NA    FALSE     FALSE
## eff.vals                     FALSE    TRUE       NA     FALSE
## path.vals                     TRUE    TRUE    FALSE        NA
## mini.genes.vals              FALSE   FALSE    FALSE     FALSE
## other.genes.vals             FALSE    TRUE    FALSE     FALSE
## genes.vals                   FALSE    TRUE    FALSE     FALSE
## eff.and.other.genes.vals     FALSE    TRUE    FALSE     FALSE
## path.and.other.genes.vals     TRUE    TRUE     TRUE     FALSE
## path.and.genes.vals          FALSE    TRUE    FALSE     FALSE
##                           mini.genes.vals other.genes.vals genes.vals
## fun.vals                            FALSE            FALSE      FALSE
## go.vals                             FALSE            FALSE      FALSE
## eff.vals                             TRUE            FALSE      FALSE
## path.vals                            TRUE            FALSE      FALSE
## mini.genes.vals                        NA            FALSE      FALSE
## other.genes.vals                     TRUE               NA      FALSE
## genes.vals                           TRUE            FALSE         NA
## eff.and.other.genes.vals             TRUE            FALSE      FALSE
## path.and.other.genes.vals            TRUE             TRUE       TRUE
## path.and.genes.vals                  TRUE            FALSE      FALSE
##                           eff.and.other.genes.vals
## fun.vals                                     FALSE
## go.vals                                      FALSE
## eff.vals                                     FALSE
## path.vals                                    FALSE
## mini.genes.vals                              FALSE
## other.genes.vals                             FALSE
## genes.vals                                   FALSE
## eff.and.other.genes.vals                        NA
## path.and.other.genes.vals                    FALSE
## path.and.genes.vals                          FALSE
##                           path.and.other.genes.vals path.and.genes.vals
## fun.vals                                      FALSE               FALSE
## go.vals                                       FALSE               FALSE
## eff.vals                                      FALSE               FALSE
## path.vals                                     FALSE               FALSE
## mini.genes.vals                               FALSE               FALSE
## other.genes.vals                              FALSE               FALSE
## genes.vals                                    FALSE               FALSE
## eff.and.other.genes.vals                      FALSE               FALSE
## path.and.other.genes.vals                        NA                TRUE
## path.and.genes.vals                           FALSE                  NA
```

## Complementary section: performance with random forest only

This section is complementary to previous section "Performance overview". The difference is that here we do not tune predictor out but fix the predictor to be random forest. Then we compare performance of different feature types wrt various criteria. Generally RF should do a fairly good job in classification while enabling simultaneously feature selection.


```r
prname <- "predictorRF"
scores <- list()
for (yname in ylist) {
  for (xname in xlist) {
    message(yname, '\t', xname)
    res.files <- list.files(path = 'Robj', 
                            pattern = paste('^ivres', xname, yname, prname, 
                                            '[[:digit:]]+', nfolds, nrepeats, 
                                            0, nfolds.inn, nrepeats.inn, 
                                            sep = '_'), 
                            full.names = TRUE)
    res <- lapply(res.files, function(f) get(load(f)))
    names(res) <- paste0("rep",1:length(res))
    # reduce to only scores and melt
    res <- lapply(res, '[', slist)
    res <- melt(res)
    colnames(res)[colnames(res) == "L1"] <- "rep"
    colnames(res)[colnames(res) == "L2"] <- "score"
    
    # record score values
    objname <- paste('res', xname, yname, prname, nfolds, nrepeats, sep = '_')
    scores[[objname]] <- data.frame(y = yname, 
                                    x = xname, 
                                    type = xlist.type[xname], 
                                    predictor = prname, 
                                    res, 
                                    row.names = NULL)
  }
}
scores <- do.call('rbind', scores)
rownames(scores) <- seq(nrow(scores))
head(scores)
```

```
##           y        x      type   predictor      value score  rep
## 1 surv.grps eff.vals path-wise predictorRF 0.89204545   acc rep1
## 2 surv.grps eff.vals path-wise predictorRF 0.01290323   fpr rep1
## 3 surv.grps eff.vals path-wise predictorRF 0.19047619   tpr rep1
## 4 surv.grps eff.vals path-wise predictorRF 0.66666667   ppv rep1
## 5 surv.grps eff.vals path-wise predictorRF 0.29629630  fval rep1
## 6 surv.grps eff.vals path-wise predictorRF 0.66451613 auroc rep1
```


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
    facet_wrap(~score, scales = "free") + 
    ggtitle(paste0("predicting ", yname, " using ", prname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
}
```

```
## Warning: Removed 32 rows containing non-finite values (stat_boxplot).
```

![plot of chunk perfRF](3_results_figure/perfRF-1.png)

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
## [1] reshape2_1.4.1 ggplot2_2.1.0 
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.1      digest_0.6.8     grid_3.2.1       plyr_1.8.3      
##  [5] gtable_0.1.2     formatR_1.3      magrittr_1.5     evaluate_0.8.3  
##  [9] scales_0.3.0     stringi_1.0-1    labeling_0.3     tools_3.2.1     
## [13] stringr_1.0.0    munsell_0.4.2    colorspace_1.2-6 knitr_1.12.3
```
