# Results
# Yunlong Jiao, 15 Apr 2016

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
## [1] "basal.grps" "surv.grps"  "tumor.grps"
```

```r
# predictors
prlist <- unique(param$prname)
prlist
```

```
## [1] "predictorGBM"        "predictorKNN"        "predictorLDA"       
## [4] "predictorLinearSVM"  "predictorLogitLasso" "predictorNB"        
## [7] "predictorRadialSVM"  "predictorRF"         "predictorSparseSVM"
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
##            y        x      type          predictor rep score     value
## 1 basal.grps eff.vals path-wise predictorLinearSVM   1   acc 0.9696970
## 2 basal.grps eff.vals path-wise predictorLinearSVM   1   fpr 0.0000000
## 3 basal.grps eff.vals path-wise       predictorLDA   1   tpr 0.8000000
## 4 basal.grps eff.vals path-wise predictorLinearSVM   1   ppv 1.0000000
## 5 basal.grps eff.vals path-wise predictorLinearSVM   1  fval 0.9189189
## 6 basal.grps eff.vals path-wise       predictorGBM   1 auroc 0.9917722
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
    facet_wrap(~score, scales = "free_y") + 
    ggtitle(paste0("predicting ", yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
}
```

![plot of chunk overview](result_figure/overview-1.pdf)

```
## Warning: Removed 15 rows containing non-finite values (stat_boxplot).
```

![plot of chunk overview](result_figure/overview-2.pdf)![plot of chunk overview](result_figure/overview-3.pdf)

## Predictor table w.r.t. AUROC

As we see that class sizes are unbalanced so that we focus on AUROC to evaluate performance. This way we also get rid of the effect from tuning cutoff parameter. 

Now we turn to look at which `predictor` suits well in making prediction for a specific group of phenotypic response `yname` using each feature `xname`. We report the frequency of such predictor is selected across 50 repeats of evaluation CV runs.


```r
# focus only auroc
sname <- "auroc"
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
## ---> 	 predicting for  basal.grps  using 	 <---
## $eff.vals
## 
##        predictorGBM  predictorRadialSVM predictorLogitLasso 
##             "62.00"             "18.00"              "8.00" 
##  predictorSparseSVM         predictorRF  predictorLinearSVM 
##              "6.00"              "4.00"              "2.00" 
## 
## $fun.vals
## 
## predictorLogitLasso  predictorRadialSVM         predictorRF 
##             "38.00"             "32.00"             "20.00" 
##  predictorLinearSVM        predictorLDA 
##              "6.00"              "4.00" 
## 
## $genes.vals
## 
## predictorLogitLasso  predictorLinearSVM         predictorRF 
##             "40.00"             "22.00"             "14.00" 
##  predictorRadialSVM        predictorGBM  predictorSparseSVM 
##             "10.00"              "8.00"              "6.00" 
## 
## $go.vals
## 
## predictorLogitLasso  predictorRadialSVM  predictorLinearSVM 
##             "50.00"             "28.00"             "18.00" 
##         predictorRF 
##              "4.00" 
## 
## $mini.genes.vals
## 
##         predictorRF        predictorGBM predictorLogitLasso 
##             "40.00"             "24.00"             "22.00" 
##  predictorRadialSVM 
##             "14.00" 
## 
## $other.genes.vals
## 
## predictorLogitLasso  predictorLinearSVM         predictorRF 
##             "54.00"             "18.00"             "10.00" 
##  predictorSparseSVM        predictorGBM  predictorRadialSVM 
##              "8.00"              "6.00"              "4.00" 
## 
## $path.vals
## 
## predictorLogitLasso        predictorGBM  predictorLinearSVM 
##             "44.00"             "40.00"              "4.00" 
##         predictorRF  predictorSparseSVM        predictorLDA 
##              "4.00"              "4.00"              "2.00" 
##  predictorRadialSVM 
##              "2.00" 
## 
## 
## ---> 	 predicting for  surv.grps  using 	 <---
## $eff.vals
## 
## predictorRadialSVM predictorLinearSVM predictorSparseSVM 
##            "86.00"            "10.00"             "4.00" 
## 
## $fun.vals
## 
## predictorRadialSVM       predictorGBM       predictorLDA 
##            "74.00"            "20.00"             "2.00" 
##        predictorRF predictorSparseSVM 
##             "2.00"             "2.00" 
## 
## $genes.vals
## 
## predictorLinearSVM predictorRadialSVM       predictorLDA 
##            "62.00"            "32.00"             "4.00" 
##        predictorRF 
##             "2.00" 
## 
## $go.vals
## 
## predictorRadialSVM predictorLinearSVM predictorSparseSVM 
##            "52.00"            "38.00"             "8.00" 
##       predictorLDA 
##             "2.00" 
## 
## $mini.genes.vals
## 
## predictorRadialSVM predictorLinearSVM       predictorGBM 
##            "52.00"            "44.00"             "2.00" 
##        predictorRF 
##             "2.00" 
## 
## $other.genes.vals
## 
## predictorLinearSVM predictorRadialSVM       predictorLDA 
##            "60.00"            "34.00"             "4.00" 
##        predictorRF 
##             "2.00" 
## 
## $path.vals
## 
## predictorRadialSVM predictorLinearSVM predictorSparseSVM 
##            "62.00"            "32.00"             "4.00" 
##       predictorLDA 
##             "2.00" 
## 
## 
## ---> 	 predicting for  tumor.grps  using 	 <---
## $eff.vals
## 
##  predictorSparseSVM predictorLogitLasso  predictorLinearSVM 
##             "40.00"             "26.00"             "18.00" 
##         predictorRF  predictorRadialSVM 
##             "10.00"              "6.00" 
## 
## $fun.vals
## 
##         predictorRF predictorLogitLasso        predictorLDA 
##             "38.00"             "14.00"             "12.00" 
##  predictorRadialSVM  predictorSparseSVM  predictorLinearSVM 
##             "12.00"             "12.00"             "10.00" 
##        predictorKNN 
##              "2.00" 
## 
## $genes.vals
## 
##         predictorRF  predictorSparseSVM predictorLogitLasso 
##             "38.00"             "38.00"             "24.00" 
## 
## $go.vals
## 
## predictorLogitLasso  predictorSparseSVM  predictorLinearSVM 
##             "58.00"             "18.00"              "8.00" 
##  predictorRadialSVM        predictorLDA         predictorRF 
##              "6.00"              "4.00"              "4.00" 
##        predictorGBM 
##              "2.00" 
## 
## $mini.genes.vals
## 
##       predictorGBM predictorSparseSVM        predictorRF 
##            "68.00"            "18.00"            "14.00" 
## 
## $other.genes.vals
## 
##  predictorSparseSVM         predictorRF predictorLogitLasso 
##             "48.00"             "32.00"             "20.00" 
## 
## $path.vals
## 
##  predictorSparseSVM  predictorLinearSVM predictorLogitLasso 
##             "46.00"             "16.00"             "14.00" 
##        predictorGBM  predictorRadialSVM         predictorRF 
##             "10.00"              "8.00"              "4.00" 
##        predictorLDA 
##              "2.00"
```

## Significance test w.r.t. AUROC

Now we look at the significance of the superiority using one feature over the other in terms of performance. For each group of phenotypic response `yname`, we compute a matrix where each entry indicates the p-value of a one-sided `t.test` testing if using the feature matrix in the row is indeed superior to using the feature matrix in the col.


```r
# significance threshold
thres <- 0.05
# focus only auroc
sname <- "auroc"
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
      tt <- t.test(x = d[[i]], y = d[[j]], alternative = 'greater', mu = 0)
      pmatrix[i,j] <- round(tt$p.value, 4)
    }
  }
  print(pmatrix)
  cat('Simplify by thresholding at p-value <', thres, '\n')
  print(pmatrix < thres)
}
```

```
## 
## ---> 	 p-value of t.test showing row superior to col 	 <---
## 
## ------> 	 predicting for  basal.grps  	 <------
##                  fun.vals go.vals eff.vals path.vals mini.genes.vals
## fun.vals           0.5000  0.8583   0.9381    0.8788          0.8377
## go.vals            0.1417  0.5000   0.7522    0.6198          0.5210
## eff.vals           0.0619  0.2478   0.5000    0.3783          0.2876
## path.vals          0.1212  0.3802   0.6217    0.5000          0.4095
## mini.genes.vals    0.1623  0.4790   0.7124    0.5905          0.5000
## other.genes.vals   0.0714  0.2715   0.5233    0.4009          0.3097
## genes.vals         0.1031  0.3367   0.5787    0.4583          0.3688
##                  other.genes.vals genes.vals
## fun.vals                   0.9286     0.8969
## go.vals                    0.7285     0.6633
## eff.vals                   0.4767     0.4213
## path.vals                  0.5991     0.5417
## mini.genes.vals            0.6903     0.6312
## other.genes.vals           0.5000     0.4439
## genes.vals                 0.5561     0.5000
## Simplify by thresholding at p-value < 0.05 
##                  fun.vals go.vals eff.vals path.vals mini.genes.vals
## fun.vals            FALSE   FALSE    FALSE     FALSE           FALSE
## go.vals             FALSE   FALSE    FALSE     FALSE           FALSE
## eff.vals            FALSE   FALSE    FALSE     FALSE           FALSE
## path.vals           FALSE   FALSE    FALSE     FALSE           FALSE
## mini.genes.vals     FALSE   FALSE    FALSE     FALSE           FALSE
## other.genes.vals    FALSE   FALSE    FALSE     FALSE           FALSE
## genes.vals          FALSE   FALSE    FALSE     FALSE           FALSE
##                  other.genes.vals genes.vals
## fun.vals                    FALSE      FALSE
## go.vals                     FALSE      FALSE
## eff.vals                    FALSE      FALSE
## path.vals                   FALSE      FALSE
## mini.genes.vals             FALSE      FALSE
## other.genes.vals            FALSE      FALSE
## genes.vals                  FALSE      FALSE
## 
## ------> 	 predicting for  surv.grps  	 <------
##                  fun.vals go.vals eff.vals path.vals mini.genes.vals
## fun.vals           0.5000  0.1194   0.8132    0.9653          0.3963
## go.vals            0.8806  0.5000   0.9809    0.9977          0.8306
## eff.vals           0.1868  0.0191   0.5000    0.8595          0.1166
## path.vals          0.0347  0.0023   0.1405    0.5000          0.0179
## mini.genes.vals    0.6037  0.1694   0.8834    0.9821          0.5000
## other.genes.vals   0.1870  0.0197   0.4956    0.8546          0.1176
## genes.vals         0.1660  0.0168   0.4569    0.8313          0.1027
##                  other.genes.vals genes.vals
## fun.vals                   0.8130     0.8340
## go.vals                    0.9803     0.9832
## eff.vals                   0.5044     0.5431
## path.vals                  0.1454     0.1687
## mini.genes.vals            0.8824     0.8973
## other.genes.vals           0.5000     0.5383
## genes.vals                 0.4617     0.5000
## Simplify by thresholding at p-value < 0.05 
##                  fun.vals go.vals eff.vals path.vals mini.genes.vals
## fun.vals            FALSE   FALSE    FALSE     FALSE           FALSE
## go.vals             FALSE   FALSE    FALSE     FALSE           FALSE
## eff.vals            FALSE    TRUE    FALSE     FALSE           FALSE
## path.vals            TRUE    TRUE    FALSE     FALSE            TRUE
## mini.genes.vals     FALSE   FALSE    FALSE     FALSE           FALSE
## other.genes.vals    FALSE    TRUE    FALSE     FALSE           FALSE
## genes.vals          FALSE    TRUE    FALSE     FALSE           FALSE
##                  other.genes.vals genes.vals
## fun.vals                    FALSE      FALSE
## go.vals                     FALSE      FALSE
## eff.vals                    FALSE      FALSE
## path.vals                   FALSE      FALSE
## mini.genes.vals             FALSE      FALSE
## other.genes.vals            FALSE      FALSE
## genes.vals                  FALSE      FALSE
## 
## ------> 	 predicting for  tumor.grps  	 <------
##                  fun.vals go.vals eff.vals path.vals mini.genes.vals
## fun.vals           0.5000  0.9239   0.9990    0.9611          0.9999
## go.vals            0.0761  0.5000   0.8495    0.5275          0.9653
## eff.vals           0.0010  0.1505   0.5000    0.0914          0.9929
## path.vals          0.0389  0.4725   0.9086    0.5000          0.9920
## mini.genes.vals    0.0001  0.0347   0.0071    0.0080          0.5000
## other.genes.vals   0.0000  0.0135   0.0000    0.0013          0.0557
## genes.vals         0.0000  0.0188   0.0000    0.0025          0.1709
##                  other.genes.vals genes.vals
## fun.vals                   1.0000     1.0000
## go.vals                    0.9865     0.9812
## eff.vals                   1.0000     1.0000
## path.vals                  0.9987     0.9975
## mini.genes.vals            0.9443     0.8291
## other.genes.vals           0.5000     0.1152
## genes.vals                 0.8848     0.5000
## Simplify by thresholding at p-value < 0.05 
##                  fun.vals go.vals eff.vals path.vals mini.genes.vals
## fun.vals            FALSE   FALSE    FALSE     FALSE           FALSE
## go.vals             FALSE   FALSE    FALSE     FALSE           FALSE
## eff.vals             TRUE   FALSE    FALSE     FALSE           FALSE
## path.vals            TRUE   FALSE    FALSE     FALSE           FALSE
## mini.genes.vals      TRUE    TRUE     TRUE      TRUE           FALSE
## other.genes.vals     TRUE    TRUE     TRUE      TRUE           FALSE
## genes.vals           TRUE    TRUE     TRUE      TRUE           FALSE
##                  other.genes.vals genes.vals
## fun.vals                    FALSE      FALSE
## go.vals                     FALSE      FALSE
## eff.vals                    FALSE      FALSE
## path.vals                   FALSE      FALSE
## mini.genes.vals             FALSE      FALSE
## other.genes.vals            FALSE      FALSE
## genes.vals                  FALSE      FALSE
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
