# Results
# Yunlong Jiao, 15 Apr 2016

This script collects results from earlier run and illustrates with tables and plots.


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, message = FALSE, fig.width = 12, fig.height = 8, dev = c("png","pdf"), fig.keep = "high", fig.path = "result_figure/", cache.path = "result_cache/")
set.seed(70236562)
source("../../src/func.R")
library(ggplot2)
```

Read in parameters.


```r
# read in parameters
param <- read.table("cluster_param.txt", header = FALSE, row.names = NULL, col.names = c("idx", "xname", "yname", "prname", "rep", "nfolds", "nrepeats"))
xlist <- unique(param$xname)
xlist
```

```
## [1] "eff.vals"         "fun.vals"         "genes.vals"      
## [4] "go.vals"          "mini.genes.vals"  "other.genes.vals"
## [7] "path.vals"
```

```r
ylist <- unique(param$yname)
ylist
```

```
## [1] "basal.grps" "surv.grps"  "tumor.grps"
```

```r
prlist <- unique(param$prname)
prlist
```

```
## [1] "predictorGBM"        "predictorKNN"        "predictorLDA"       
## [4] "predictorLinearSVM"  "predictorLogitLasso" "predictorNB"        
## [7] "predictorRadialSVM"  "predictorRF"         "predictorSparseSVM"
```

```r
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
slist <- c("acc","fpr","tpr","ppv","fval","auroc")
slist
```

```
## [1] "acc"   "fpr"   "tpr"   "ppv"   "fval"  "auroc"
```

Read in CV folds results.


```r
# gather scores
scores <- list()
for (xname in xlist) {
  for (yname in ylist){
    for (prname in prlist) {
      message(".", appendLF = FALSE)
      objname <- paste('res', xname, yname, prname, nfolds, nrepeats, sep = '_')
      # read and combine cv results from cluster runs
      cvres.files <- list.files(path = "Robj", pattern = objname, full.names = TRUE)
      nfiles <- length(cvres.files) # number of cv jobs successfully done
      cvres <- lapply(cvres.files, function(f) get(load(f)))
      names(cvres) <- cvres.files
      assign(objname, cvres)
      # plot ROC
      cvres <- crossValidationCombineResults(cvres)
      if (!dir.exists('figures')) dir.create('figures')
      plotROCcv(res = cvres, savepath = paste0('figures/', objname, '.pdf'))
      # write in scores
      ss <- cvres[slist]
      ss[sapply(ss, is.null)] <- NA
      scores[[objname]] <- data.frame(x = xname, 
                                      y = yname, 
                                      predictor = prname, 
                                      score = slist, 
                                      value = unlist(ss), 
                                      n.avg.cvfolds = nfiles,
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
##          x          y    predictor score      value n.avg.cvfolds
## 1 eff.vals basal.grps predictorGBM   acc 0.94103512            50
## 2 eff.vals basal.grps predictorGBM   fpr 0.01022524            50
## 3 eff.vals basal.grps predictorGBM   tpr 0.73268066            50
## 4 eff.vals basal.grps predictorGBM   ppv 0.94597975            50
## 5 eff.vals basal.grps predictorGBM  fval 0.82207255            50
## 6 eff.vals basal.grps predictorGBM auroc 0.96812517            50
```

# Overview

One plot corresponds to predicting for one type of phenotypic responses. In each plot, along x-axis we have different feature matrices, along y-axis we have different evaluation measures, the barplot shows the maximum CV scores from different predictors indicating the best predictor suited for such features.

Notes:

- `cutoff` threshold for predicted probability is always set constant 0.5, although it might be interesting to tune it in case that the class sizes in training data are unbalanced.
- The CV evaluated scores are overfitting typically for models with tuning parameters. A decent solution would be to set `prname` as one of the tuning parameters in nested CV call, or simply use another independent dataset to evaluate the scores.


```r
# plot each grps in a separate figure
for (yname in ylist) {
  d <- subset(scores, scores$y == yname)
  p1 <- ggplot(d, aes(x = x, y = value)) + 
    stat_summary_bin(aes(fill = x), fun.y = "max", geom = "bar") + 
    stat_summary(aes(label = round(..y.., 3)), fun.y = "max", geom = "text", size = 4,
                 colour = "black", position = position_dodge(0.9), vjust = -0.3) + 
    facet_wrap(~score) + 
    coord_cartesian(ylim = c(0, 1)) + 
    ggtitle(paste0("predicting ", yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
}
```

![plot of chunk overview](result_figure/overview-1.png)![plot of chunk overview](result_figure/overview-2.png)![plot of chunk overview](result_figure/overview-3.png)

# AUROC

As we see that class sizes are unbalanced so that we focus on AUROC to evaluate performance.

One plot corresponds to predicting for one type of phenotypic responses. In each plot, along x-axis we have different feature matrices, along y-axis we have CV AUROC scores, the boxplot shows the variance from different CV folds indicating the performance when data is randomly split.


```r
# focus only one type of scores
sname <- "auroc"
# plot each grps in a separate figure
for (yname in ylist) {
  # gather scores
  d <- list()
  for (xname in xlist) {
    for (prname in prlist) {
      objname <- paste('res', xname, yname, prname, nfolds, nrepeats, sep = '_')
      cvres <- get(objname)
      ss <- lapply(cvres, function(u) u[[sname]])
      ss[sapply(ss, is.null)] <- NA
      d[[objname]] <- data.frame(x = xname, 
                                 y = yname, 
                                 predictor = prname, 
                                 score = sname, 
                                 value = unlist(ss), 
                                 rep = 1:length(ss), 
                                 row.names = NULL)
    }
  }
  d <- do.call('rbind', d)
  rownames(d) <- seq(nrow(d))
  # plot
  p1 <- ggplot(d, aes(x = x, y = value)) + 
    geom_boxplot(aes(fill = x), alpha = 0.8) + 
    facet_wrap(~predictor) + 
    ggtitle(paste0(sname, " for predicting ", yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
  ### TODO t-test on boxplot to quantify superiority of one feature type over another?
}
```

![plot of chunk auroc](result_figure/auroc-1.png)![plot of chunk auroc](result_figure/auroc-2.png)![plot of chunk auroc](result_figure/auroc-3.png)

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
## other attached packages:
## [1] ggplot2_2.1.0
## 
## loaded via a namespace (and not attached):
##  [1] labeling_0.3     colorspace_1.2-6 scales_0.3.0     plyr_1.8.3      
##  [5] magrittr_1.5     formatR_1.3      tools_3.2.1      gtable_0.1.2    
##  [9] Rcpp_0.12.1      stringi_1.0-1    grid_3.2.1       knitr_1.12.3    
## [13] digest_0.6.8     stringr_1.0.0    munsell_0.4.2    evaluate_0.8.3
```
