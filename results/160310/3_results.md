# Results
# Yunlong Jiao, 10 March 2016

This script collects results from earlier run and illustrates with tables and plots.


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, message = FALSE, fig.width = 12, fig.height = 8, dev = c("png","pdf"), fig.keep = "high", fig.path = "result_figure/", cache.path = "result_cache/")
set.seed(70236562)
source("../../src/func.R")
library(ggplot2)
```

# Read in scores


```r
param <- read.table("cluster_param.txt", header = FALSE, row.names = NULL, col.names = c("idx", "xname", "yname", "prname", "rep", "nfolds", "nrepeats"))
xlist <- unique(param$xname)
ylist <- unique(param$yname)
prlist <- unique(param$prname)
nfolds <- unique(param$nfolds); stopifnot(length(nfolds) == 1)
nrepeats <- unique(param$nrepeats); stopifnot(length(nrepeats) == 1)
slist <- c("acc","fpr","tpr","ppv","fval","concordance.index","auroc")

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
      assign(objname, crossValidationCombineResults(cvres))
      # plot ROC
      if (!dir.exists('figures')) dir.create('figures')
      plotROCcv(res = get(objname), savepath = paste0('figures/', objname, '.pdf'))
      # write in scores
      ss <- get(objname)[slist]
      ss[sapply(ss, is.null)] <- NA
      scores[[objname]] <- data.frame(x = xname, 
                                      y = yname, 
                                      predictor = prname, 
                                      score = slist, 
                                      value = unlist(ss), 
                                      n.cv.folds = nfiles,
                                      row.names = NULL)
    }
  }
}
scores <- do.call('rbind', scores)
rownames(scores) <- seq(nrow(scores))
# write out
write.table(scores, file = "scores.txt", row.names = TRUE, col.names = TRUE, sep = '\t')
# a bit more pruning for plotting
scores <- subset(scores, scores$score != "concordance.index")
scores$value <- round(scores$value, 3)
head(scores)
```

```
##          x          y    predictor score value n.cv.folds
## 1 eff.vals basal.grps predictorGBM   acc 0.941         50
## 2 eff.vals basal.grps predictorGBM   fpr 0.010         50
## 3 eff.vals basal.grps predictorGBM   tpr 0.733         50
## 4 eff.vals basal.grps predictorGBM   ppv 0.946         50
## 5 eff.vals basal.grps predictorGBM  fval 0.822         50
## 7 eff.vals basal.grps predictorGBM auroc 0.968         50
```

# Overview

We show the score values for different each feature types (with variance of boxplots across different methods in boxplot).


```r
# plot each grps in a separate figure
for (yname in ylist) {
  d <- subset(scores, scores$y == yname)
  p1 <- ggplot(d, aes(x = x, y = value)) + 
    geom_boxplot(aes(fill = x), alpha = 0.8) + 
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


```r
# focus only one type of scores
key <- "auroc"
# plot each grps in a separate figure
for (yname in ylist) {
  d <- subset(scores, scores$score == key & scores$y == yname)
  yrange <- c(min(d$value)*0.95, 1)
  p1 <- ggplot(d, aes(x = x, y = value)) + 
    geom_bar(aes(fill = x), stat = "identity", position = "dodge") + 
    geom_text(aes(label = value), vjust = -0.3, colour = "black", 
              position = position_dodge(0.9), size = 4) + 
    facet_wrap(~predictor) + 
    coord_cartesian(ylim = yrange) + 
    ggtitle(paste0(key, " for predicting ", yname)) + 
    theme(axis.text.x = element_blank())
  plot(p1)
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
## [1] ROCR_1.0-7    gplots_2.17.0 ggplot2_1.0.1 knitr_1.12.3 
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.1        gtools_3.5.0       digest_0.6.8      
##  [4] bitops_1.0-6       MASS_7.3-44        grid_3.2.1        
##  [7] plyr_1.8.3         gtable_0.1.2       formatR_1.3       
## [10] magrittr_1.5       evaluate_0.8.3     scales_0.3.0      
## [13] KernSmooth_2.23-15 stringi_1.0-1      reshape2_1.4.1    
## [16] gdata_2.17.0       labeling_0.3       proto_0.3-10      
## [19] tools_3.2.1        stringr_1.0.0      munsell_0.4.2     
## [22] colorspace_1.2-6   caTools_1.17.1
```
