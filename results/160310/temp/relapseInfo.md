# Process ICGC BRCA data
# Yunlong Jiao, 13 Apr 2016

This script shows how incomplete relapse info from ICGC BRCA dataset is so that survival prediction cannot proceed on this dataset...


```r
knitr::opts_chunk$set(error = FALSE)
datapath <- '../../../data/BRCA_saved_data/'
options(stringsAsFactors = FALSE)
```


```r
# samples
load(paste0(datapath, '110_genes_vals.RData'))
samplelist <- colnames(genes.vals)
# align
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
alignid <- subset(alignid, !is.na(alignid$donor))
id <- match(alignid$donor, y$icgc_donor_id)
alignid <- alignid[!is.na(id), ]
id <- id[!is.na(id)]
y <- y[id, ]
```


```r
# number of patients with clinical info available 
# after alignment between info_donor and info_sample
nrow(y)
```

```
## [1] 879
```

```r
# everyone has donor vital status
table(y$donor_vital_status)
```

```
## 
##    alive deceased 
##      755      124
```

```r
# only two patients have concrete survival time
table(y$donor_survival_time)
```

```
## 
## 315 412 
##   1   1
```

```r
# most people have relapse status but data is highly unbalanced
table(y$disease_status_last_followup)
```

```
## 
##                    complete remission        progression 
##                123                680                 76
```

```r
# followup period ranges wildly between patients
# (data is not fair to make binary prediction on binary relapse status?)
quantile(y$donor_interval_of_last_followup, na.rm = T)
```

```
##     0%    25%    50%    75%   100% 
##   -7.0  120.5  412.5 1154.5 6796.0
```

```r
# vital status is overall survival but not disease-specific?
table(y[,c("donor_vital_status","disease_status_last_followup")])
```

```
##                   disease_status_last_followup
## donor_vital_status     complete remission progression
##           alive     93                642          20
##           deceased  30                 38          56
```

```r
# 5-year true relapses is even rare in data
i <- which(y$donor_interval_of_last_followup < 365*5)
length(i) # number of followup less than 5 years
```

```
## [1] 683
```

```r
table(y$disease_status_last_followup[i])
```

```
## 
##                    complete remission        progression 
##                 81                584                 18
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
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## loaded via a namespace (and not attached):
## [1] magrittr_1.5  formatR_1.2.1 tools_3.2.3   stringi_1.0-1 knitr_1.12.3 
## [6] methods_3.2.3 stringr_1.0.0 evaluate_0.8
```
