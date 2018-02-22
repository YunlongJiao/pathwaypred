# Process ICGC BRCA data
# Yunlong Jiao, 13 Apr 2016, UPDATED 22 FEB 2018

This script shows how incomplete relapse info from ICGC BRCA dataset is so that survival prediction cannot proceed on this dataset...

NOTE survival info has been completed in later releases of the data by ICGC.


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
alignmat <- read.table(paste0(datapath, 'sample.BRCA-US.tsv'), sep = "\t", quote = "\"", header = TRUE)
y <- read.table(paste0(datapath, 'donor.BRCA-US.tsv'), sep = "\t", quote = "\"", header = TRUE)
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
sum(!is.na(y$donor_interval_of_last_followup))
```

```
## [1] 755
```

```r
quantile(y$donor_survival_time, na.rm = T)
```

```
##      0%     25%     50%     75%    100% 
##  158.00  851.25 1563.00 2384.00 4456.00
```

```r
# followup period ranges wildly between patients
# (data is not fair to make binary prediction on binary relapse status?)
sum(!is.na(y$donor_survival_time))
```

```
## [1] 124
```

```r
quantile(y$donor_interval_of_last_followup, na.rm = T)
```

```
##   0%  25%  50%  75% 100% 
##   -7  120  411 1155 7067
```

```r
# most people have relapse status but data is highly unbalanced
table(y$disease_status_last_followup)
```

```
## 
##                    complete remission        progression 
##                200                658                 21
```

```r
# vital status is overall survival but not disease-specific?
table(y[,c("donor_vital_status","disease_status_last_followup")])
```

```
##                   disease_status_last_followup
## donor_vital_status     complete remission progression
##           alive     76                658          21
##           deceased 124                  0           0
```

```r
# 5-year true relapses is even rare in data
i <- which(y$donor_interval_of_last_followup < 365*5)
length(i) # number of followup less than 5 years
```

```
## [1] 682
```

```r
table(y$disease_status_last_followup[i])
```

```
## 
##                    complete remission        progression 
##                 65                599                 18
```

# Session info


```r
sessionInfo()
```

```
## R version 3.4.3 (2017-11-30)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Sierra 10.12.6
## 
## Matrix products: default
## BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## loaded via a namespace (and not attached):
## [1] compiler_3.4.3  magrittr_1.5    tools_3.4.3     yaml_2.1.16    
## [5] stringi_1.1.6   knitr_1.18      stringr_1.2.0   evaluate_0.10.1
```
