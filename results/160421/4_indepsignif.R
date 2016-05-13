# this script is used to run double nested CV runs to later select predictors
# run this script on crom01 for prediction results

# read cluster parameters from bash command line
ags <- commandArgs(trailingOnly = TRUE)
stopifnot(length(ags) == 8)
xname <- as.character(ags[1]) # feature matrix
yname <- as.character(ags[2]) # response groups
i.fold <- as.integer(ags[3]) # cv fold index
nfolds <- as.integer(ags[4]) # number of folds
nrepeats <- as.integer(ags[5]) # number of repeats
pthres <- as.numeric(ags[6]) # p-value threshold
test <- as.character(ags[7]) # test method that chooses from c("t.test", "wilcox.test", ...)
method <- as.character(ags[8]) # multiple testing correction method, default should be "none"

# start! ------------------------------------------------------------------

message("\nRunning ...\n", paste(ags, collapse = "\t"),"\n")

objname <- paste0('plist_', paste(ags, collapse = '_'))
objpath <- paste0('Robj/', objname, '.RData')
if (file.exists(objpath)) {
  message('job already done !!')
  quit(save = 'no')
}


message("loading features and groups ...")
load('dat.RData')
source("../../src/func.R")

# get R objects
xtr <- get(xname)
ytr <- get(yname)
# keep only samples for which label info is available
samplelist <- intersect(rownames(xtr), names(ytr))
message('Sample size = ', length(samplelist))
xtr <- xtr[samplelist,]
ytr <- ytr[samplelist]

# rm other un-used objects
rm(list = ls(pattern = '[.](grps|vals|kmat)$'))
gc()

message("creating data split for CV ...")
set.seed(94151402)
foldIndices <- caret::createMultiFolds(1:nrow(xtr), k=nfolds, times=nrepeats)
fold <- foldIndices[[i.fold]]
# set train-test split
train.fold <- fold
test.fold <- setdiff(1:nrow(xtr), train.fold)


message('indep signif test ... ')

res <- featselectIndepSignif(xtr = xtr[train.fold, , drop=F], ytr = ytr[train.fold], 
                             pthres = pthres, test = test, method = method, 
                             keep.signif = FALSE) # return all features
assign(objname, res)
if (!dir.exists('Robj')) dir.create('Robj')
save(list = objname, file = objpath)
message('new job saved up !!')
quit(save = 'no')
