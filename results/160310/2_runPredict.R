# this script is used to run double nested CV runs to later select predictors
# run this script on crom01 for prediction results

# read cluster parameters from bash command line
ags <- commandArgs(trailingOnly = TRUE)
stopifnot(length(ags) == 9)
xname <- as.character(ags[1]) # feature matrix
yname <- as.character(ags[2]) # response groups
prname <- as.character(ags[3]) # predictor
i.fold <- as.integer(ags[4]) # cv fold index
nfolds <- as.integer(ags[5]) # number of folds
nrepeats <- as.integer(ags[6]) # number of repeats
i.fold.inn <- as.integer(ags[7]) # inner cv fold index
nfolds.inn <- as.integer(ags[8]) # inner number of folds
nrepeats.inn <- as.integer(ags[9]) # inner number of repeats

message("\nRunning ...\n", paste(ags, collapse = "\t"),"\n")

if (i.fold.inn == 0) {
  objname <- paste0('res_', paste(ags, collapse = '_'))
} else {
  objname <- paste0('cvres_', paste(ags, collapse = '_'))
}
objpath <- paste0('Robj/', objname, '.RData')
if (file.exists(objpath)) {
  message('job already done !!')
  quit(save = 'no')
}

# start! ------------------------------------------------------------------

message("loading features and groups ...")
load('dat.RData')
source("../../src/func.R") # complementary functions for kernel kmeans for top-k rankings

# get R objects
xtr <- get(xname)
ytr <- get(yname)
# rm other un-used objects
rm(list = ls(pattern = '[.](grps|vals)$'))
# keep only samples for which label info is available
samplelist <- intersect(rownames(xtr), names(ytr))
message('Sample size = ', length(samplelist))
xtr <- xtr[samplelist,]
ytr <- ytr[samplelist]

# create cv folds for model evaluation
set.seed(94151402)
foldIndices <- caret::createMultiFolds(1:nrow(xtr), k=nfolds, times=nrepeats)
fold <- foldIndices[[i.fold]]
# create inner cv folds to select predictor
set.seed(19817919)
foldIndices.inn <- caret::createMultiFolds(fold, k=nfolds.inn, times=nrepeats.inn)
# set train-test split
if (i.fold.inn == 0) {
  train.fold <- fold
  test.fold <- setdiff(1:nrow(xtr), train.fold)
} else {
  train.fold <- fold[foldIndices.inn[[i.fold.inn]]]
  test.fold <- setdiff(fold, train.fold)
}

message('train and predict ... ')
assign(objname, 
       indepValidation(xtr = xtr[train.fold, , drop=F], ytr = ytr[train.fold], 
                       xtst = xtr[test.fold, , drop=F], ytst = ytr[test.fold], 
                       predictor = prname))
if (!dir.exists('Robj')) dir.create('Robj')
save(list = objname, file = objpath)
message('new job saved up !!')
quit(save = 'no')
