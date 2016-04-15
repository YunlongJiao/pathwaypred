# this script is used to run outter loop of CV runs to evaluate model performance
# run this script on crom01 for prediction results

# read cluster parameters from bash command line
ags <- commandArgs(trailingOnly = TRUE)
stopifnot(length(ags) == 6)
xname <- as.character(ags[1]) # feature matrix
yname <- as.character(ags[2]) # response groups
prname <- as.character(ags[3]) # predictor
i.fold <- as.integer(ags[4]) # cv fold index
nfolds <- as.integer(ags[5]) # number of folds
nrepeats <- as.integer(ags[6]) # number of repeats

message("\nRunning ...\n", paste(ags, collapse = "\t"),"\n")

objname <- paste(ags, collapse = '_')
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

# create cv folds
set.seed(94151402)
foldIndices <- caret::createMultiFolds(1:nrow(xtr), k=nfolds, times=nrepeats)
train.fold <- foldIndices[[i.fold]]
test.fold <- setdiff(1:nrow(xtr), train.fold)

message('train and predict ... ')
assign(objname, 
       indepValidation(xtr = xtr[train.fold, , drop=F], ytr = ytr[train.fold], 
                       xtst = xtr[test.fold, , drop=F], ytst = ytr[test.fold], 
                       predictor = prname))
if (!dir.exists('Robj')) dir.create('Robj')
save(list = objname, file = objpath)
message('new job saved up !!')
quit(save = 'no')
