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

sep <- "[.]fs[.]"
stopifnot(grepl(paste0("^.+[.]vals",sep,".+[.]vals$"), xname)) # NO-PEN features and PEN features are separated by "sep"
stopifnot(prname == "predictorLogitLasso")

# start! ------------------------------------------------------------------

message("\nRunning ...\n", paste(ags, collapse = "\t"),"\n")

if (i.fold.inn == 0) {
  objname <- paste0('ivres_', paste(ags, collapse = '_'))
} else {
  objname <- paste0('cvres_', paste(ags, collapse = '_'))
}
objpath <- paste0('Robj/', objname, '.RData')
if (file.exists(objpath)) {
  message('job already done !!')
  quit(save = 'no')
}


message("loading features and groups ...")
load('dat.RData')
source("../../src/func.R")

# load feature sets
xname.sp <- unlist(strsplit(xname, sep))
xname.nopen <- xname.sp[1]
xname.pen <- xname.sp[2]

# get pre-selected features
featlist.nopen <- get(load(gsub(xname, xname.nopen, objpath)))[["featlist.short"]]
# get R objects
xtr <- cbind(get(xname.nopen)[,featlist.nopen], removeConst(get(xname.pen))) # rownames were well aligned already
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
penalty.factor <- rep(1, ncol(xtr))
if (length(featlist.nopen) > 0) 
  penalty.factor[1:length(featlist.nopen)] <- 0 # no penalty on pre-selected features
res <- indepValidation(xtr = xtr[train.fold, , drop=F], ytr = ytr[train.fold], 
                       xtst = xtr[test.fold, , drop=F], ytst = ytr[test.fold], 
                       predictor = prname, 
                       remove.const = FALSE, # constant already removed
                       penalty.factor = penalty.factor) # no penalty on pre-selected features
message('feature selection ... ')
fsname <- gsub("^predictor", "featselect", prname)
res$featlist.short <- names(get(fsname, mode = "function")(model = res$model))

assign(objname, res)
if (!dir.exists('Robj')) dir.create('Robj')
save(list = objname, file = objpath)
message('new job saved up !!')
quit(save = 'no')
