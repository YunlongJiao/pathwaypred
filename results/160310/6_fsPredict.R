# this script is used to run double nested CV runs with a specific prediction manner
# 1) pathway features are pre-selected by logitLasso
# 2) build a logit model with pre-selected pathways always included in the model while genes 
#    that provide othogonal information are selected
# NOTE for computational reason, penalty.factor is set 1 to 1e-6 instead of 1 to 0 exactly;
#      see glmnet() for more details

# run this script on crom01 for prediction results

# read cluster parameters from bash command line
ags <- commandArgs(trailingOnly = TRUE)
stopifnot(length(ags) == 9)
xname <- as.character(ags[1]) # feature matrix
yname <- as.character(ags[2]) # response groups
prname <- as.character(ags[3]) # predictor # NOTE "prname" corresponds to "3"-rd entry DO NOT CHANGE FOR SAKE OF lamags
i.fold <- as.integer(ags[4]) # cv fold index
nfolds <- as.integer(ags[5]) # number of folds
nrepeats <- as.integer(ags[6]) # number of repeats
i.fold.inn <- as.integer(ags[7]) # inner cv fold index # NOTE "i.fold.inn" corresponds to "7"-th entry DO NOT CHANGE FOR SAKE OF lamags
nfolds.inn <- as.integer(ags[8]) # inner number of folds
nrepeats.inn <- as.integer(ags[9]) # inner number of repeats

lamags <- ags
lamags[3] <- "predictorLogitLasso"
lamags[7] <- as.integer(0) # for later paste into lam-name,lam-path to retrieve lambda list so that it always corresponds to i.fold.inn=0
cvags <- ags
cvags[7] <- "[[:digit:]]+"
sep <- "[.]fs[.]"
stopifnot(grepl(paste0("^.+[.]vals",sep,".+[.]vals$"), xname)) # NO-PEN features and PEN features are separated by "sep"
stopifnot(prname == "predictorLogitLasso2StepFS")
fsname <- "featselectLogitLasso" # for the final featselect

# start! ------------------------------------------------------------------

message("\nRunning ...\n", paste(ags, collapse = "\t"),"\n")

if (i.fold.inn == 0) {
  objname <- paste0('ivres_', paste(ags, collapse = '_'))
} else {
  objname <- paste0('cvres_', paste(ags, collapse = '_'))
}
objpath <- paste0('Robj/', objname, '.RData')
if (file.exists(objpath) && 
    (grepl('^cvres',objname) || 
     (grepl('^ivres',objname) && !inherits(objobj <- try(get(load(objpath))), "try-error") && !is.null(objobj$test.class)))) {
  message('job already done !!')
  quit(save = 'no')
}

message("loading features and groups ...")
load('dat.RData')
source("../../src/func.R")

# specific predictor function TO BE added to func.R file ...
predictorLogitLasso2StepFS <- function(xtr, xtst, ytr, alpha = 1, cutoff = 0.5, do.normalize = TRUE, 
                                       lam.nopen, lam.pen, # for specific 2-step training
                                       i.fold.inn, # deciding if it is CV run or not
                                       cv.patt, # for reading CV results
                                       penalty.factor.ratio.min = 1e-6, penalty.factor.ratio.max = 1e6, ...)
{
  # lam.nopen,lam.pen are two lists of lambda's to fit model with in either step respectively
  # i.fold.inn is set for cluster run which is first to run CV folds and later select lambdas and make prediction
  # NOTE predictorLogitLasso2StepFS() calls glmnet() for path run, instead of cv.glmnet() in predictorLogitLasso()
  # TODO in this implementation lambda's are selected by lambda.min (instead of lambda.1se) wrt CV classification accuracy
  
  library(glmnet)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  lam.nopen <- sort(lam.nopen, decreasing = TRUE)
  lam.pen <- sort(lam.pen, decreasing = TRUE)
  
  message("Training Step I ...")
  model.nopen <- glmnet(x = as.matrix(xtr), y = as.factor(as.character(ytr)), family = "multinomial", 
                        lambda = lam.nopen, 
                        alpha = alpha, standardize = do.normalize, ...)
  message("Feat selection Step I ...")
  featlist.nopen <- featselectLogitLasso2StepFS(model = model.nopen, keep.signif = FALSE) # fliter-out later
  
  message("Training Step II ...")
  penalty.factor <- rep(1, ncol(xtr))
  names(penalty.factor) <- colnames(xtr)
  model <- lapply(1:length(lam.nopen), function(i){
    message(".", appendLF = FALSE)
    featlist.nopen <- featlist.nopen[[i]]
    # no penalty on pre-selected features
    penalty.factor[names(featlist.nopen[featlist.nopen > 0])] <- penalty.factor.ratio.min
    # exclude other features
    penalty.factor[names(featlist.nopen[featlist.nopen == 0])] <- penalty.factor.ratio.max
    model.pen <- glmnet(x = as.matrix(xtr), y = as.factor(as.character(ytr)), family = "multinomial", 
                        lambda = lam.pen, 
                        alpha = alpha, standardize = do.normalize, ...)
    message("+", appendLF = FALSE)
    model.pen$featlist.short <- featselectLogitLasso2StepFS(model = model.pen, keep.signif = TRUE)
    return(model.pen)
  })
  if (is.null(names(lam.nopen)))
    names(lam.nopen) <- paste0("s",0:(length(lam.nopen)-1))
  names(model) <- names(lam.nopen)
  model$lam.nopen <- lam.nopen
  message(" done!")
  
  # should not return path run predictions if this is CV run
  pred.class <- NULL
  pred.prob <- NULL
  # combine CV results if this is NOT CV run
  cv.files <- list.files(path = 'Robj', pattern = cv.patt, full.names = TRUE)
  if (i.fold.inn == 0) {
    message("combining CV runs ...")
    # get CV results from previous cluster runs
    cv.res <- lapply(cv.files, function(f) try(get(load(f))))
    cv.res <- cv.res[sapply(cv.res, function(x) !inherits(x, "try-error"))]
    # get CV performance
    message("getting cv.perf ...")
    cv.perf <- lapply(cv.res, function(res){ # for each CV fold
      message(".", appendLF = FALSE)
      # get CV test samples
      cv.ytst <- res$true.class
      stopifnot(all(names(cv.ytst) %in% rownames(xtr)))
      # get CV trained models
      model.nopen <- res$model
      stopifnot(all(model.nopen$lam.nopen == lam.nopen))
      
      # get feat number
      cv.n.featlist.short <- lapply(model.nopen[names(model.nopen$lam.nopen)], function(model.pen){
        # get totally selected featlist for each lam.pen
        return(sapply(model.pen$featlist.short, length))
      })
      cv.n.featlist.short <- do.call('rbind', cv.n.featlist.short) # "rbind" leads to lam.nopen in row, lam.pen in col
      
      # get CV acc
      cv.acc <- lapply(model.nopen[names(model.nopen$lam.nopen)], function(model.pen){ # for each lam.nopen
        # predict for each lam.pen
        stopifnot(all(model.pen$lambda == lam.pen))
        cv.pred <- predict(object = model.pen, newx = as.matrix(xtr[names(cv.ytst), ]), type = "class", ...)
        return(apply(cv.pred, 2, function(pp) mean(pp == cv.ytst)))
      })
      cv.acc <- do.call('rbind', cv.acc) # "rbind" leads to lam.nopen in row, lam.pen in col
      
      return(list(cv.n.featlist.short=cv.n.featlist.short, cv.acc=cv.acc))
    })
    # average fold n.feat
    cv.n.featlist.short <- unlist(lapply(cv.perf, '[[', 'cv.n.featlist.short'))
    dim(cv.n.featlist.short) <- c(length(lam.nopen), length(lam.pen), length(cv.res))
    cv.n.featlist.short <- apply(cv.n.featlist.short, c(1,2), mean)
    # average fold scores
    cv.acc <- unlist(lapply(cv.perf, '[[', 'cv.acc'))
    dim(cv.acc) <- c(length(lam.nopen), length(lam.pen), length(cv.res))
    cv.acc <- apply(cv.acc, c(1,2), mean)
    # lambda selection rule!! because we set lam.nopen and lam.pen to be decreasing and cv.acc is created by "rbind", we now select lambda's confining to 
    #   1) larger CV acc; 
    #   2) then less total number of feat (lambda); 
    #   3) then more number of nopen feat 
    ind.acc <- which(cv.acc == max(cv.acc), arr.ind = TRUE)
    ind.acc.n.feat <- which(cv.n.featlist.short == min(cv.n.featlist.short[ind.acc]), arr.ind = TRUE)
    #   2) then less total number of feat (lambda); 
    #   3) then more number of nopen feat 
    ind.acc <- which(cv.acc == max(cv.acc), arr.ind = TRUE)
    ind.acc.n.feat <- which(cv.n.featlist.short == min(cv.n.featlist.short[ind.acc]), arr.ind = TRUE)
    ind <- ind.acc.n.feat[which.max(ind.acc.n.feat[ ,"row"]), ] # which.max returns the first hit
    best.lam.nopen <- lam.nopen[ind["row"]]
    best.lam.pen <- lam.pen[ind["col"]]
    
    # get trained model
    message("getting trained model ...")
    model <- model[[names(best.lam.nopen)]]
    model$best.lam.nopen <- best.lam.nopen
    model$best.lam.pen <- best.lam.pen
    model$featlist.short <- NULL # annul this previously registered entry
    # predict with seleted lambda's
    message("making prediction ...")
    pred <- drop(predict(object = model, newx = as.matrix(xtst), type = "response", 
                         s = best.lam.pen, ...))
    pred.class <- colnames(pred)[max.col(pred)]
    names(pred.class) <- rownames(xtst)
    pred.prob <- matrix(0, nrow = nrow(xtst), ncol = length(classes), 
                        dimnames = list(rownames(xtst), classes))
    pred.prob[ , colnames(pred)] <- pred
    message(" done done!")
  }
  
  res <- list(model=model, class=pred.class, prob=pred.prob, cutoff=cutoff)
  return(res)
}

featselectLogitLasso2StepFS <- function(model = NULL, ..., keep.signif = TRUE)
{
  # computes absolute value of coefficients for each feature
  # returns LIST of only those with non-zero contribution to at least one of the phenotypes ordered by decreasing sum contribution
  # NOTE model must be fit with glmnet() for path run and only makes sense when do.normalize set TRUE as we do sum contribution
  
  library(glmnet)
  
  if (is.null(model))
    model <- predictorLogitLasso2StepFS(...)$model
  s <- predict(object = model, newx = NULL, type = "coefficients", ...)
  dims <- c(nrow(s[[1]]), ncol(s[[1]]), length(s))
  dims.names <- list(rownames(s[[1]]), colnames(s[[1]]), names(s))
  s <- as.vector(do.call('cbind', s))
  dim(s) <- dims
  dimnames(s) <- dims.names
  fs <- lapply(1:dim(s)[2], function(j){
    u <- s[ ,j, ] # reduce to 2-dim matrix
    u <- rowSums(abs(u[-1, ]))
    if (keep.signif)
      u <- u[u > 0] # can be strongly positive or negative related!
    u <- sort(u, decreasing = TRUE)
    return(u)
  })
  names(fs) <- dimnames(s)[[2]]
  
  return(fs)
}


# load feature sets
xname.sp <- unlist(strsplit(xname, sep))
xname.nopen <- xname.sp[1]
xname.pen <- xname.sp[2]

# get lambda lists from result lists of previous predictions made on individual feature types
lampath <- paste0('Robj/ivres_', paste(lamags, collapse = '_'), '.RData')
lam.nopen <- get(load(gsub(xname, xname.nopen, lampath)))$model$lambda
lam.pen <- get(load(gsub(xname, xname.pen, lampath)))$model$lambda
# get R objects
xtr <- cbind(removeConst(get(xname.nopen)), removeConst(get(xname.pen))) # rownames were well aligned already
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
res <- indepValidation(xtr = xtr[train.fold, , drop=F], ytr = ytr[train.fold], 
                       xtst = xtr[test.fold, , drop=F], ytst = ytr[test.fold], 
                       predictor = prname, 
                       remove.const = FALSE, # constant already removed
                       lam.nopen = lam.nopen, lam.pen = lam.pen, i.fold.inn = i.fold.inn, cv.patt = paste0('^cvres_', paste(cvags, collapse = '_')))
message('feature selection ... ')
if (!is.null(res$model$best.lam.pen))
  res$featlist.short <- names(get(fsname, mode = "function")(model = res$model, 
                                                             s = res$model$best.lam.pen)) # get feats for only one lambda value

assign(objname, res)
