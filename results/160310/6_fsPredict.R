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
                                       featlist.nopen, # indices of feats to include in the 1-step fs
                                       i.fold.inn, # deciding if it is CV run or not
                                       cv.patt, # for reading CV results stored in folder Robj/
                                       penalty.factor.ratio = 10, # for 2-step penalty.factor
                                       lam.pen.ratio = penalty.factor.ratio*penalty.factor.ratio, ...) # multiplied to lam.pen to accommodate 2-step fs
{
  # lam.nopen,lam.pen are two lists of lambda's to fit model with in either step respectively
  # i.fold.inn is set for cluster run which is first to run CV folds and later select lambdas and make prediction
  # NOTE predictorLogitLasso2StepFS() calls glmnet() for path run, instead of cv.glmnet() in predictorLogitLasso()
  # TODO in this implementation lambda's are selected by lambda.min (instead of lambda.1se) wrt CV classification accuracy
  
  library(glmnet)
  
  classes <- names(tab <- table(ytr))
  classes.eff <- classes[tab > 0]
  stopifnot(length(classes.eff) >= 2)
  
  # sort and name lambda lists
  lam.nopen <- sort(lam.nopen, decreasing = TRUE)
  if (is.null(names(lam.nopen)))
    names(lam.nopen) <- paste0("s",0:(length(lam.nopen)-1))
  lam.pen <- sort(lam.pen*lam.pen.ratio, decreasing = TRUE)
  if (is.null(names(lam.pen)))
    names(lam.pen) <- paste0("s",0:(length(lam.pen)-1))
  
  message("Training Step I ...")
  model.nopen <- glmnet(x = as.matrix(xtr[ ,featlist.nopen,drop=FALSE]), y = as.factor(as.character(ytr)), family = "multinomial", 
                        lambda = lam.nopen, 
                        alpha = alpha, standardize = do.normalize, ...)
  message("Feat selection Step I ...")
  featlist.nopen <- featselectLogitLasso2StepFS(model = model.nopen, s = lam.nopen, keep.signif = FALSE) # fliter-out later
  
  message("Training Step II ...")
  penalty.factor <- rep(1, ncol(xtr))
  names(penalty.factor) <- colnames(xtr)
  model <- lapply(1:length(lam.nopen), function(i){
    message(".", appendLF = FALSE)
    featlist.nopen <- featlist.nopen[[i]]
    # no penalty on pre-selected features
    penalty.factor[names(featlist.nopen)[featlist.nopen > 0]] <- 1/penalty.factor.ratio
    # exclude other features
    penalty.factor[names(featlist.nopen)[featlist.nopen == 0]] <- penalty.factor.ratio
    model.pen <- try(glmnet(x = as.matrix(xtr), y = as.factor(as.character(ytr)), family = "multinomial", 
                            lambda = NULL, penalty.factor = penalty.factor,
                            alpha = alpha, standardize = do.normalize, ...))
    message("+", appendLF = FALSE)
    if (!inherits(model.pen, "try-error")) {
      if (min(model.pen$lambda) > max(lam.pen) || max(model.pen$lambda) < min(lam.pen))
        message("To retrain 2-step model.pen for lam.nopen = ", lam.nopen[i], "is necessary !")
      model.pen$featlist.short <- featselectLogitLasso2StepFS(model = model.pen, s = lam.pen, keep.signif = TRUE)
    }
    return(model.pen)
  })
  names(model) <- names(lam.nopen)
  model$lam.nopen <- lam.nopen
  model$lam.pen <- lam.pen
  message(" done!")
  
  # should not return path run predictions if this is CV run
  pred.class <- NULL
  pred.prob <- NULL
  # combine CV results if this is NOT CV run
  if (i.fold.inn == 0 && length(cv.files <- list.files(path = 'Robj', pattern = cv.patt, full.names = TRUE)) > 0) {
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
      cv.n.featlist.short <- lapply(model.nopen[names(model.nopen$lam.nopen)], function(model.pen){ # for each lam.nopen
        # get totally selected featlist for all lam.pen
        if (inherits(model.pen, "try-error"))
          return(rep(NA, length(lam.pen)))
        else
          return(sapply(model.pen$featlist.short, length))
      })
      cv.n.featlist.short <- do.call('rbind', cv.n.featlist.short) # "rbind" leads to lam.nopen in row, lam.pen in col
      
      # get CV acc
      cv.acc <- lapply(model.nopen[names(model.nopen$lam.nopen)], function(model.pen){ # for each lam.nopen
        # predict for all lam.pen
        if (inherits(model.pen, "try-error"))
          return(rep(NA, length(lam.pen)))
        else {
          cv.pred <- predict(object = model.pen, newx = as.matrix(xtr[names(cv.ytst), ]), s = lam.pen, type = "class", ...)
          return(apply(cv.pred, 2, function(pp) mean(pp == cv.ytst)))
        }
      })
      cv.acc <- do.call('rbind', cv.acc) # "rbind" leads to lam.nopen in row, lam.pen in col
      
      return(list(cv.n.featlist.short=cv.n.featlist.short, cv.acc=cv.acc))
    })
    # average fold n.feat
    cv.n.featlist.short <- unlist(lapply(cv.perf, '[[', 'cv.n.featlist.short'))
    dim(cv.n.featlist.short) <- c(length(lam.nopen), length(lam.pen), length(cv.res))
    cv.n.featlist.short <- apply(cv.n.featlist.short, c(1,2), mean, na.rm = TRUE)
    # average fold scores
    cv.acc <- unlist(lapply(cv.perf, '[[', 'cv.acc'))
    dim(cv.acc) <- c(length(lam.nopen), length(lam.pen), length(cv.res))
    cv.acc <- apply(cv.acc, c(1,2), mean, na.rm = TRUE)
    # lambda selection rule!! because we set lam.nopen and lam.pen to be decreasing and cv.acc is created by "rbind", we now select lambda's confining to 
    #   1) larger CV acc; 
    #   2) then less total number of feat (lambda); 
    #   3) then more number of nopen feat 
    ind.acc <- which(cv.acc == max(cv.acc), arr.ind = TRUE)
    ind.acc.n.feat <- which(cv.n.featlist.short == min(cv.n.featlist.short[ind.acc]), arr.ind = TRUE)
    ind <- ind.acc.n.feat[which.max(ind.acc.n.feat[ ,"row"]), ] # which.max returns the first hit
    best.lam.nopen <- lam.nopen[ind["row"]]
    best.lam.pen <- lam.pen[ind["col"]]
    message("FOUND best.lam.nopen = \t", best.lam.nopen, "\tbest.lam.pen = \t", best.lam.pen)
    
    # get trained model
    message("getting trained model ...")
    model <- model[[names(best.lam.nopen)]]
    stopifnot(!inherits(model, "try-error"))
    model$best.lam.nopen <- best.lam.nopen
    model$best.lam.pen <- best.lam.pen
    model$featlist.nopen <- model$featlist.short[[names(best.lam.pen)]]
    model$featlist.short <- NULL # annul this previously registered entry
    # predict with seleted lambda's
    message("making prediction ...")
    pred <- predict(object = model, newx = as.matrix(xtst), type = "response", s = best.lam.pen, ...)
    pred <- drop(pred)
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

featselectLogitLasso2StepFS <- function(model = NULL, s = model$lambda, ..., keep.signif = TRUE)
{
  # computes absolute value of coefficients for each feature
  # returns LIST of only those with non-zero contribution to at least one of the phenotypes ordered by decreasing sum contribution
  # NOTE model must be fit with glmnet() for path run and only makes sense when do.normalize set TRUE as we do sum contribution
  
  library(glmnet)
  
  if (is.null(model))
    model <- predictorLogitLasso2StepFS(...)$model
  coefs <- predict(object = model, newx = NULL, s = s, type = "coefficients", ...)
  dims <- c(nrow(coefs[[1]]), ncol(coefs[[1]]), length(coefs))
  dims.names <- list(rownames(coefs[[1]]), colnames(coefs[[1]]), names(coefs))
  coefs <- as.vector(do.call('cbind', coefs))
  dim(coefs) <- dims
  dimnames(coefs) <- dims.names
  feats <- lapply(1:dim(coefs)[2], function(j){
    u <- coefs[ ,j, ] # reduce to 2-dim matrix
    u <- rowSums(abs(u[-1, ]))
    if (keep.signif)
      u <- u[u > 0] # can be strongly positive or negative related!
    u <- sort(u, decreasing = TRUE)
    return(u)
  })
  names(feats) <- names(s)
  
  return(feats)
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
xtr.nopen <- removeConst(get(xname.nopen))
xtr.pen <- removeConst(get(xname.pen))
xtr <- cbind(xtr.nopen, xtr.pen) # rownames were well aligned already
featlist.nopen <- 1:ncol(xtr.nopen) # indices to include in 1-step fs
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
                       lam.nopen = lam.nopen, featlist.nopen = featlist.nopen, 
                       lam.pen = lam.pen, 
                       i.fold.inn = i.fold.inn, cv.patt = paste0('^cvres_', paste(cvags, collapse = '_')))
message('feature selection ... ')
if (!is.null(res$model$best.lam.pen))
  res$featlist.short <- names(get(fsname, mode = "function")(model = res$model, 
                                                             s = res$model$best.lam.pen)) # get feats for only one lambda value

assign(objname, res)
if (!dir.exists('Robj')) dir.create('Robj')
save(list = objname, file = objpath)
message('new job saved up !!')
quit(save = 'no')

