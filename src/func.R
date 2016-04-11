# # generic libraries needed for general purpose
# # from CRAN
# library(knitr)
# library(devtools)
# library(caret)
# library(ROCR)
# # from Bioconductor
# library(survcomp)
# # already exists for 3.2.1-atalas on crom01
# library(survival)
# library(ggplot2)

# # import classifier libraries in those predictors called
# # from CRAN
# library(glmnet) # L1-GLM
# library(LiblineaR) # L1reg-L2loss-SVM
# library(gbm) # GBM
# library(kernlab) # estimate gamma for Gaussian RBF SVM
# # already exists for 3.2.1-atalas on crom01
# library(MASS) # LDA
# library(e1071) # NB, SVM and KNN tuning
# library(class) # KNN predicting
# library(randomForest) # RF
# library(rpart) # DT # tooooo many errors T T
# library(nnet) # NNet # tooooo many errors T T

library(methods) # manually import important BASE packages otherwise error in calling glmnet()
options(stringsAsFactors = FALSE)


# data processing ---------------------------------------------------------


normalizeData <- function (xtr, xtst=NULL, do.center=TRUE, do.scale=TRUE, ...)
{
	# normalizeData() normalizes col-wisely to mean 0 (do.center=T) and sd 1 (do.scale=T)
	# NOTE: xtst is normalized wrt xtr
	
	if(is.vector(xtr)){
		warning("vector features are coerced to a matrix!")
		xtr <- matrix(xtr, ncol=1)
		if(!is.null(xtst)){
			xtst <- matrix(xtst, ncol=1)
		}
	}
	
	numid <- sapply(as.data.frame(xtr), is.numeric)
	t <- scale(xtr[,numid], center=do.center, scale=do.scale)
	xtr[,numid] <- t
	
	if(is.null(xtst)){
		return(xtr)
	} else{
		if(do.center){
			xtst[,numid] <- sweep(xtst[,numid], 2, attr(t,"scaled:center"), FUN="-")
		}
		if(do.scale){
			xtst[,numid] <- sweep(xtst[,numid], 2, attr(t,"scaled:scale"), FUN="/")
		}
		return(list(xtr=xtr, xtst=xtst))
	}
}

removeConst <- function (xtr, xtst=NULL, tol=1e-6)
{
  numid <- which(sapply(as.data.frame(xtr), is.numeric))
  constid <- apply(xtr[,numid], 2, function(u){
    (mean(u, na.rm = TRUE) < tol && diff(range(u, na.rm = TRUE)) < tol) || (diff(range(u, na.rm = TRUE)) / mean(u, na.rm = TRUE) < tol)
  })
  xtr <- xtr[ , numid[!constid], drop = F]
  if(is.null(xtst)){
    return(xtr)
  } else{
    xtst <- xtst[ , numid[!constid], drop = F]
    return(list(xtr=xtr, xtst=xtst))
  }
}


# validation --------------------------------------------------------------


evaluatePred <- function(pred, ytst, ysurv = NULL, pos.label = names(table(ytst))[2], savepath = NULL, beta = 1, ...)
{
	# pred results from calling predictorXXX(), ess. a list containing $class and $prob of which class must be factor or character
	# ytst is true binary class indicators
	# ysurv is of class Surv containing true surv.time and surv.event, defaulted unavailable
	# NOTE large pred$prob should be associated with short survival time!!
	# NOTE F(beta=1) is used by default while F(beta=0.5) can emphasize on precision more than recall, i.e. reduce FPR
  
  library(ROCR)
  library(survival)
  library(survcomp)
  
	if(length(pred$class) != length(ytst)){
		stop("Inconsistent prediction length!")
	}
	
	if(NA %in% ytst){
		id <- is.na(ytst)
		ytst <- ytst[!id]
		pred$class <- pred$class[!id]
		pred$prob <- pred$prob[!id]
	}
	
	if(length(unique(ytst)) > 2){
		warning("Multi-class true response in test vector!")
	} else if(length(unique(ytst)) == 1){
		warning("Single-class true response in test vector!")
	}
	
	if(is.factor(pred$class) || is.character(pred$class)){
		# pos.label is the char for positive label
		posid <- ytst==pos.label
		# quantitative criteria
		acc <- mean(pred$class==ytst)
		ppv <- sum(pred$class==pos.label & ytst==pos.label)/sum(pred$class==pos.label)
		fpr <- sum(pred$class==pos.label & ytst!=pos.label)/sum(ytst!=pos.label)
		tpr <- sum(pred$class==pos.label & ytst==pos.label)/sum(ytst==pos.label)
		f_measure <- function(b, precision, recall){
			return( (1+b*b)*precision*recall/(b*b*precision+recall) )
		}
		fval <- f_measure(b=beta,precision=ppv,recall=tpr)
		
		# concordance index for survival risk prediction
		if(is.null(ysurv)){
			ci <- NA
		} else{
			ci <- concordance.index(pred$prob, surv.time=ysurv[,1], surv.event=ysurv[,2])$c.index
		}
		
		# auroc
		predROCR <- prediction(pred$prob, ytst)
		auroc <- as.numeric(performance(predROCR, "auc")@y.values)
		if(length(auroc)>1){
			warning("Average AUC score computed!")
			auroc <- mean(auroc, na.rm=TRUE)
		}
		
		# ROC curve
		if(!is.null(savepath)){
		  if (!grepl("[.]pdf$", savepath)) savepath <- paste0(savepath, "_ROC_.pdf")
			plotROCcv(res=list(test.prob=pred$prob, true.class=ytst), savepath=savepath)
		}
		
		# return results
		return(list(acc=acc, fpr=fpr, tpr=tpr, ppv=ppv, fval=fval, concordance.index=ci, auroc=auroc))
	} else{
		stop("Indefinite predicted response type!")
	}
}


indepValidation <- function(xtr, ytr, xtst, ytst, predictor, ysurv = NULL, 
                            remove.const = TRUE, 
                            pthres = 0, method = "wilcox.test", 
                            ..., save.model = FALSE, seed = 53359292)
{
	# xtr, xtst are n*p feature matrices and ytr, ytst are label vectors
	# predictor is a char denoting the predictor function name (starting with predictorxxx)
  # remove.const defaulted TRUE to remove constant feature columns in xtr
  # pthres defaulted 0 does no indep signif feature reduction
	# indepValidation() returns model performances validated against test set which have been trained on training set
  
  library(survival)
  
  if (!is.null(seed)) set.seed(seed)
  
  classes <- sort(unique(as.character(ytr)))
  neg.label <- classes[1]; pos.label <- classes[2]
	
	if(!is.character(predictor)){
		stop("Specify predictor with a character!")
	}
	
	if (!is.null(dimnames(xtr)) && !is.null(dimnames(xtst)) && !is.null(names(ytr)) && !is.null(names(ytst))) {
	  # adjust correp rows and cols position in case they do not match
	  ytr <- ytr[rownames(xtr)]
	  ytst <- ytst[rownames(xtst)]
	  xtst <- xtst[,colnames(xtr),drop=F]
	}
	
  # remove constant features
  if (remove.const) {
    d <- removeConst(xtr=xtr, xtst=xtst)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
	# indep signif test for feature reduction
	if(pthres > 0){
		plist <- indepSignif(xtr=xtr, ytr=ytr, method = method, ...)
		featlist <- names(plist)[plist < pthres]
		if(length(featlist) == 0){
			warning("No features survived after independence significance test hence no feature reduction is reserved!")
		  featlist <- colnames(xtr)
		}
		xtr <- xtr[,featlist,drop=F]
		xtst <- xtst[,featlist,drop=F]
	} else{
		featlist <- colnames(xtr)
	}
	
	# model and predict
  FUN <- get(predictor, mode = "function")
	pred <- FUN(xtr=xtr, xtst=xtst, ytr=ytr, ...)
	
	# evaluation and result list
	res <- c(list(validation="seqVal", predictor=predictor, cutoff=pred$cutoff, featlist=featlist, 
	              true.surv=ysurv, true.class=ytst, test.class=pred$class, test.prob=pred$prob), 
	         evaluatePred(pred=pred, ytst=ytst, ysurv=ysurv, pos.label = pos.label, ...))
	
	if(save.model){
	  res$model <- pred$model
	}
	
	return(res)
}



crossValidation <- function(xtr, ytr, ..., seed=61215942, nfolds=5, nrepeats=10)
{
	# nfolds, nrepeats are used to generate nrepeats simultaneous training splits each being nfolds
	# crossValidation() returns cross validated results from a single dataset
  
  library(caret)
  
	# data split
  if (!is.null(seed)) set.seed(seed)
	foldIndices <- createMultiFolds(1:nrow(xtr), k=nfolds, times=nrepeats)
	
	message(nrepeats," repeated experiments of ",nfolds, "-fold cross validation!")
	foldres <- lapply(foldIndices, function(fold){
		message("New fold started...")
		ivres <- indepValidation(xtr=xtr[fold,,drop=F], ytr=ytr[fold], xtst=xtr[-fold,,drop=F], ytst=ytr[-fold], ..., save.model=FALSE)
		return(ivres)
	})
	
	cvres <- crossValidationCombineResults(foldres)
	
	return(cvres)
}



crossValidationCombineResults <- function(foldres, 
                                          elem_ave=c("acc","fpr","tpr","ppv","fval","concordance.index","auroc"), 
                                          elem_keepfold=c("true.surv","true.class","test.class","test.prob","featlist"),
                                          elem_inherit=c("predictor", "cutoff"))
{
  # elem_ave are those scores to average over cv fold
  # elem_keepfold are those to keep as they come from each fold
  # elem_inherit are those to inherit some parameters from res of subfold as in ensembleResults()
  
  if (is.null(foldres) || (is.list(foldres) && length(foldres) == 0))
    return(list())
  
  cvres <- list(validation=paste0(nrepeats,"_repeats_",nfolds,"_fold_cross_validation"))
  cvres <- c(cvres,foldres[[1]][elem_inherit])
  for(elem in elem_keepfold){
    cvres[[elem]] <- lapply(foldres, function(u){u[[elem]]})
  }
  for(elem in elem_ave){
    cvres[[elem]] <- mean(sapply(foldres, function(u){u[[elem]]}), na.rm=TRUE)
  }
  
  return(cvres)
}



plotROCcv <- function(res, savepath, beta = 1)
{
	# res is the result list from calling indepValidation() or crossValidation(), ess having as elements list of $true.class and $test.prob (can be cross-fold result given by list)
  
  library(ROCR)
  
  if (is.null(res) || (is.list(res) && length(res) == 0))
    return(NULL)
  
	alpha <- function(b){return(1/(1+b*b))}
	
	if (!grepl("[.]pdf$", savepath)) savepath <- paste0(savepath, "_ROC_.pdf")
	
	if (!is.null(savepath)) pdf(savepath, height=10, width=15)
	
	par(mfrow=c(2,2))
	predROCR <- prediction(res$test.prob, res$true.class)
	
	acccurve <- performance(predROCR, "acc")
	plot(acccurve, lwd=2, main="Classification accuracy")
	grid()
	
	fcurve <- performance(predROCR, "f", alpha=alpha(b=beta))
	plot(fcurve, lwd=2, main=paste0("F-measure (beta=",beta,")"))
	grid()
	
	roccurve <- performance(predROCR, "tpr", "fpr")
	plot(roccurve, lwd=2, main="ROC curve")
	grid()
	
# 	tprcurve <- performance(predROCR, "tpr")
# 	plot(tprcurve, lwd=2, main="True Positive Rate (Recall or Sensitivity)")
# 	grid()
	
# 	fprcurve <- performance(predROCR, "fpr")
# 	plot(fprcurve, lwd=2, main="False Positive Rate (Fallout)")
# 	grid()
	
# 	liftcurve <- performance(predROCR, "lift", "rpp")
# 	plot(liftcurve, lwd=2, main="Lift Chart")
# 	grid()
	
# 	ppvcurve <- performance(predROCR, "ppv")
# 	plot(ppvcurve, lwd=2, main="Positive Predictive Value (Precision)")
# 	grid()
	
	prcurve <- performance(predROCR, "prec", "rec")
	plot(prcurve, lwd=2, main="Precision/Recall Graph")
	grid()
	
	if (!is.null(savepath)) dev.off()
	
	return(predROCR)
}


# feature reduction -------------------------------------------------------


indepSignif <- function(xtr, ytr, method = "wilcox.test", ...)
{
	# indepSignif() runs feature-by-feature independence significance test between contrasting classes
	# xtr should be assigned col names and ytr is a vector of labels
	# method denotes a char from c("t.test","wilcox.test",...)
	
	featlist <- colnames(xtr)
	plist <- sapply(featlist, function(featname){
		fo <- paste0(featname," ~ y")
		t <- get(method)(as.formula(fo), data=data.frame(y=as.factor(ytr),as.data.frame(xtr)), ...)
		return(t$p.value)
	})
	
	return(plist)
}


# classifiers -------------------------------------------------------------
# in this section
# xtr/xtst have ntr/ntst samples in rows, p explanatory variables in cols
# ytr have ntr (binary) indicators


predictorLDA <- function(xtr, xtst, ytr, cutoff = 0, do.normalize = TRUE, ...){
  # no feature selection, no parameter tuning
	
  library(MASS)
  
	classes <- sort(unique(as.character(ytr)))
	
	if(do.normalize){
		d <- normalizeData(xtr, xtst, ...)
		xtr <- d$xtr
		xtst <- d$xtst
	}
	
	model <- lda(xtr, as.factor(ytr), ...)
	
	pred <- predict(model, xtst)
	pred_prob <- pred$x[,1]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



# predictorLogit <- function(xtr, xtst, ytr, do.stepAIC = FALSE, cutoff = 0, do.normalize = TRUE, ...){
#   	# MANY BUGS FOR LARGE FEATURE MATRIX ...
# 	
# 	classes <- sort(unique(as.character(ytr)))
# 	if(length(classes) < 2){
# 		warning("Singe-class for training set!")
# 		res <- list(model="Single-class classifier", class=rep(classes,nrow(xtst)), prob=rep(1,nrow(xtst)))
# 		return(res)
# 	} else if (length(classes) > 2){
# 		stop("Multi-class setting...")
# 	}
# 	
# 	if(do.normalize){
# 		d <- normalizeData(xtr, xtst, ...)
# 		xtr <- d$xtr
# 		xtst <- d$xtst
# 	}
# 	
# 	model <- glm(y~., family=binomial(link="logit"), data=data.frame(y=as.factor(ytr),as.as.data.frame(xtr)), ...)
# 	if(do.stepAIC){
# 		message("Choosing a model by stepwise AIC...")
# 		model <- stepAIC(model, trace=0)
# 	}
# 	
# 	pred_prob <- predict(model, as.data.frame(xtst), type="link")
# 	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
# 	
# 	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
# 	return(res)
# }



predictorLogitLasso <- function(xtr, xtst, ytr, cutoff = 0, do.normalize = TRUE, ...){
	# predictorLogit() with Lasso for feature selection with parameter tuning
	
	library(glmnet)
  
	classes <- sort(unique(as.character(ytr)))
	if(length(classes) < 2){
		warning("Singe-class for training set!")
		res <- list(model="Single-class classifier", class=rep(classes,nrow(xtst)), prob=rep(1,nrow(xtst)))
		return(res)
	} else if (length(classes) > 2){
		stop("Multi-class setting...")
	}
	
	model <- cv.glmnet(as.matrix(xtr), ytr, family="binomial", standardize = do.normalize, ...)
	
	pred_prob <- predict(model, as.matrix(xtst), type="link")[,1]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorLinearSVM <- function(xtr, xtst, ytr, kernel = "linear", cost = 10^(-3:3), cutoff = 0, do.normalize = TRUE, ...){
	# cost is a vector of C parameter range to tune with, no feature selection
	# NOTE cross validation is set to minimize misclassification error by default of tune.svm arguments
	
	library(e1071)
  
	classes <- sort(unique(as.character(ytr)))
	
	o <- order(ytr, decreasing = TRUE) # otherwise decision values might be inversed
	model <- tune.svm(xtr[o,], as.factor(ytr[o]), kernel=kernel, cost=cost, scale=do.normalize, type="C-classification", tunecontrol=tune.control(sampling="cross", cross=5), ...)
	
	pred <- predict(model$best.model, xtst, decision.values=T)
	pred_prob <- attr(pred,"decision.values")[,1]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model$best.model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorRadialSVM <- function(xtr, xtst, ytr, do.normalize = TRUE, ...){
  # predictorLinearSVM with linear kernel replaced by Gaussian RBF kernel
  # gamma for radial kernel is set with median trick implemented by kernlab:sigest()
  
  gamma <- kernlab::sigest(xtr, scaled=do.normalize)['50%']
  res <- predictorLinearSVM(xtr = xtr, xtst = xtst, ytr = ytr, kernel = "radial", gamma = gamma, do.normalize = do.normalize, ...)
  return(res)
}



predictorKNN <- function(xtr, xtst, ytr, k = seq(1,10,2), do.normalize = TRUE, cutoff = 0.5, ...){
  # no feature selection, with parameter tuning for k
  # NOTE k takes only odd values to avoid confusion from ties
  # NOTE cross validation is set to minimize misclassification error by default of tune.knn arguments
	
	library(class)
	library(e1071)
	
	classes <- sort(unique(as.character(ytr)))
	
	if(do.normalize){
		d <- normalizeData(xtr, xtst, ...)
		xtr <- d$xtr
		xtst <- d$xtst
	}
	
	model <- tune.knn(xtr, ytr, k = k, tunecontrol=tune.control(sampling="cross", cross=5), ...)
	k <- model$best.parameters$k
	pred <- knn(xtr, xtst, ytr, k = k, prob = TRUE, ...)
	pred_prob <- attr(pred,"prob"); pred_prob[pred == classes[1]] <- 1 - pred_prob[pred == classes[1]]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=paste0(k,"-nearest neighbour classification"), class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorNB <- function(xtr, xtst, ytr, do.normalize = TRUE, cutoff = 0.5, ...){
  # no feature selection, no parameter tuning
  
  library(e1071)
  
  classes <- sort(unique(as.character(ytr)))
  
  if(do.normalize){
    d <- normalizeData(xtr, xtst, ...)
    xtr <- d$xtr
    xtst <- d$xtst
  }
  
  model <- naiveBayes(xtr, ytr, ...)
  pred_prob <- predict(model, xtst, type = "raw", ...)[,classes[2]]
  pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
  
  res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
  return(res)
}



predictorGBM <- function(xtr, xtst, ytr, n.trees = 1500, shrinkage = 0.002, interaction.depth = 6, bag.fraction = 1, cutoff = 0, ...){
  # no feature selection, no parameter tuning
	
	library(gbm)
  
	classes <- sort(unique(as.character(ytr)))
	ytr <- as.numeric(as.factor(ytr))-1
	
	model <- gbm.fit(xtr, ytr, distribution="bernoulli", interaction.depth=interaction.depth, n.trees=n.trees, shrinkage=shrinkage, bag.fraction=bag.fraction, verbose=FALSE, ...)
	
	pred_prob <- predict(model, xtst, model$n.trees, type="link")
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorRF <- function(xtr, xtst, ytr, ntrees = 500, cutoff = 0.5, ...){
  # no feature selection, no parameter tuning
	
	library(randomForest)
  
	classes <- sort(unique(as.character(ytr)))
	
	model <- randomForest(xtr, as.factor(ytr), ntree=ntrees, importance=TRUE, proximity=TRUE, do.trace=FALSE, ...)
	
	pred_prob <- predict(model, xtst, type="prob")[,classes[2]]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



# predictorDT <- function(xtr, xtst, ytr, cutoff = 0.5, do.prune = FALSE, ...){
# 	# TOO MANY BUGS ...
# 
# 	library(rpart)
#   
# 	classes <- sort(unique(as.character(ytr)))
# 	
# 	model <- rpart(y ~ ., data=data.frame(y=ytr,as.data.frame(xtr)), method="class", ...)
# 	if(do.prune){
# 		message("Pruning trees...")
# 		model <- prune(model, cp=model$cptable[which.min(model$cptable[,"xerror"]),"CP"])
# 	}
# 	
# 	pred_prob <- predict(model, as.data.frame(xtst), type="prob")[,classes[2]]
# 	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
# 	
# 	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
# 	return(res)
# }



# predictorNNet <- function(xtr, xtst, ytr, size = 2, cutoff = 0.5, do.normalize = TRUE, ...){
# 	# TOO MANY BUGS ...
# 	
# 	library(nnet)
#   
# 	classes <- sort(unique(as.character(ytr)))
# 	targets <- class.ind(as.factor(ytr))
# 	
# 	if(do.normalize){
# 		d <- normalizeData(xtr, xtst, ...)
# 		xtr <- d$xtr
# 		xtst <- d$xtst
# 	}
# 	
# 	model <- nnet(xtr, targets, size=size, ...)
# 	
# 	pred_prob <- predict(model, xtst, type="raw")[,classes[2]]
# 	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
# 	
# 	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
# 	return(res)
# }



predictorSparseSVM <- function(xtr, xtst, ytr, cost = 1, cutoff = 0, do.normalize = TRUE, ...){
  # cost is a constant C parameter so no parameter tuning
  # L1-regularized L2-loss support vector classification (i.e. with feature selection)
	
	library(LiblineaR)
  
	classes <- sort(unique(as.character(ytr)))
	
	if(do.normalize){
		d <- normalizeData(xtr, xtst, ...)
		xtr <- d$xtr
		xtst <- d$xtst
	}
	
	o <- order(ytr, decreasing=TRUE)
	model <- LiblineaR(xtr[o,], ytr[o], type=5, cost=cost, ...)
	
	pred <- predict(model, xtst, decisionValues = TRUE)
	pred_prob <- pred$decisionValues[,1] # though col 1 stands for classes[2]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}


# utiles ------------------------------------------------------------------


# NA proportion function
NArate <- function(v) mean(is.na(v))



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
# Taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


# Taken from http://www.r-bloggers.com/memory-management-in-r-a-few-tips-and-tricks/
# improved list of objects

.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}

# shorthand
lsos <- function(..., n=10) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}


