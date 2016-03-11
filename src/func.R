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
# library(LiblineaR) # L1-SVM
# library(gbm) # GBM
# # already exists for 3.2.1-atalas on crom01
# library(MASS) # LDA, logit
# library(e1071) # SVM
# library(class) # KNN
# library(randomForest) # RF
# library(rpart) # DT
# library(nnet) # NNet

options(stringsAsFactors = FALSE)

########################
###
### data processing
###
########################

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
	
	numid <- sapply(as.data.frame(xtr),is.numeric)
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



########################
###
### validation
###
########################

evaluatePred <- function(pred, ytst, ysurv = NULL, pos.label = names(table(ytst))[2], savepath = NULL, beta = 0.5, ...)
{
	# pred results from calling predictorXXX(), ess. a list containing $class and $prob of which class must be factor or character
	# ytst is true binary class indicators
	# ysurv is of class Surv containing true surv.time and surv.event, defaulted unavailable
	# NOTE large pred$prob should be associated with short survival time!!
	# NOTE F(beta=0.5) or equiv F(alpha=0.8) is used by default since we emphasize on precision more than recall, i.e. reduce FPR
  
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
                            pthres = 0, method = "wilcox.test", 
                            ..., save.model = FALSE, seed = 53359292)
{
	# xtr, xtst are n*p feature matrices and ytr, ytst are label vectors
	# predictor is a char denoting the predictor function name (starting with predictorxxx)
  # pthres defaulted 0 does no indep signif feature reduction
	# indepValidation() returns model performances validated against test set which have been trained on training set
  
  library(survival)
  
  classes <- sort(unique(as.character(ytr)))
  neg.label <- classes[1]; pos.label <- classes[2]
	
	if (!is.null(seed)) set.seed(seed)
	if(!is.character(predictor)){
		stop("Specify predictor with a character!")
	}
	
	if (!is.null(dimnames(xtr)) && !is.null(dimnames(xtst)) && !is.null(names(ytr)) && !is.null(names(ytst))) {
	  # adjust correp rows and cols in case they do not match
	  ytr <- ytr[rownames(xtr)]
	  ytst <- ytst[rownames(xtst)]
	  xtst <- xtst[,colnames(xtr),drop=F]
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



crossValidation <- function(xtr, ytr, ..., seed=61215942, nfolds=5, nrepeats=10, 
                            elem_ave=c("acc","fpr","tpr","ppv","fval","concordance.index","auroc"), 
                            elem_keepfold=c("true.surv","true.class","test.class","test.prob","featlist"),
                            elem_inherit=c("predictor", "cutoff"))
{
	# nfolds, nrepeats are used to generate nrepeats simultaneous training splits each being nfolds
	# elem_ave are those scores to average over cv fold
	# elem_keepfold are those to keep as they come from each fold
	# elem_inherit are those to inherit some parameters from res of subfold as in ensembleResults()
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



plotROCcv <- function(res, savepath)
{
	# res is the result list from calling indepValidation() or crossValidation(), ess having as elements list of $true.class and $test.prob (can be cross-fold result given by list)
  
  library(ROCR)
  
	alpha <- function(b){return(1/(1+b*b))}
	
	if (!is.null(savepath)) pdf(savepath, height=10, width=15)
	
	par(mfrow=c(2,3))
	predROCR <- prediction(res$test.prob, res$true.class)
	
	roccurve <- performance(predROCR, "tpr", "fpr")
	plot(roccurve, lwd=2, main="ROC curve")
	grid()
	tprcurve <- performance(predROCR, "tpr")
	plot(tprcurve, lwd=2, main="True Positive Rate (Recall or Sensitivity)")
	grid()
	fprcurve <- performance(predROCR, "fpr")
	plot(fprcurve, lwd=2, main="False Positive Rate (Fallout)")
	grid()
	
	# liftcurve <- performance(predROCR, "lift", "rpp")
	# plot(liftcurve, lwd=2, main="Lift Chart")
	# grid()
	prcurve <- performance(predROCR, "prec", "rec")
	plot(prcurve, lwd=2, main="Precision/Recall Graph")
	grid()
	fcurve <- performance(predROCR, "f", alpha=alpha(b=0.5))
	plot(fcurve, lwd=2, main="F-measure (beta=0.5)")
	grid()
	ppvcurve <- performance(predROCR, "ppv")
	plot(ppvcurve, lwd=2, main="Positive Predictive Value (Precision)")
	grid()
	
	if (!is.null(savepath)) dev.off()
	
	return(predROCR)
}



########################
###
### feature reduction
###
########################

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



# featIllustrate_static <- function(xtr, ytr, method = "wilcox.test", featlist=NULL, pthres = 0.05, savepath = NULL, cols=3, ...)
# {
# 	# featlist denotes the idx or names of col to illustrate differences
# 	# savepath needs to end in ".pdf"
# 	# featIllustrate_static() defaults to run indepSignif() with inter-class histograms illustrating difference
# 	
#   
#   library(ggplot2)
#   
# 	classes <- sort(unique(as.character(ytr)))
# 	
# 	if(is.null(featlist)){
# 		plist <- indepSignif(xtr=xtr, ytr=ytr, method = method, ...)
# 		featlist <- names(plist)[plist < pthres]
# 		if(length(featlist) == 0){
# 			warning("No feature survived after independence significance test!")
#       plist <- NA
#       featlist <- colnames(xtr)
# 		}
# 	} else{
# 		plist <- NA
# 		featlist <- colnames(xtr[,featlist,drop=F])
# 	}
# 	
# 	plots <- list()
# 	for(featname in featlist){
# 		d <- data.frame(feat=xtr[[featname]], grp=ytr)
# 		grp_mean <- tapply(d$feat, d$grp, mean, na.rm=TRUE)
# 		grp_mean <- data.frame(grp=names(grp_mean), m=grp_mean)
# 		
# 		p1 <- ggplot(d) + xlab(featname) + 
# 			geom_density(aes(x=feat, color=grp), size=1) + 
# 			geom_histogram(aes(x=feat, fill=grp, y=..density..), alpha=0.5, position="identity") + 
# 			guides(color=FALSE, fill=guide_legend(title=NULL)) + 
# 			geom_vline(data=grp_mean, aes(xintercept=m, color=grp), linetype="dashed", size=1)
# 		plots[[featname]] <- p1
# 	}
# 	
# 	if(is.null(savepath)){
# 		return(plots)
# 	} else{
#     if (length(featlist) < 20) {
#       if (!grepl("[.]pdf$", savepath)) savepath <- paste0(savepath, "_features.pdf")
#       nw <- min(cols, length(plots))
#       nh <- ceiling(length(plots)/nw)
#       pdf(file=savepath, width=8*nw, height=5*nh)
#       multiplot(plotlist=plots, cols=cols)
#       dev.off()
#     } else {
#       for (featname in featlist) {
#         savepath <- paste(savepath, featname, ".pdf", sep = "_")
#         pdf(file=savepath, width=8, height=5)
#         plot(plots[[featname]])
#         dev.off()
#       }
#     }
# 		return(plist)
# 	}
# }



########################
###
### classifiers
###
########################

### in this section
### xtr/xtst have ntr/ntst samples in rows, p explanatory variables in cols
### ytr have ntr (binary) indicators

predictorLDA <- function(xtr, xtst, ytr, cutoff = 0, do.normalize = TRUE, ...){
	
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



predictorLogit <- function(xtr, xtst, ytr, do.stepAIC = FALSE, cutoff = 0, do.normalize = TRUE, ...){
	
	library(MASS)
  
	classes <- sort(unique(as.character(ytr)))
	if(length(classes) < 2){
		warning("Singe-class for training set!")
		res <- list(model="Single-class classifier", class=rep(classes,nrow(xtst)), prob=rep(1,nrow(xtst)))
		return(res)
	} else if (length(classes) > 2){
		stop("Multi-class setting...")
	}
	
	if(do.normalize){
		d <- normalizeData(xtr, xtst, ...)
		xtr <- d$xtr
		xtst <- d$xtst
	}
	
	model <- glm(y ~ ., data=data.frame(y=as.factor(ytr),as.data.frame(xtr)), family="binomial", ...)
	if(do.stepAIC){
		message("Choosing a model by stepwise AIC...")
		model <- stepAIC(model, trace=0)
	}
	
	pred_prob <- predict(model, as.data.frame(xtst), type="link")
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorLogitLasso <- function(xtr, xtst, ytr, cutoff = 0, do.normalize = TRUE, ...){
	# predictorLogit() with Lasso
	
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



predictorSVM <- function(xtr, xtst, ytr, kernel = "linear", cost = 1, cutoff = 0, do.normalize = TRUE, ...){
	# cost is defaulted constant value without tuning by cross-validation
	# NOTE cross validation is set to minimize misclassification error by default of tune.svm arguments
	
	library(e1071)
  
	classes <- sort(unique(as.character(ytr)))
	
	o <- order(ytr, decreasing = TRUE) # otherwise decision values might be inversed
	model <- tune.svm(xtr[o,], as.factor(ytr[o]), kernel=kernel, cost=cost, scale=do.normalize, type="C-classification", ...)
	
	pred <- predict(model$best.model, xtst, decision.values=T)
	pred_prob <- attr(pred,"decision.values")[,1]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model$best.model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorKNN <- function(xtr, xtst, ytr, k = 5, do.normalize = TRUE, cutoff = 0.5, ...){
	
	library(class)
  
	classes <- sort(unique(as.character(ytr)))
	
	if(do.normalize){
		d <- normalizeData(xtr, xtst, ...)
		xtr <- d$xtr
		xtst <- d$xtst
	}
	
	pred <- knn(xtr, xtst, ytr, k, prob = TRUE, ...)
	
	pred_prob <- attr(pred,"prob"); pred_prob[pred == classes[1]] <- 1 - pred_prob[pred == classes[1]]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=paste0(k,"-nearest neighbour classification"), class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorGBM <- function(xtr, xtst, ytr, n.trees = 100, shrinkage = 0.001, interaction.depth = 3, cutoff = 0, ...){
	
	library(gbm)
  
	classes <- sort(unique(as.character(ytr)))
	ytr <- as.numeric(as.factor(ytr))-1
	
	model <- gbm.fit(xtr, ytr, distribution="bernoulli", interaction.depth=interaction.depth, n.trees=n.trees, shrinkage=shrinkage, verbose=TRUE, ...)
	
	pred_prob <- predict(model, xtst, model$n.trees, type="link")
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorRF <- function(xtr, xtst, ytr, ntrees = 1000, cutoff = 0.5, ...){
	
	library(randomForest)
  
	classes <- sort(unique(as.character(ytr)))
	
	model <- randomForest(xtr, as.factor(ytr), ntree=ntrees, importance=TRUE, proximity=TRUE, do.trace=FALSE, ...)
	
	pred_prob <- predict(model, xtst, type="prob")[,classes[2]]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorDT <- function(xtr, xtst, ytr, cutoff = 0.5, do.prune = FALSE, ...){
	
	library(rpart)
  
	classes <- sort(unique(as.character(ytr)))
	
	model <- rpart(y ~ ., data=data.frame(y=ytr,as.data.frame(xtr)), method="class", ...)
	if(do.prune){
		message("Pruning trees...")
		model <- prune(model, cp=model$cptable[which.min(model$cptable[,"xerror"]),"CP"])
	}
	
	pred_prob <- predict(model, as.data.frame(xtst), type="prob")[,classes[2]]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorNNet <- function(xtr, xtst, ytr, size = 2, cutoff = 0.5, do.normalize = TRUE, ...){
	
	library(nnet)
  
	classes <- sort(unique(as.character(ytr)))
	targets <- class.ind(as.factor(ytr))
	
	if(do.normalize){
		d <- normalizeData(xtr, xtst, ...)
		xtr <- d$xtr
		xtst <- d$xtst
	}
	
	model <- nnet(xtr, targets, size=size, ...)
	
	pred_prob <- predict(model, xtst, type="raw")[,classes[2]]
	pred_class <- rep(classes[1],nrow(xtst)); pred_class[pred_prob > cutoff] <- classes[2]
	
	res <- list(model=model, class=pred_class, prob=pred_prob, cutoff=cutoff)
	return(res)
}



predictorSparseSVM <- function(xtr, xtst, ytr, cost = 1, cutoff = 0, do.normalize = TRUE, ...){
	# L1-regularized L2-loss support vector classification
	
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



########################
###
### others
###
########################

# NA proportion function
NArate <- function(v){mean(is.na(v))}



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


