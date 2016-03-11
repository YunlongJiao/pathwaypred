# run this script on crom01 for prediction results

# read nclust from bash command line
ags <- commandArgs(trailingOnly = TRUE)
xname <- ags[1]
yname <- ags[2]
predictor <- ags[3]

# start! ------------------------------------------------------------------

message("loading features and groups and predictors ...")
load('dat.RData')
source("../../src/func.R") # complementary functions for kernel kmeans for top-k rankings

xtr <- get(xname)
ytr <- get(yname)
objname <- paste('res', xname, yname, predictor, sep = '_')

message('train and predict ... ')
assign(objname, 
       crossValidation(xtr = xtr, ytr = ytr, predictor = predictor, seed = 94151402))

message('plot ROC ... ')
if (!dir.exists('figures')) dir.create('figures')
plotROCcv(res = get(objname), savepath = paste0('figures/', objname, '.pdf'))

message('write results ... ')
score <- c("acc","fpr","tpr","ppv","fval","concordance.index","auroc")
dd <- data.frame(xname, yname, predictor, score, unlist(get(objname)[score]))
if (!file.exists('scores.txt')) file.create('scores.txt')
write.table(dd, file = "scores.txt", row.names = FALSE, col.names = FALSE, append = TRUE)

message('save up !!')
if (!dir.exists('Robj')) dir.create('Robj')
save(list = objname, file = paste0('Robj/', objname, '.RData'))
