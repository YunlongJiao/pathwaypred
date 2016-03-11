# run this script on crom01 for prediction results

# read nclust from bash command line
ags <- commandArgs(trailingOnly = TRUE)
xname <- ags[1]
yname <- ags[2]
predictor <- ags[3]

# start! ------------------------------------------------------------------

message("loading features and groups and predictors ...")
load('Robj/dat.RData')
source("../../src/func.R") # complementary functions for kernel kmeans for top-k rankings

xtr <- get(xname)
ytr <- get(yname)
savepath <- paste('figures/ROC', xname, yname, predictor, '.pdf', sep = '_')
objname <- paste('res', xname, yname, predictor, sep = '_')
if (!dir.exists('figures')) dir.create('figures')

message('train and predict ... ')
assign(objname, 
       crossValidation(xtr = xtr, ytr = ytr, predictor = predictor, savepath = savepath, seed = 94151402))

message('write results ... ')
score <- c("acc","fpr","tpr","ppv","fval","concordance.index","auroc")
dd <- data.frame(xname, yname, predictor, score, unlist(get(objname)[score]))
if (!file.exists('2_scores.txt')) file.create('2_scores.txt')
write.table(dd, file = "2_scores.txt", row.names = FALSE, col.names = FALSE, append = TRUE)

message('save up !!')
if (!dir.exists('Robj/2_runPredict')) dir.create('Robj/2_runPredict')
save(list = objname, file = paste0('Robj/2_runPredict/', objname, '.RData'))
