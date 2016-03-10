# run on cipf cluster with R version=3.2.1

# -- our name ---
#$ -N processData
#$ -S /bin/sh
# Make sure that the .e and .o file arrive in the
# working directory
#$ -cwd
# Merge the standard out and standard error to one file
#$ -o logs/
#$ -e logs/

source $HOME/.bashrc
R CMD BATCH /home/yjiao/pathwaypred/results/160310/1_processData.R
