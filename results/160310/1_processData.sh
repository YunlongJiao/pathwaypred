#!/bin/bash

## run on cipf cluster with R version=3.2.1

## -- our name ---
#$ -N processData
## Make sure that the .e and .o file arrive in the
## working directory
#$ -cwd
## Merge the standard out and standard error to one file
#$ -o logs/
#$ -e logs/

source $HOME/.bashrc
Rscript /home/yjiao/pathwaypred/results/160310/1_processData.R
