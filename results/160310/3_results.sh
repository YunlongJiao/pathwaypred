#!/bin/bash

## run on crom01 with R version=3.2.1

## -- our name ---
#$ -N results
#$ -l h_vmem=32G
## Make sure that the .e and .o file arrive in the
## working directory
#$ -cwd
## Merge the standard out and standard error to one file in one folder
#$ -j y
#$ -o logs
#$ -e logs
## Long time execution
#$ -l long=true
#$ -q long.q

source $HOME/.bashrc
Rscript -e "knitr::knit('$PWD/3_results.Rmd')"
