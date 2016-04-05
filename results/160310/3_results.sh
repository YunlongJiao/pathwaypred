#!/bin/bash

## run on crom01 with R version=3.2.1

## -- our name ---
#$ -N processData
## Make sure that the .e and .o file arrive in the
## working directory
#$ -cwd
## Merge the standard out and standard error to one file
#$ -o logs
#$ -e logs

source $HOME/.bashrc
Rscript $PWD/3_results.R
