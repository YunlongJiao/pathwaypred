#!/bin/bash

## run on crom01 with R version=3.2.1

## -- our name ---
#$ -N pathwayPred
#$ -pe smp 1
#$ -l h_vmem=2G
## Make sure that the .e and .o file arrive in the
## working directory
#$ -cwd
## Merge the standard out and standard error to one file
#$ -o logs
#$ -e logs

## set up distributed jobs for nclust range
#$ -t 1-9500
## limit the number of simultaneous jobs
#$ -tc 80

source $HOME/.bashrc

FILE=$PWD/cluster_param.txt
FIELD1=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f2)
FIELD2=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f3)
FIELD3=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f4)
FIELD4=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f5)
FIELD5=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f6)
FIELD6=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f7)
Rscript $PWD/2_runPredict.R $FIELD1 $FIELD2 $FIELD3 $FIELD4 $FIELD5 $FIELD6
