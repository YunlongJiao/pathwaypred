#!/bin/bash

## run on crom01 with R version=3.2.1

## -- our name ---
#$ -N pathwayPred
#$ -pe smp 1
#$ -l h_vmem=8G
## Make sure that the .e and .o file arrive in the
## working directory
#$ -cwd
## Merge the standard out and standard error to one file in one folder
#$ -j y
#$ -o logs
#$ -e logs

## set up distributed jobs for nclust range
#$ -t 1-47250
## limit the number of simultaneous jobs
#$ -tc 100

source $HOME/.bashrc

FILE=$PWD/2_param.txt
FIELD1=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f2)
FIELD2=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f3)
FIELD3=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f4)
FIELD4=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f5)
FIELD5=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f6)
FIELD6=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f7)
FIELD7=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f8)
FIELD8=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f9)
FIELD9=$(grep "^$SGE_TASK_ID " $FILE | cut -d' ' -f10)
Rscript $PWD/2_runPredict.R $FIELD1 $FIELD2 $FIELD3 $FIELD4 $FIELD5 $FIELD6 $FIELD7 $FIELD8 $FIELD9
