#!/bin/bash

## run on crom01 with R version=3.2.1

## -- our name ---
#$ -N selectPredict
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
#$ -t 1-9450
## limit the number of simultaneous jobs
#$ -tc 100

source $HOME/.bashrc

SEEDFILE=$PWD/3_selectPredict.txt
FIELD1=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f1)
FIELD2=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f2)
FIELD3=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f3)
FIELD4=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f4)
FIELD5=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f5)
FIELD6=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f6)
Rscript $PWD/3_selectPredict.R $FIELD1 $FIELD2 $FIELD3 $FIELD4 $FIELD5 $FIELD6
