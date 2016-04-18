#!/bin/bash

## run on crom01 with R version=3.2.1

## -- our name ---
#$ -N runPredict
#$ -pe smp 1
#$ -l h_vmem=8G
## Make sure that the .e and .o file arrive in the
## working directory
#$ -cwd
## Merge the standard out and standard error to one file in one folder
#$ -j y
## Disable log output
#$ -o /dev/null
#$ -e /dev/null

## set up distributed jobs for nclust range
#$ -t 1-56700
## limit the number of simultaneous jobs
#$ -tc 100

source $HOME/.bashrc

SEEDFILE=$PWD/2_runPredict.txt
FIELD1=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f1)
FIELD2=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f2)
FIELD3=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f3)
FIELD4=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f4)
FIELD5=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f5)
FIELD6=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f6)
FIELD7=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f7)
FIELD8=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f8)
FIELD9=$(cat $SEEDFILE | head -n $SGE_TASK_ID | tail -n 1 | cut -d' ' -f9)
Rscript $PWD/2_runPredict.R $FIELD1 $FIELD2 $FIELD3 $FIELD4 $FIELD5 $FIELD6 $FIELD7 $FIELD8 $FIELD9
