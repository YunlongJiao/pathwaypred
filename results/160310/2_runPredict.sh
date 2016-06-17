#!/bin/bash

## split job runs, check *.core.sh for details

qsub -t 1-50000 2_runPredict.core.sh
qsub -t 50001-79200 2_runPredict.core.sh
