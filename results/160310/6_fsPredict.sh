#!/bin/bash

## split job runs, check *.core.sh for details

# submit first job named "xxx"
qsub -N CVfsPredict -t 401-2400 -o /dev/null 6_fsPredict.core.sh

# now using -hold_jid second job doesnt start until the last job named "xxx" finishes
qsub -hold_jid CVfsPredict -N fsPredict -t 1-400 -o /dev/null 6_fsPredict.core.sh

# now gather results with 7_fsResults.Rmd
qsub -hold_jid fsPredict 7_fsResults.sh
