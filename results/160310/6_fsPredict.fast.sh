#!/bin/bash

## split job runs, check *.core.sh for details

rm logs/*fs*
rm Robj/*fs*

# submit first job named "xxx"
qsub -N CVfsPredict -t 401-409 -o logs 6_fsPredict.core.sh

# now using -hold_jid second job doesnt start until the last job named "xxx" finishes
qsub -hold_jid CVfsPredict -N fsPredict -t 1-9 -o logs 6_fsPredict.core.sh
