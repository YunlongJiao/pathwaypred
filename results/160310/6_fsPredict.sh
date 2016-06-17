#!/bin/bash

## split job runs, check *.core.sh for details

qsub -N CVfsPredict -t 401-2400 6_fsPredict.core.sh

# now using -hold_jid second job doesnt start until job named "job_name_1" finishes
qsub -hold_jid CVfsPredict -N fsPredict -t 1-400 6_fsPredict.core.sh
