#!/bin/bash
qsub -q fas_high -l nodes=6:ppn=8,walltime=01:00:00 cluster.sh
