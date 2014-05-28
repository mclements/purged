#!/bin/bash
qsub -q fas_high -l nodes=2:ppn=8,walltime=0:05:00 cluster.sh
