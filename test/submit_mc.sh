#!/bin/bash
qsub -q fas_high -l nodes=1:ppn=8,walltime=0:05:00 mc.sh
