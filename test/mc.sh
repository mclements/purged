#!/bin/bash
#PBS -o MC
#PBS -j oe
module load Apps/R/3.0.2
cd $PBS_O_WORKDIR
R --slave -f mc.R
