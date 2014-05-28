#!/bin/bash
#PBS -o CLUSTER
#PBS -j oe
module load Apps/R/3.0.2
module load Rpkgs/DOSNOW
module load Rpkgs/RMPI
cd $PBS_O_WORKDIR
mpirun -n 1 R --slave -f cluster.R
