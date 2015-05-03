#!/bin/bash

## add dependencies
module add gcc/4.6.0 
module add R/3.0.2
module add gsl/1.16

processes_per_node=8
total_processes=`expr $processes_per_node \* $SP_NODES`

## initiating one processor with mpirun and the others within R
mpirun -np 1 --hostfile $SP_HOSTFILE R --slave -f $1
