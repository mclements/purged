#!/bin/bash
module add gcc/4.6.0 
module add R/3.0.2
module add gsl/1.16
R --slave -f $1
