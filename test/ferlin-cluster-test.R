
## test
library(Rmpi)
library(snow)
library(parallel)
## set up nodes and cores
NumberOfNodes <- 2
CoresPerNode <- 8
mc.cores <- max(1, NumberOfNodes*CoresPerNode-1) # minus one for master
system.time(lapply(1:100, function(i) mean(1:1e8)))
cl <- makeMPIcluster(mc.cores)
system.time(out <- clusterMap(cl, function(i) mean(1:1e8), 1:100))
stopCluster(cl)
mpi.quit()
