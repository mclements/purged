library(Rmpi)
library(snow)
library(parallel)
require(purged)
##load(file="~/Documents/clients/ted/smokingList-20140528.RData")
load(file="/home/fas/holford/mc2495/src/R/purged/test/smokingList-20140528.RData")
strata <- transform(expand.grid(from=seq(1890,1990,by=5),sex=1:2),
                    to=from+4)
smokingSex <- sapply(smokingList,function(obj) obj$sex)
smokingCohort <- sapply(smokingList,function(obj) obj$cohort)
stratifiedData <- lapply(1:nrow(strata), function(i) {
    smokingList[strata$sex[i]==smokingSex &
                smokingCohort>=strata$from[i] & smokingCohort<=strata$to[i]]
})
cl <- makeMPIcluster(max(1,mpi.universe.size() - 1))
cat(sprintf("Running with %d workers\n", length(cl)))
clusterCall(cl, function() { library(purged); NULL })
do.all <- function(strata) {
    optimObjective <- function(stratum) {
        objective <- function(beta) {
            int01 <- beta[i <- 1]
            int12 <- beta[i <- i+1]
            beta01 <- beta[(i+1):(i <- i+2)]
            beta12 <- beta[(i+1):(i <- i+2)]
            beta20 <- beta[i <- i+1]
            do.call("sum",
                    lapply(stratum, function(obj) {
                        smoking <- obj$smoking
                        mort <- obj$mort
                        .Call("gsl_main2Reclassified",
                              smoking$smkstat - 1, 
                              as.integer(smoking$recall),
                              smoking$agestart,
                              smoking$agequit,
                              smoking$age, # age observed
                              smoking$freq, # frequency (as double)
                              int01,
                              int12,
                              c(10,20,30), # knots01
                              c(20,40,60), # knots12
                              beta01,
                              beta12,
                              beta20,
                              mort$age,
                              mort$Never,
                              mort$Current,
                              mort$Former,
                              FALSE, # debug
                              package="purged")
                    }))
        }
        init <- c(-1.93636374725587, -3.4001339744127, 3.01006761619528, -3.02733667354355, 
                  1.01152056888029, 0.394911534677803, -3.98158254291634)
        try(optim(init,objective,control=list(maxit=500),hessian=TRUE))
    }
    clusterApply(cl,strata, optimObjective)
}
out <- do.all(stratifiedData)
save(out,file="/home/fas/holford/mc2495/src/R/purged/test/out-20140528.RData")
stopCluster(cl)
mpi.quit()
