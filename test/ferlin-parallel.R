
## test
library(parallel)
library(purged)
## data
load(file="~/src/R/purged/test/smokingList-20140528.RData")
load("~/src/R/purged/test/out-20140528.RData") # previous run
strata <- transform(expand.grid(from=seq(1890,1990,by=5),sex=1:2),
                                        to=from+4)
smokingSex <- sapply(smokingList,function(obj) obj$sex)
smokingCohort <- sapply(smokingList,function(obj) obj$cohort)
stratifiedData <- lapply(1:nrow(strata), function(i) {
      smokingList[strata$sex[i]==smokingSex &
                                  smokingCohort>=strata$from[i] & smokingCohort<=strata$to[i]]
    })
inits <- lapply(out, function(obj) obj$par)
optimObjective <- function(stratum,init,hessian=FALSE) {
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
  ## init <- c(-1.91166383109674, -4.86139268487933, 3.84976690180782, 0.445314035540679,
  ## -0.344483039469035, 0.979805469221358, -3.53232595501467)
  ## out <- nlminb(init,objective,control=list(iter.max=5))
  out <- nlminb(init,objective,control=list(trace=1))
  out$objective <- objective
  out$coefficients <- out$par
  if (hessian) {
    out$hessian <- optimHess(coef(out), objective)
  }
  return(out)
}
system.time(out <- mclapply(1:42, function(i) {cat("Object",i,"\n"); try(optimObjective(stratifiedData[[i]], inits[[i]]))}, mc.cores=7))
save(out,file="~/src/R/purged/test/out-20140528_D.RData")
