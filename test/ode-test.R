refresh
## read in the data
if (doOnce <- FALSE) {
    setwd("~/Documents/clients/ted")
    files <- within(expand.grid(sex=1:2,type=c("Never","Former","Current")),
                    { filename <- sprintf("all_final_results_%s%s.csv",type,sex) })
    mort <- do.call("rbind",
                    lapply(1:nrow(files),function(i) {
                        data <- read.csv(files$filename[i])
                        data$sex <- files$sex[i]
                        data$smkstat <- switch(as.character(files$type[i]),
                                               Never=1,Current=2,Former=3)
                        data
                    }))
    mort <- transform(transform(mort, cohort=year),year=cohort+age)
    require(foreign)
    smoking <- within(read.dta("cisnet_apc_mar19_agg.dta"),
                      { cohort <- year-age })
    Lab.palette <- colorRampPalette(c("white", "darkblue"), space = "Lab")
    x <- subset(smoking,select=c(year,age))
    jpeg("~/src/R/purged/test/nhis_density.jpg",width=480*4,height=480*4,pointsize=12*4)
    smoothScatter(x, colramp = Lab.palette,xlab="Calendar year",ylab="Age (years)")
    dev.off()
    pdf("~/src/R/purged/test/nhis_density.pdf")
    smoothScatter(x, colramp = Lab.palette,xlab="Calendar year",ylab="Age (years)")
    dev.off()
    ## recall: 0=current status, 1=recall, 2=age quit only
    ##
    cohorts <- unique(subset(smoking,cohort>=1890,select=c(sex,cohort)))
    smokingList <-
        lapply(1:nrow(cohorts), function(i) {
            .sex <- cohorts$sex[i]
            .cohort <- cohorts$cohort[i]
            list(sex=.sex,
                 cohort=.cohort,
                 smoking=subset(smoking,cohort==.cohort & sex==.sex),
                 mort=data.frame(age=subset(mort,cohort==.cohort & sex==.sex & smkstat == 1)$age,
                 Never=subset(mort,cohort==.cohort & sex==.sex & smkstat == 1)$final,
                 Current=subset(mort,cohort==.cohort & sex==.sex & smkstat == 2)$final,
                 Former=subset(mort,cohort==.cohort & sex==.sex & smkstat == 3)$final))
        })
    cohorts <- unique(subset(mort,cohort>=1890,select=c(sex,cohort)))
    mortList <-
        lapply(1:nrow(cohorts), function(i) {
            .sex <- cohorts$sex[i]
            .cohort <- cohorts$cohort[i]
            list(sex=.sex,cohort=.cohort,
                 data.frame(age=subset(mort,cohort==.cohort & sex==.sex & smkstat == 1)$age,
                 Never=subset(mort,cohort==.cohort & sex==.sex & smkstat == 1)$final,
                 Current=subset(mort,cohort==.cohort & sex==.sex & smkstat == 2)$final,
                 Former=subset(mort,cohort==.cohort & sex==.sex & smkstat == 3)$final))
        })
    ## cohorts <- data.frame(t(sapply(smokingList,function(obj) c(sex=obj$sex,cohort=obj$cohort))))
    ## smokingList <- smokingList[with(cohorts,order(sex,cohort))]
    class(smokingList) <- c("lookupList","list")
    subset.lookupList <- function(obj,subset) {
        e <- substitute(subset)
        obj[sapply(obj, function(x) eval(e,x,parent.frame()))]
    }
    save(smokingList,mortList,subset.lookupList,file="~/Documents/clients/ted/smokingList-20140528.RData")
}
require(parallel)
require(purged)
require(survival)
load(file="~/src/R/purged/test/smokingList-20140528.RData")
strata <- transform(expand.grid(from=seq(1890,1990,by=5),sex=1:2),
                    to=from+4)
smokingSex <- sapply(smokingList,function(obj) obj$sex)
smokingCohort <- sapply(smokingList,function(obj) obj$cohort)
stratifiedData <- lapply(1:nrow(strata), function(i) {
    smokingList[strata$sex[i]==smokingSex &
                smokingCohort>=strata$from[i] & smokingCohort<=strata$to[i]]
})
## simplify the data: combine the 5 single-year cohorts
stratifiedData2 <- lapply(stratifiedData, function(obj)
                 list(sex=obj[[1]]$sex, cohort=obj[[3]]$cohort, mort=obj[[3]]$mort,
                      smoking=do.call("rbind", lapply(obj, function(elt) elt$smoking))))

## testing and optimisation for ps
test <- function(callName="call_purged_ps",init,stratum,sp01=0.1,sp12=1,output_type="negll",debug=FALSE) {
    nterm01 <- nterm12 <- 5
    smoking <- stratum$smoking
    mort <- stratum$mort
    args <- list(finalState=as.integer(smoking$smkstat - 1), 
                 recall=as.integer(smoking$recall),
                 time1=smoking$agestart,
                 time2=smoking$agequit,
                 time3=smoking$age, # age observed
                 freq=smoking$freq, # frequency (as double)
                 lower01=10,
                 upper01=40,
                 nterm01=as.integer(nterm01),
                 pmatrix01=attr(pspline(c(10,40),nterm=nterm01),"pparm"),
                 sp01=sp01,
                 lower12=10,
                 upper12=70,
                 nterm12=as.integer(nterm12),
                 pmatrix12=attr(pspline(c(10,70),nterm=nterm12),"pparm"),
                 sp12=sp12,
                 ages0=mort$age,
                 mu0=mort$Never,
                 mu1=mort$Current,
                 mu2=mort$Former,
                 N=nrow(smoking),
                 nages0=length(mort$age),
                 nbeta01=as.integer(nterm01+2),
                 nbeta12=as.integer(nterm12+2),
                 output_type=output_type,
                 debug=debug)
    parseBeta <- function(beta) {
        args$int01 <- beta[i <- 1]
        args$beta01 <- beta[(i+1):(i <- i+nterm01+2)]
        args$int12 <- beta[i <- i+1]
        args$beta12 <- beta[(i+1):(i <- i+nterm12+2)]
        args$beta20 <- beta[i <- i+1]
        args
    }
    objective <- function(beta) {
        args <- parseBeta(beta)
        args$output_type="pnegll"
        value <- .Call(callName, args, package="purged")
        cat(value); cat("\n")
        value
    }
    gradient <- function(beta) {
        args <- parseBeta(beta)
        args$output_type="pnegll_gradient"
        .Call(callName, args, package="purged")
    }
    fun <- function(beta) {
        args <- parseBeta(beta)
        args$output_type <- output_type
        .Call(callName, args, package="purged")
    }
    if (output_type == "optim") {
        fit <- nlminb(init,objective,gradient)
        fit$hessian <- optimHess(fit$par,objective,gradient)
        fit$vcov <- solve(fit$hessian)
        return(fit)
    }
    return(fun(init))
}
init <- c(-3,
          rep(0.1,5+2),
          -4,
          rep(0.1,5+2),
          log(0.01))
init <- c(-2.2084,
          1.7865, 3.9082, 4.2282, 2.6938, 1.7596, 2.0776, 0.8403, -0.5627,
          -5.2397,
          -0.9299, -1.105, -0.6869, 0.1301, 1.1833, 2.0894, 
          -4.3141)
test("call_purged_ps",init,sp01=0.1,sp12=1,stratum=stratifiedData2[[10]],output_type="negll",debug=FALSE) 

test("call_purged_ps",init,sp01=0.1,sp12=1,stratum=stratifiedData2[[10]],output_type="pnegll_gradient")
test("call_purged_ps",init,sp01=0.1,sp12=1,stratum=stratifiedData2[[10]],output_type="negll",debug=FALSE) 
test("call_purged_ps",init,sp01=0.1,sp12=1,stratum=stratifiedData2[[10]],output_type="pnegll",debug=FALSE) 
## temp <- stratifiedData[[10]][[1]]
## temp$smoking[1+c(1253:1259),]
##debug(test)
dbeta <- function(x,i,scale,eps=1e-4) { x[i] <- x[i]+scale*eps; x }
dtest <- function(...,beta=init,i=1,eps=1e-4)
    (test(...,init=dbeta(beta,i,scale=1,eps=eps)) -
     test(...,init=dbeta(beta,i,scale=-1,eps=eps))) / (2*eps)
test("call_purged_ps",init,sp01=0.1,sp12=1,stratum=stratifiedData2[[10]],output_type="pnegll_gradient")
zapsmall(sapply(1:length(init), function(i)
                dtest("call_purged_ps",beta=init,stratum=stratifiedData2[[10]],sp01=0.1,sp12=1,
                      output_type="pnegll",i=i, eps=1e-6)))
test("call_purged_ps",init,sp01=0.1,sp12=1,stratum=stratifiedData2[[10]],output_type="negll_gradient")
zapsmall(sapply(1:length(init), function(i)
                dtest("call_purged_ps",beta=init,stratum=stratifiedData2[[10]],sp01=0.1,sp12=1,
                      output_type="negll",i=i, eps=1e-6)))
system.time(fit <- test("call_purged_ps",init,sp01=0.1,sp12=1,stratum=stratifiedData2[[10]],
                        output_type="optim",debug=FALSE))


##
## compare finite-differences with calculated gradients
dbeta <- function(x,i,scale,eps=1e-4) { x[i] <- x[i]+scale*eps; x }
dtest <- function(...,beta=init,i=1,eps=1e-4)
    (test(...,init=dbeta(beta,i,scale=1,eps=eps)) -
     test(...,init=dbeta(beta,i,scale=-1,eps=eps))) / (2*eps)
## sapply(1:7, function(i)
##        dtest("call_purged_ns",stratum=stratifiedData[[10]],i=i))
## test("call_purged_ns",init=init,stratum=stratifiedData[[10]],output_type="negll_gradient")
## FAIL: estimates are different
test("call_purged_ns",init=init,stratum=stratifiedData[[10]],output_type="dP_ij") /
    zapsmall(sapply(1:7, function(i)
                    dtest("call_purged_ns",beta=init,stratum=stratifiedData[[10]],
                          output_type="P_ij",i=i)))
test("call_purged_ns",init=init,stratum=stratifiedData[[10]],output_type="negll_gradient") /
    zapsmall(sapply(1:7, function(i)
                    dtest("call_purged_ns",beta=init,stratum=stratifiedData[[10]],
                          output_type="negll",i=i)))



## simple testing for ns (OK on valgrind)
test <- function(callName="call_purged_ns",init,stratum,output_type="negll",debug=FALSE) {
    obj <- stratum #[[1]]
    nterm01 <- 2
    nterm12 <- 2
    P <- function(beta) {
        int01 <- beta[i <- 1]
        beta01 <- beta[(i+1):(i <- i+nterm01)]
        int12 <- beta[i <- i+1]
        beta12 <- beta[(i+1):(i <- i+nterm12)]
        beta20 <- beta[i <- i+1]
        smoking <- obj$smoking[1:5,]
        mort <- obj$mort
        args <- list(finalState=as.integer(smoking$smkstat - 1), 
                     recall=as.integer(smoking$recall),
                     time1=smoking$agestart,
                     time2=smoking$agequit,
                     time3=smoking$age, # age observed
                     freq=smoking$freq, # frequency (as double)
                     int01=int01,
                     int12=int12,
                     beta01=beta01,
                     beta12=beta12,
                     beta20=beta20,
                     ages0=mort$age,
                     mu0=mort$Never,
                     mu1=mort$Current,
                     mu2=mort$Former,
                     N=nrow(smoking),
                     knots01=c(10,20,30),
                     knots12=c(20,40,60),
                     nages0=length(mort$age),
                     nbeta01=length(beta01),
                     nbeta12=length(beta12),
                     output_type=output_type,
                     debug=debug)
        .Call(callName, args, package="purged")
    }
    return(P(init))
}
init <- c(-1.91166383109674, 3.84976690180782, 0.445314035540679,
          -4.86139268487933, -0.344483039469035, 0.979805469221358,
          -3.53232595501467)
##debug(test)
test("call_purged_ns",init=init,stratum=stratifiedData2[[10]], output_type="negll_gradient")
test("call_purged_ns",init=init,stratum=stratifiedData2[[10]], output_type="negll")


## optimize ns
require(bbmle)
optimObjective <- function(obj,init,debug=FALSE) {
    smoking <- obj$smoking
    mort <- obj$mort
    args <- list(finalState=as.integer(smoking$smkstat - 1), 
                 recall=as.integer(smoking$recall),
                 time1=smoking$agestart,
                 time2=smoking$agequit,
                 time3=smoking$age, # age observed
                 freq=smoking$freq, # frequency (as double)
                 ages0=mort$age,
                 mu0=mort$Never,
                 mu1=mort$Current,
                 mu2=mort$Former,
                 N=nrow(smoking),
                 knots01=c(10,20,30),
                 knots12=c(20,40,60),
                 nages0=length(mort$age),
                 nbeta01=2L,
                 nbeta12=2L,
                 debug=FALSE)
    parseBeta <- function(beta) {
        args$int01 <- beta[i <- 1]
        args$beta01 <- beta[(i+1):(i <- i+2)]
        args$int12 <- beta[i <- i+1]
        args$beta12 <- beta[(i+1):(i <- i+2)]
        args$beta20 <- beta[i <- i+1]
        args
    }
    objective <- function(beta) {
        args <- parseBeta(beta)
        args$output_type <- "negll"
        value <- .Call("call_purged_ns",args,package="purged")
        if (debug) cat(sprintf("value=%f\n",value))
        value
    }
    gradient <- function(beta) {
        args <- parseBeta(beta)
        args$output_type <- "negll_gradient"
        value <- .Call("call_purged_ns",args,package="purged")
        if (debug) { cat("gradient="); cat(value); cat("\n") }
        value
    }
    ##optim(init,objective,control=list(trace=2)) # slow
    ## optim(init,objective,gradient,method="BFGS") # FAILS
    fit <- nlminb(init,objective,gradient)
    fit$hessian <- optimHess(fit$par,objective,gradient)
    fit$vcov <- solve(fit$hessian)
    fit
}
init <- c(-1.91166383109674,
          3.84976690180782, 0.445314035540679, 
          -4.86139268487933,
          -0.344483039469035, 0.979805469221358,
          -3.53232595501467)
system.time(optim1 <- optimObjective(stratifiedData2[[1]], init*1.1))

## simple testing for ns
test <- function(callName="call_purged_ns",init,stratum,output_type="negll",debug=FALSE) {
    obj <- stratum #[[1]]
    nterm01 <- 2
    nterm12 <- 2
    P <- function(beta) {
        int01 <- beta[i <- 1]
        beta01 <- beta[(i+1):(i <- i+nterm01)]
        int12 <- beta[i <- i+1]
        beta12 <- beta[(i+1):(i <- i+nterm12)]
        beta20 <- beta[i <- i+1]
        smoking <- obj$smoking[1:5,]
        mort <- obj$mort
        args <- list(finalState=as.integer(smoking$smkstat - 1), 
                     recall=as.integer(smoking$recall),
                     time1=smoking$agestart,
                     time2=smoking$agequit,
                     time3=smoking$age, # age observed
                     freq=smoking$freq, # frequency (as double)
                     int01=int01,
                     int12=int12,
                     beta01=beta01,
                     beta12=beta12,
                     beta20=beta20,
                     ages0=mort$age,
                     mu0=mort$Never,
                     mu1=mort$Current,
                     mu2=mort$Former,
                     N=nrow(smoking),
                     knots01=c(10,20,30),
                     knots12=c(20,40,60),
                     nages0=length(mort$age),
                     nbeta01=length(beta01),
                     nbeta12=length(beta12),
                     output_type=output_type,
                     debug=debug)
        .Call(callName, args, package="purged")
    }
    return(P(init))
}
init <- c(-1.91166383109674, 3.84976690180782, 0.445314035540679,
          -4.86139268487933, -0.344483039469035, 0.979805469221358,
          -3.53232595501467)
##debug(test)
test("call_purged_ns",init=init,stratum=stratifiedData2[[10]])
##
## compare finite-differences with calculated gradients
dbeta <- function(x,i,scale,eps=1e-4) { x[i] <- x[i]+scale*eps; x }
dtest <- function(...,beta=init,i=1,eps=1e-4)
    (test(...,init=dbeta(beta,i,scale=1,eps=eps)) -
     test(...,init=dbeta(beta,i,scale=-1,eps=eps))) / (2*eps)
## sapply(1:7, function(i)
##        dtest("call_purged_ns",stratum=stratifiedData[[10]],i=i))
## test("call_purged_ns",init=init,stratum=stratifiedData[[10]],output_type="negll_gradient")
## FAIL: estimates are different
test("call_purged_ns",init=init,stratum=stratifiedData[[10]],output_type="dP_ij") /
    zapsmall(sapply(1:7, function(i)
                    dtest("call_purged_ns",beta=init,stratum=stratifiedData[[10]],
                          output_type="P_ij",i=i)))
test("call_purged_ns",init=init,stratum=stratifiedData[[10]],output_type="negll_gradient") /
    zapsmall(sapply(1:7, function(i)
                    dtest("call_purged_ns",beta=init,stratum=stratifiedData[[10]],
                          output_type="negll",i=i)))


## simple testing for ns
test <- function(callName="call_purged_ns",init,stratum,output_type="negll",debug=FALSE) {
    obj <- stratum #[[1]]
    nterm01 <- 2
    nterm12 <- 2
    P <- function(beta) {
        int01 <- beta[i <- 1]
        beta01 <- beta[(i+1):(i <- i+nterm01)]
        int12 <- beta[i <- i+1]
        beta12 <- beta[(i+1):(i <- i+nterm12)]
        beta20 <- beta[i <- i+1]
        ## do.call("sum",
        ##         mclapply(stratum, function(obj) {
                    smoking <- obj$smoking[1:5,]
                    mort <- obj$mort
                    .Call(callName,
                          list(finalState=as.integer(smoking$smkstat - 1), 
                               recall=as.integer(smoking$recall),
                               time1=smoking$agestart,
                               time2=smoking$agequit,
                               time3=smoking$age, # age observed
                               freq=smoking$freq, # frequency (as double)
                               int01=int01,
                               int12=int12,
                               beta01=beta01,
                               beta12=beta12,
                               beta20=beta20,
                               ages0=mort$age,
                               mu0=mort$Never,
                               mu1=mort$Current,
                               mu2=mort$Former,
                               N=nrow(smoking),
                               knots01=c(10,20,30),
                               knots12=c(20,40,60),
                               nages0=length(mort$age),
                               nbeta01=length(beta01),
                               nbeta12=length(beta12),
                               output_type=output_type,
                               debug=debug), 
                          package="purged")
                ## }, mc.cores=2))
    }
    ## return(pnegll(init))
    return(P(init))
}
init <- c(-1.91166383109674, 3.84976690180782, 0.445314035540679,
          -4.86139268487933, -0.344483039469035, 0.979805469221358,
          -3.53232595501467)
##debug(test)
## system.time(temp3 <- test("call_purged_ns",init=init,stratum=stratifiedData[[10]]))
## compare finite-differences with calculated gradients
dbeta <- function(x,i,scale,eps=1e-4) { x[i] <- x[i]+scale*eps; x }
dtest <- function(...,beta=init,i=1,eps=1e-4)
    (test(...,init=dbeta(beta,i,scale=1,eps=eps)) -
     test(...,init=dbeta(beta,i,scale=-1,eps=eps))) / (2*eps)
## sapply(1:7, function(i)
##        dtest("call_purged_ns",stratum=stratifiedData2[[10]],i=i))
## test("call_purged_ns",init=init,stratum=stratifiedData2[[10]],output_type="negll_gradient")
## FAIL: estimates are different
## test("call_purged_ns",init=init,stratum=stratifiedData2[[10]],output_type="dP_ij") /
## zapsmall(sapply(1:7, function(i)
##                 dtest("call_purged_ns",beta=init,stratum=stratifiedData2[[10]],
##                       output_type="P_ij",i=i)))
##
test("call_purged_ns",init=init,stratum=stratifiedData2[[10]],output_type="negll")
test("call_purged_ns",init=init,stratum=stratifiedData2[[10]],output_type="negll_gradient")
zapsmall(sapply(1:7, function(i)
                dtest("call_purged_ns",beta=init,stratum=stratifiedData2[[10]],
                      output_type="negll",i=i)))


##system.time(out <- mclapply(1:42, function(i) {cat("Object",i,"\n"); try(optimObjective(stratifiedData[[i]], inits[[i]]))}, mc.cores=2))
##save(out,file="~/src/R/purged/test/out-20150502-full.RData")


optimObjective <- function(stratum,sp01=0.1,sp12=1,hessian=FALSE) {
    nterm01 <- nterm12 <- 5
    args <- 
    pnegll <- function(beta) {
        int01 <- beta[i <- 1]
        int12 <- beta[i <- i+1]
        beta01 <- beta[(i+1):(i <- i+nterm01+2)]
        beta12 <- beta[(i+1):(i <- i+nterm12+2)]
        beta20 <- beta[i <- i+1]
        ## print(beta,digits=3)
        do.call("sum",
                mclapply(stratum, function(obj) {
                    smoking <- obj$smoking
                    mort <- obj$mort
                    .Call("gsl_main2ReclassifiedPS",
                          smoking$smkstat - 1, 
                          as.integer(smoking$recall),
                          smoking$agestart,
                          smoking$agequit,
                          smoking$age, # age observed
                          smoking$freq, # frequency (as double)
                          int01,
                          int12,
                          10, 40, nterm01, attr(pspline(c(10,40),nterm=nterm01),"pparm"), sp01,
                          10, 70, nterm12, attr(pspline(c(10,70),nterm=nterm12),"pparm"), sp12,
                          beta01,
                          beta12,
                          beta20,
                          mort$age,
                          mort$Never,
                          mort$Current,
                          mort$Former,
                          FALSE, # debug
                          package="purged")$pnegll
                }, mc.cores=2))
    }
    negll <- function(beta) {
        int01 <- beta[i <- 1]
        int12 <- beta[i <- i+1]
        beta01 <- beta[(i+1):(i <- i+nterm01+2)]
        beta12 <- beta[(i+1):(i <- i+nterm12+2)]
        beta20 <- beta[i <- i+1]
        do.call("sum",
                mclapply(stratum, function(obj) {
                    smoking <- obj$smoking
                    mort <- obj$mort
                    .Call("gsl_main2ReclassifiedPS",
                          smoking$smkstat - 1, 
                          as.integer(smoking$recall),
                          smoking$agestart,
                          smoking$agequit,
                          smoking$age, # age observed
                          smoking$freq, # frequency (as double)
                          int01,
                          int12,
                          10, 40, nterm01, attr(pspline(c(10,40),nterm=nterm01),"pparm"), sp01,
                          10, 70, nterm12, attr(pspline(c(10,70),nterm=nterm12),"pparm"), sp12,
                          beta01,
                          beta12,
                          beta20,
                          mort$age,
                          mort$Never,
                          mort$Current,
                          mort$Former,
                          FALSE, # debug
                          package="purged")$negll
                }, mc.cores=2))
    }
    ## return(pnegll(init))
    out <- nlminb(init,pnegll,control=list(trace=1,rel.tol=1e-7))
    out$negll <- negll # function 
    out$pnegll <- pnegll # function
    out$coefficients <- out$par
    if (hessian) {
        out$hessian <- optimHess(coef(out), pnegll)
        ## Hl <- numDeriv::hessian(llike,coef)
        Hl <- optimHess(coef(out), negll)
        Hinv <- solve(out$hessian)
        out$trace <- sum(diag(Hinv %*% Hl))
    }
    return(out)
}
init <- c(-3,
          -4,
          rep(0.1,5+2),
          rep(0.1,5+2),
          log(0.01))
## init <- c(-2.2084, -5.2397, 1.7865, 3.9082, 4.2282, 2.6938, 1.7596, 2.0776, 
##           0.8403, -0.5627, -0.9299, -1.105, -0.6869, 0.1301, 1.1833, 2.0894, 
##           -4.3141)
init <- c(-2.2403,
          -5.2734,
          1.7934, 3.9026, 4.1966, 2.6475, 1.7147, 2.032, 0.7888,
          -0.1281, -0.2629, -0.4078, 0.0032, 0.8361, 1.9041, 2.8168, 
          -4.4402)
system.time(fit1890.1 <- optimObjective(stratifiedData[[1]],sp01=0.1,sp12=1,hessian=TRUE))
str(fit1890.1)
with(fit1890.1, pnegll(par))
with(fit1890.1, dput(round(par,4)))
diag(solve(fit1890.1$hessian))
##
## coef.nlminb <- coef(fit1890.1)
coef.nlminb <- c(-2.30426533639024, -5.59189828926138, 1.7591741291095, 3.94699944347449, 
4.05923551226308, 2.39407490385323, 1.39327811089255, 1.92043040616437, 
1.37925846582657, 0.705860430956955, 0.138782801435078, -0.970265049666329, 
0.0159183612281592, 0.0379401144182051, 0.100106028532845, 0.263734898066751, 
0.388138157647989, 0.201325379289236, -0.059788133365412, -0.13877113685924, 
-0.0194342640363815, 0.44876575129484, -4.71103159635782)
with(fit1890.1, pnegll(coef.nlminb))




test <- function(callName="gsl_main2ReclassifiedPS",init,stratum,sp01=0.1,sp12=1) {
    obj <- stratum[[1]]
    nterm01 <- nterm12 <- 5
    P <- function(beta) {
        int01 <- beta[i <- 1]
        int12 <- beta[i <- i+1]
        beta01 <- beta[(i+1):(i <- i+nterm01+2)]
        beta12 <- beta[(i+1):(i <- i+nterm12+2)]
        beta20 <- beta[i <- i+1]
        ## do.call("sum",
        ##         mclapply(stratum, function(obj) {
                    smoking <- obj$smoking[1:5,]
                    mort <- obj$mort
                    .Call(callName,
                          as.integer(smoking$smkstat - 1), 
                          as.integer(smoking$recall),
                          smoking$agestart,
                          smoking$agequit,
                          smoking$age, # age observed
                          smoking$freq, # frequency (as double)
                          int01,
                          int12,
                          10, 40, as.integer(nterm01), attr(pspline(c(10,40),nterm=nterm01),"pparm"), sp01,
                          10, 70, as.integer(nterm12), attr(pspline(c(10,70),nterm=nterm12),"pparm"), sp12,
                          beta01,
                          beta12,
                          beta20,
                          mort$age,
                          mort$Never,
                          mort$Current,
                          mort$Former,
                          TRUE, # debug
                          package="purged")
                ## }, mc.cores=2))
    }
    ## return(pnegll(init))
    return(P(init))
}
init <- c(-3,
          -4,
          rep(0.1,5+2),
          rep(0.1,5+2),
          log(0.01))
init <- c(-2.2084, -5.2397, 1.7865, 3.9082, 4.2282, 2.6938, 1.7596, 2.0776, 
          0.8403, -0.5627, -0.9299, -1.105, -0.6869, 0.1301, 1.1833, 2.0894, 
          -4.3141)
system.time(temp3 <- test("gsl_main2ReclassifiedPS2",init,stratifiedData[[10]],sp01=0.1,sp12=1))
system.time(temp2 <- test("gsl_main2ReclassifiedPS",init,stratifiedData[[10]],sp01=0.1,sp12=1))
temp2$negll
temp3$negll
temp2$pnegll
temp3$pnegll
head(stratifiedData[[10]][[1]]$smoking)


## Reimplement the example in R
require(splines)
require(deSolve)
require(survival)
parms <- list(beta=init)
d <- stratifiedData[[1]][[1]]
mort <- d$mort
mu0 <- splinefun(mort$age, mort$Never)
mu1 <- splinefun(mort$age, mort$Current)
mu2 <- splinefun(mort$age, mort$Former)
predict.pspline <-
function (object, newx, ...) 
{
    if (missing(newx)) 
        return(object)
    a <- c(list(x = newx, penalty = FALSE), attributes(object)[c("degree", "df", 
        "Boundary.knots", "intercept", "nterm", "eps", "method")])
    newobj <- do.call("pspline", a)
    newobj
}
nsfun <- function(knots,centre=knots[1]) {
    nsobj <- ns(knots,df=length(knots)-1)
    nsref <- predict(nsobj,centre)
    fun <- function(t,beta) apply(predict(nsobj,t),1,
                                  function(row) exp(beta[1]+sum((row-nsref) * beta[-1])))
    fun
}
psfun <- function(Boundary.knots,nterm,centre=Boundary.knots[1]) {
    obj <- pspline(Boundary.knots,nterm=nterm)
    ref <- as.vector(predict(obj,centre))
    ## function(newx) predict(obj,c(Boundary.knots[1]+1.2345e-7,newx))[-1,,drop=FALSE]
    fun <- function(newx,beta) apply(predict(obj,newx),1,
                                  function(row) exp(beta[1]+sum((row-ref) * beta[-1])))
    fun
}
##undebug(psfun)
alpha01 <- psfun(c(10,40),5,20)
alpha12 <- psfun(c(10,70),5,40)
##
with(list(a=1), {
    nterm01 <- nterm12 <- 5
    int01 <- parms$beta[i <- 1]
    int12 <- parms$beta[i <- i+1]
    beta01 <- parms$beta[(i+1):(i <- i+nterm01+2)]
    beta12 <- parms$beta[(i+1):(i <- i+nterm12+2)]
    beta23 <- parms$beta[i <- i+1]
    t <- c(0,seq(.5,71.5,by=1))
    a01 <- alpha01(t,c(int01,beta01))
    a12 <- alpha12(t,c(int12,beta12))
    a23 <- exp(beta23)
    par(mfrow=1:2)
    plot(t,a01,type="l")
    plot(t,a12,type="l")
})
##
bounds <- function(t,lower,upper) pmax(pmin(t,upper),lower)
Never <- 1; Current <- 2; Former <- 3; Reclassified <- 4; Death <- 5
odefunc <- function(t,y,parms) {
    nterm01 <- nterm12 <- 5
    int01 <- parms$beta[i <- 1]
    int12 <- parms$beta[i <- i+1]
    beta01 <- parms$beta[(i+1):(i <- i+nterm01+2)]
    beta12 <- parms$beta[(i+1):(i <- i+nterm12+2)]
    beta23 <- parms$beta[i <- i+1]
    a01 <- alpha01(t,c(int01,beta01))
    a12 <- alpha12(t,c(int12,beta12))
    a23 <- exp(beta23)
    m0 <- mu0(bounds(t,0.5,99.5))
    m1 <- mu1(bounds(t,0.5,99.5))
    m2 <- mu2(bounds(t,0.5,99.5))
    list(c(-y[Never]*a01,
           -y[Current]*(m1-m0+a12)+y[Never]*a01,
           -y[Former]*(m2-m0+a23)+y[Current]*a12,
           y[Former]*a23,y[Current]*(m1-m0)+y[Former]*(m2-m0)),
         c(a01=a01,a12=a12,a23=a23,m0=m0,m1=m1,m2=m2))
}     
Ps <- function(i,t,parms,size=5) {
    y <- rep(0,size)
    y[i] <- 1
    ode(y=y, t=t, parms, func = odefunc, rtol=1e-8, atol=1e-8)[-1,2:6]
}
require(abind)
out <- do.call("abind",c(lapply(1:5,function(i) Ps(i,c(0,seq(.5,71.5,by=1)), parms)),list(along=3)))
out <- aperm(out,3:1)
log(sum(out[1,c(1,4),71])/sum(out[1,1:4,71]))*80
log(out[1,2,71]/sum(out[1,1:4,71]))*81
log(out[1,3,72]/sum(out[1,1:4,72]))*8
##
temp <- temp3$P
dim(temp) <- c(5,5,ncol(temp))
v <- temp[,,26+1]; 68*log((v[1]+v[1+(4-1)*5])/(v[1]+v[1+(2-1)*5]+v[1+(3-1)*5]+v[1+(4-1)*5]))
mat <- solve(temp[,,17+1]) %*% temp[,,26+1]
v <- temp[,,26+1]; 68*log((v[1]+v[1+(4-1)*5])/(v[1]+v[1+(2-1)*5]+v[1+(3-1)*5]+v[1+(4-1)*5]))



temp <- temp3$P
dim(temp) <- c(5,5,ncol(temp))
temp <- zapsmall(temp) # and recalculate
temp <- apply(temp,3,function(m) cbind(m[,-5],1-rowSums(m[,-5])))
dim(temp) <- c(5,5,ncol(temp))
inverses <- apply(temp,3,solve)
dim(temp2$inverse) <- dim(inverses) <- c(5,5,dim(temp)[3])
inverses[,,1]/temp2$inverse[,,1]
inverses[,,10]/temp2$inverse[,,10]
##
fun <- function(P,s,t) solve(P[,,s],P[,,t])
fun(temp,10,60)
fun2 <- function(P,inverses,s,t) inverses[,,s] %*% P[,,t]
fun2(temp,inverses,10,60)
fun(temp,10,60)/fun2(temp,inverses,10,60)
fun3 <- function(P,inverses,s,t,i,j) sum(inverses[i,,s] * P[,j,t])
fun3(temp,inverses,10,60,1,4)


optimObjective <- function(stratum,sp01=0.1,sp12=1,hessian=FALSE) {
    nterm01 <- nterm12 <- 5
    pnegll <- function(beta) {
        int01 <- beta[i <- 1]
        int12 <- beta[i <- i+1]
        beta01 <- beta[(i+1):(i <- i+nterm01+2)]
        beta12 <- beta[(i+1):(i <- i+nterm12+2)]
        beta20 <- beta[i <- i+1]
        ## print(beta,digits=3)
        do.call("sum",
                mclapply(stratum, function(obj) {
                    smoking <- obj$smoking
                    mort <- obj$mort
                    .Call("gsl_main2ReclassifiedPS",
                          smoking$smkstat - 1, 
                          as.integer(smoking$recall),
                          smoking$agestart,
                          smoking$agequit,
                          smoking$age, # age observed
                          smoking$freq, # frequency (as double)
                          int01,
                          int12,
                          10, 40, nterm01, attr(pspline(c(10,40),nterm=nterm01),"pparm"), sp01,
                          10, 70, nterm12, attr(pspline(c(10,70),nterm=nterm12),"pparm"), sp12,
                          beta01,
                          beta12,
                          beta20,
                          mort$age,
                          mort$Never,
                          mort$Current,
                          mort$Former,
                          FALSE, # debug
                          package="purged")$pnegll
                }, mc.cores=2))
    }
    negll <- function(beta) {
        int01 <- beta[i <- 1]
        int12 <- beta[i <- i+1]
        beta01 <- beta[(i+1):(i <- i+nterm01+2)]
        beta12 <- beta[(i+1):(i <- i+nterm12+2)]
        beta20 <- beta[i <- i+1]
        do.call("sum",
                mclapply(stratum, function(obj) {
                    smoking <- obj$smoking
                    mort <- obj$mort
                    .Call("gsl_main2ReclassifiedPS",
                          smoking$smkstat - 1, 
                          as.integer(smoking$recall),
                          smoking$agestart,
                          smoking$agequit,
                          smoking$age, # age observed
                          smoking$freq, # frequency (as double)
                          int01,
                          int12,
                          10, 40, nterm01, attr(pspline(c(10,40),nterm=nterm01),"pparm"), sp01,
                          10, 70, nterm12, attr(pspline(c(10,70),nterm=nterm12),"pparm"), sp12,
                          beta01,
                          beta12,
                          beta20,
                          mort$age,
                          mort$Never,
                          mort$Current,
                          mort$Former,
                          FALSE, # debug
                          package="purged")$negll
                }, mc.cores=2))
    }
    ## return(pnegll(init))
    out <- nlminb(init,pnegll,control=list(trace=1,rel.tol=1e-7))
    out$negll <- negll # function 
    out$pnegll <- pnegll # function
    out$coefficients <- out$par
    if (hessian) {
        out$hessian <- optimHess(coef(out), pnegll)
        ## Hl <- numDeriv::hessian(llike,coef)
        Hl <- optimHess(coef(out), negll)
        Hinv <- solve(out$hessian)
        out$trace <- sum(diag(Hinv %*% Hl))
    }
    return(out)
}
init <- c(-3,
          -4,
          rep(0.1,5+2),
          rep(0.1,5+2),
          log(0.01))
## init <- c(-2.2084, -5.2397, 1.7865, 3.9082, 4.2282, 2.6938, 1.7596, 2.0776, 
##           0.8403, -0.5627, -0.9299, -1.105, -0.6869, 0.1301, 1.1833, 2.0894, 
##           -4.3141)
init <- c(-2.2403, -5.2734, 1.7934, 3.9026, 4.1966, 2.6475, 1.7147, 2.032, 
          0.7888, -0.1281, -0.2629, -0.4078, 0.0032, 0.8361, 1.9041, 2.8168, 
          -4.4402)
system.time(fit1890.1 <- optimObjective(stratifiedData[[1]],sp01=0.1,sp12=1,hessian=TRUE))
str(fit1890.1)
with(fit1890.1, pnegll(par))
with(fit1890.1, dput(round(par,4)))
diag(solve(fit1890.1$hessian))
##
## coef.nlminb <- coef(fit1890.1)
coef.nlminb <- c(-2.30426533639024, -5.59189828926138, 1.7591741291095, 3.94699944347449, 
4.05923551226308, 2.39407490385323, 1.39327811089255, 1.92043040616437, 
1.37925846582657, 0.705860430956955, 0.138782801435078, -0.970265049666329, 
0.0159183612281592, 0.0379401144182051, 0.100106028532845, 0.263734898066751, 
0.388138157647989, 0.201325379289236, -0.059788133365412, -0.13877113685924, 
-0.0194342640363815, 0.44876575129484, -4.71103159635782)
with(fit1890.1, pnegll(coef.nlminb))


## test MPI
library(Rmpi)
library(snow)
library(parallel)
library(purged)
## set up nodes and cores
NumberOfNodes <- 1
CoresPerNode <- 2
mc.cores <- max(1, NumberOfNodes*CoresPerNode-1) # minus one for master
cl <- makeMPIcluster(mc.cores)
cat(sprintf("Running with %d workers\n", length(cl)))
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
    out <- nlminb(init,objective)
    out$objective <- objective
    out$coefficients <- out$par
    if (hessian) {
        out$hessian <- optimHess(coef(out), objective)
    }
    return(out)
}
clusterCall(cl, function() {
    library(purged)
    NULL })
out <- clusterMap(cl, optimObjective, stratifiedData[1], inits[1])
save(out,file="~/src/R/purged/test/out-20140528_B.RData")
stopCluster(cl)
mpi.quit()

## ns models
refresh
library(parallel)
library(purged)
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

## range for the smoking initiation and cessation times
sapply(stratifiedData,
       function(obj) range(sapply(obj, function(obji) with(obji$smoking,agestart[recall==1])),na.rm=TRUE))
sapply(stratifiedData,
       function(obj) range(sapply(obj, function(obji) with(obji$smoking,agequit[recall==1])),na.rm=TRUE))
lapply(stratifiedData,
       function(obj) {
           data <- do.call("rbind", lapply(obj, function(obji) subset(obji$smoking,recall==1)))
           i <- rep(1:nrow(data), data$freq)
           list(start=table(data$agestart[i],useNA="no"),
                knots1=quantile(data$agestart[i],c(0.025,0.5,0.975),na.rm=TRUE),
                quit=table(data$agequit[i],useNA="no"),
                knots2=quantile(data$agequit[i],c(0.025,0.5,0.975),na.rm=TRUE))
       })

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
                          TRUE, # debug
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
##system.time(out <- mclapply(1:42, function(i) {cat("Object",i,"\n"); try(optimObjective(stratifiedData[[i]], inits[[i]]))}, mc.cores=2))
##save(out,file="~/src/R/purged/test/out-20150502-full.RData")
load("~/src/R/purged/test/out-20150502-full.RData")

## Continuing from above
load("~/src/R/purged/test/out-20150502-full.RData")
predictPij <- function(obj,beta) {
        int01 <- beta[i <- 1]
        int12 <- beta[i <- i+1]
        beta01 <- beta[(i+1):(i <- i+2)]
        beta12 <- beta[(i+1):(i <- i+2)]
        beta20 <- beta[i <- i+1]
        smoking <- obj$smoking
        mort <- obj$mort
        .Call("gsl_predReclassified",
                      max(smoking$age),
                          smoking$smkstat - 1, 
                          smoking$age, # age observed
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
                          TRUE, # debug
                          package="purged")
    }

with(list(age=10:80), {
    P <- matrix(predictPij(list(mort=stratifiedData[[1]][[3]]$mort,smoking=expand.grid(smkstat=1:4,age=age)),out[[1]]$par),ncol=4,byrow=TRUE)
    P <- t(apply(P,1,function(row) row/sum(row)))
    P <- cbind(P[,1]+P[,4],P)
    matplot(age,P,type="l")
    legend("topright",legend=c("Never (self-report)","Never (actual)","Current","Former","Reclassified"),col=1:5,lty=1:5,bty="n")
})

require(akima)
require(dplyr)
prev <- with(list(age=10:50),
             lapply(1:42, function(i) {
                 P <- matrix(predictPij(list(mort=stratifiedData[[i]][[3]]$mort,smoking=expand.grid(smkstat=1:4,age=age)),out[[i]]$par),ncol=4,byrow=TRUE)
                 P <- t(apply(P,1,function(row) row/sum(row)))
                 data.frame(cohort=stratifiedData[[i]][[3]]$cohort,ages=age,sex=ifelse(i<=21,1,2),current=P[,2],former=P[,3])
             }))
prev.1 <- do.call("rbind",prev) %>% mutate(year = cohort+ages) %>% filter(sex==1 & year<=2010)
fld <- with(prev.1, interp(x = ages, y = year, z = log(current)))
fld$z <- exp(fld$z)
pdf("~/Downloads/us-smoking-current-males.pdf")
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               main="Current smoking prevalence, males",
               xlab="Age (years)", ylab="Calendar period",
               color.palette =
                 colorRampPalette(c("white","blue")))
dev.off()
prev.2 <- do.call("rbind",prev) %>% mutate(year = cohort+ages) %>% filter(sex==2 & year<=2010)
fld <- with(prev.2, interp(x = ages, y = year, z = log(current)))
fld$z <- exp(fld$z)
pdf("~/Downloads/us-smoking-current-females.pdf")
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               main="Current smoking prevalence, females",
               xlab="Age (years)", ylab="Calendar period",
               color.palette =
                 colorRampPalette(c("white","blue")))
dev.off()



## delta method
require(numdelta)
males <- mapply(function(obj,data)
                {
                    obj$Sigma <- solve(obj$hessian)
                    obj$data <- data
                    structure(obj,class="purged")
                },
                out[1:21],
                stratifiedData[1:21],
                SIMPLIFY=FALSE)
vcov.purged <- function(obj) obj$Sigma
coef.purged <- function(obj) obj$par
`coef<-.purged` <- function(object,value) { object$par <- value; object }
## trans <- logit <- function(x) log(x/(1-x))
## itrans <- expit <- function(x) 1/(1+exp(-x))
trans <- function(x) log(-log(x))
itrans <- function(x) exp(-exp(x))
## itrans(trans(.9))
FUN <- function(object,age=10:80)
    trans(predictPij(list(mort=object$data[[3]]$mort,smoking=expand.grid(smkstat=3,age=age)),object$par))
with(list(age=10:80), {
    pred1 <- numdelta:::predictnl(males[[1]],FUN)
    matplot(age,itrans(cbind(pred1$fit,confint(pred1))),type="l")
})

## Use a log-log transformation for confidence intervals for the predicted probabilities
## (a logit transform did not perform as well as the log-log transform)
## Assessed using sampling from the posterior
require(mvtnorm)
sims <- lapply(1:100, function(i) {
    object <- males[[1]]
    object$par <- rmvnorm(1,object$par,object$Sigma)
    FUN(object)
})
sims2 <- do.call("rbind",sims)
sum2 <- t(apply(sims2,2,quantile,c(0.025,0.975)))
with(list(age=10:80), {
    pred1 <- predictnl(males[[1]],FUN)
    matplot(age,itrans(cbind(pred1$fit,confint(pred1))),type="l")
    lines(age,itrans(sum2[,1]),lwd=3)
    lines(age,itrans(sum2[,2]),lwd=3)
})

## Can we smooth the different results?
## Should we smooth the parameters, the rates or the predicted probabilities?
coefs <- sapply(out[1:21], function(obj) obj$par)
## The issue here is that not all of the cohorts have sensible estimates for parameters outside of the observed data (although such parameters should also have large standard errors)
ses <- sapply(males, function(obj) sqrt(diag(obj$Sigma))) # oops

## How about a quick way to pull out the transition intensities?
require(rstpm2)
load("~/src/R/purged/test/out-20150502-full.RData")
lambdaij <- function(obj,age=seq(10,80,length=301),beta=NULL) {
    beta <- obj$par
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- beta[(i+1):(i <- i+2)]
    beta12 <- beta[(i+1):(i <- i+2)]
    beta20 <- beta[i <- i+1]
    data <- data.frame(age=age)
    ns01 <- model.matrix(~1+nsx(age,knots=20,centre=20,Boundary.knots=c(10,30)),data)
    ns12 <- model.matrix(~1+nsx(age,knots=40,centre=40,Boundary.knots=c(20,60)),data)
    data.frame(age=age,
          alpha01=exp(ns01 %*% c(int01,beta01)),
          alpha12=exp(ns12 %*% c(int12,beta12)))
    }
with(lambdaij(out[[1]]), matplot(age,cbind(alpha01,alpha12),type="l"))




getStratifiedData <- function(cohort,sex) {
    offset <- 21
    index <- (cohort-1885)/5
    stratifiedData[index+offset*(sex==2)]
}
##
(fit1890.1 <- optimObjective(stratifiedData[[1]],sp01=0.1,sp12=1,maxit=50,hessian=TRUE))
##
##
stratifiedData[[21]]
(fit1890.2 <- optimObjective(stratifiedData[[1+21]],sp01=0.01,sp12=1))



gcv<-function(optim.fit){
  like <- optim.fit@like
  Hl <- numDeriv::hessian(like,coef(pstpm2.fit))
  if (any(is.na(Hl)))
      Hl <- optimHess(coef(pstpm2.fit), like)
  Hinv <- -vcov(pstpm2.fit)
  trace <- sum(diag(Hinv%*%Hl))
  l <- like(coef(pstpm2.fit))
  structure(-l+trace,negll=-l,trace=trace)
}


## display results
require(survival)
predict.pspline <- 
function (object, newx, ...) 
{
    if (missing(newx)) 
        return(object)
    a <- c(list(x = newx, penalty = FALSE), attributes(object)[c("intercept", "Boundary.knots","nterm")])
    do.call("pspline", a)
}
display <- function(fit,df01=7,df12=7) {
    beta <- fit$par
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- c(int01,beta[(i+1):(i <- i+df01)])
    beta12 <- c(int12,beta[(i+1):(i <- i+df12)])
    beta20 <- beta[i <- i+1]
    s01 <- pspline(c(10,40),nterm=df01-2)
    s12 <- pspline(c(10,60),nterm=df12-2)
    if (length(fit$hessian)>0) {
        sigma <- solve(fit$hessian)
        sigma01 <- sigma[i <- c(1,3:(2+df01)),i]
        i <-  c(2,(3+df01):(2+df01+df12))
        sigma12 <- sigma[i,i]
        sigma20 <- sigma[7,7] # ??
    }
    sfun <- function(obj,centre) {
        nsref <- predict(obj,centre)
        fun <- function(t,beta) apply(predict(obj,t),1,
                                      function(row) exp(beta[1]+sum((row-nsref) * beta[-1])))
        fun
    }
    nsXfun <- function(obj,centre) {
        nsref <- predict(obj,centre)
        fun <- function(t) apply(predict(obj,t),1,
                                      function(row) c(1,row-nsref))
        fun
    }
    alpha01 <- sfun(s01,centre=20)
    alpha12 <- sfun(s12,centre=40)
    X01 <- nsXfun(s01,centre=20)
    X12 <- nsXfun(s12,centre=40)
    ages=seq(0,100,length=301)
    a01 <- alpha01(ages,beta01)
    a12 <- alpha12(ages,beta12)
    if (length(fit$hessian)>0) {
        sd01 <- as.vector(sqrt(colSums(X01(ages)* (sigma01 %*% X01(ages)))))
        matplot(ages,a01*exp(cbind(0,-1.96*sd01,1.96*sd01)),type="l",xlim=c(0,80),
                xlab="Age (years)", ylab="Hazard", main="Smoking initiation")
        a12 <- alpha12(ages,beta12)
        sd12 <- as.vector(sqrt(colSums(X12(ages)* (sigma12 %*% X12(ages)))))
        matplot(ages,a12*exp(cbind(0,-1.96*sd12,1.96*sd12)),type="l",xlim=c(0,80),
                xlab="Age (years)", ylab="Hazard", ylim=c(0,min(1,max(a12))),
                main="Smoking cessation")
        alpha20 <- list(est=exp(beta20),lower=exp(beta20-1.96*sqrt(sigma20)),
                        upper=exp(beta20+1.96*sqrt(sigma20)))
        cat(sprintf("alpha20=%f (95%% CI: %f, %f)\n",alpha20$est,alpha20$lower,alpha20$upper))
        invisible(list(ages=ages,a01=a01,sd01=sd01,a12=a12,sd12=sd12,alpha20=alpha20))
    } else {
        plot(ages,a01,type="l",xlim=c(0,80),
             xlab="Age (years)", ylab="Hazard", main="Smoking initiation")
        plot(ages,a12,type="l",xlim=c(0,80),
             xlab="Age (years)", ylab="Hazard", ylim=c(0,min(1,max(a12))),
             main="Smoking cessation")
    }
}
par(mfrow=c(2,2))
##debug(display)
display(fit1890.1)
display(fit1890.2)





library(Rmpi)
library(snow)
library(parallel)
require(purged)
##load(file="~/Documents/clients/ted/smokingList-20140528.RData")
load(file="~/src/R/purged/test/smokingList-20140528.RData")
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
        optim(init,objective,control=list(trace=0),hessian=TRUE)
    }
    clusterApply(cl,strata, optimObjective)
}
out <- do.all(stratifiedData)
save("~/src/R/purged/test/out-20140528.RData")
##objective(init)
stopCluster(cl)
##mpi.quit()


library(Rmpi)
library(snow)
library(parallel)
require(purged)
##load(file="~/Documents/clients/ted/smokingList-20140528.RData")
load(file="~/src/R/purged/test/smokingList-20140528.RData")
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
    optimObjective <- function(data) {
        objective <- function(beta) {
            int01 <- beta[i <- 1]
            int12 <- beta[i <- i+1]
            beta01 <- beta[(i+1):(i <- i+2)]
            beta12 <- beta[(i+1):(i <- i+2)]
            beta20 <- beta[i <- i+1]
            do.call("sum",
                    clusterApply(cl, data, function(obj) {
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
        optim(init,objective,control=list(trace=2,maxit=150),hessian=TRUE)
    }
    lapply(strata, optimObjective)
}
out <- do.all(stratifiedData)
out1990 <- do.all(stratifiedData[c(21,42)])
##out32 <- do.all(stratifiedData[32])
##out11 <- do.all(stratifiedData[11])
##out1910 <- do.all(stratifiedData[c(5,26)])
save(out32,out11,out1910,"~/src/R/purged/test/out-20140528.RData")
##objective(init)
stopCluster(cl)
##mpi.quit()

## optim1 <- optim(init,objective,control=list(trace=2),hessian=TRUE) # SLOW



## Display the estimated transition intensities
require(splines)
display <- function(fit,df01=2,df12=2) {
    beta <- fit$par
    sigma <- solve(fit$hessian)
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- c(int01,beta[(i+1):(i <- i+df01)])
    beta12 <- c(int12,beta[(i+1):(i <- i+df12)])
    beta20 <- beta[i <- i+1]
    knots01 <- c(10,20,30)
    knots12 <- c(20,40,60)
    sigma01 <- sigma[c(1,3,4),c(1,3,4)]
    sigma12 <- sigma[c(2,5,6),c(2,5,6)]
    sigma20 <- sigma[7,7]
    nsfun <- function(knots,centre=knots[1]) {
        nsobj <- ns(knots,df=length(knots)-1)
        nsref <- predict(nsobj,centre)
        fun <- function(t,beta) apply(predict(nsobj,t),1,
                                      function(row) exp(beta[1]+sum((row-nsref) * beta[-1])))
        fun
    }
    nsXfun <- function(knots,centre=knots[1]) {
        nsobj <- ns(knots,df=length(knots)-1)
        nsref <- predict(nsobj,centre)
        fun <- function(t,beta) apply(predict(nsobj,t),1,
                                      function(row) c(1,row-nsref))
        fun
    }
    alpha01 <- nsfun(knots01,centre=20)
    alpha12 <- nsfun(knots12,centre=40)
    X01 <- nsXfun(knots01,centre=20)
    X12 <- nsXfun(knots12,centre=40)
    ages=seq(0,100,length=301)
    a01 <- alpha01(ages,beta01)
    sd01 <- as.vector(sqrt(colSums(X01(ages)* (sigma01 %*% X01(ages)))))
    matplot(ages,a01*exp(cbind(0,-1.96*sd01,1.96*sd01)),type="l",xlim=c(0,40),
            xlab="Age (years)", ylab="Hazard", main="Smoking initiation")
    a12 <- alpha12(ages,beta12)
    sd12 <- as.vector(sqrt(colSums(X12(ages)* (sigma12 %*% X12(ages)))))
    matplot(ages,a12*exp(cbind(0,-1.96*sd12,1.96*sd12)),type="l",xlim=c(0,80),
            xlab="Age (years)", ylab="Hazard", ylim=c(0,min(1,max(a12))),
            main="Smoking cessation")
    alpha20 <- list(est=exp(beta20),lower=exp(beta20-1.96*sqrt(sigma20)),
                 upper=exp(beta20+1.96*sqrt(sigma20)))
    cat(sprintf("alpha20=%f (95%% CI: %f, %f)\n",alpha20$est,alpha20$lower,alpha20$upper))
    invisible(list(ages=ages,a01=a01,sd01=sd01,a12=a12,sd12=sd12,alpha20=alpha20))
}
## display results from HPC
## system("ssh mc2495@omega.hpc.yale.edu `cd ~/src/R/purged/test; ./submit_cluster.sh`")
## system("scp mc2495@omega.hpc.yale.edu:~/src/R/purged/test/out-20140528.RData ~/src/R/purged/test/")
##load("~/src/R/purged/test/out-20140528.RData")
##load("~/src/R/purged/test/out-20150502-full.RData")
load("~/src/R/purged/test/out-20150818-ns.RData")
out <- lapply(1:42, function(j) { new <- lapply(1:8,function(i) out[i,j][[1]]); names(new) <- rownames(out); new})
load(file="~/src/R/purged/test/smokingList-20140528.RData")
strata <- transform(expand.grid(from=seq(1890,1990,by=5),sex=1:2),
                    to=from+4)
smokingSex <- sapply(smokingList,function(obj) obj$sex)
smokingCohort <- sapply(smokingList,function(obj) obj$cohort)
## cbind(strata,t(sapply(out,function(obj) display(obj)$alpha20)))
summ <- do.call("rbind",
        lapply(1:length(out),function(i)
               with(display(out[[i]]),
                    data.frame(one=1,
                               sex=strata$sex[i],
                               cohort=strata$from[i],
                               ages=ages,
                               a01=a01,
                               lower01=a01*exp(-1.96*sd01),
                               upper01=a01*exp(1.96*sd01),
                               a12=a12,
                               lower12=a12*exp(-1.96*sd12),
                               upper12=a12*exp(1.96*sd12)))))
summ.1 <- subset(summ,sex==1)
summ.2 <- subset(summ,sex==2)
summ.1 <- subset(summ,sex==1 & cohort<=1950)
summ.2 <- subset(summ,sex==2 & cohort<=1950)
par(mfrow=c(2,2))
display(out[[5]])
display(out[[21+5]])


require(lattice)
require(mgcv)
require(dplyr)
d <- mutate(summ.1,year=cohort+ages,loga01=log(a12)) %>% filter(year<=2010 & ages>=10)
fit1 <- gam(loga01~s(ages,year),data=d,sp=1e-2); fit1$sp
plot(fit1, pers=TRUE)
plot(fit1,se=FALSE)
pred1 <- data.frame(d, pred=predict(fit1,type="response"))

library(akima)

# interpolation
d <- mutate(summ.1,year=cohort+ages) %>% filter(year<=2010 & ages>=10 & ages<95)
fld <- with(subset(d, ages<45), interp(x = ages, y = year, z = log(a01)))
fld$z <- exp(fld$z)
pdf("~/Downloads/us-smoking-uptake-males.pdf")
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               main="Smoking uptake rates, males",
               xlab="Age (years)", ylab="Calendar period",
               color.palette =
                 colorRampPalette(c("white","blue")))
dev.off()
fld <- with(d, interp(x = ages, y = year, z = log(a12)))
## fld$z <- exp(fld$z)
pdf("~/Downloads/us-smoking-cessation-males-log.pdf")
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               main="Smoking cessation log(rates), males",
               color.palette =
                 colorRampPalette(c("white", "blue")))
dev.off()

d <- mutate(summ.2,year=cohort+ages) %>% filter(year<=2010 & ages>=10 & ages<95)
fld <- with(subset(d, ages<45), interp(x = ages, y = year, z = log(a01)))
fld$z <- exp(fld$z)
pdf("~/Downloads/us-smoking-uptake-females.pdf")
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               main="Smoking uptake rates, females",
               xlab="Age (years)", ylab="Calendar period",
               color.palette =
                 colorRampPalette(c("white","blue")))
dev.off()
fld <- with(d, interp(x = ages, y = year, z = log(a12)))
##fld$z <- exp(fld$z)
pdf("~/Downloads/us-smoking-cessation-females-log.pdf")
filled.contour(x = fld$x,
               y = fld$y,
               z = fld$z,
               main="Smoking cessation log(rates), females",
               color.palette =
                 colorRampPalette(c("white", "blue")))
dev.off()

                    
pdf("~/src/R/purged/test/coplot-1.pdf")
coplot(one ~ ages | factor(cohort),
         data=summ.1,
         subscripts=T,
         panel=function(x,y,subscripts,...){
           z <- summ.1[subscripts,]
           lines(z$ages,z$a01, ...)
           lines(z$ages,z$lower01, lty=2, ...)
           lines(z$ages,z$upper01, lty=2, ...)
         }, type="l", ylim=c(0,.2), xlim=c(10,35), xlab="Age (years)",
         ylab="Hazard")
dev.off()
pdf("~/src/R/purged/test/coplot-2.pdf")
coplot(one ~ ages | factor(cohort),
         data=summ.2,
         subscripts=T,
         panel=function(x,y,subscripts,...){
           z <- summ.2[subscripts,]
           lines(z$ages,z$a01, ...)
           lines(z$ages,z$lower01, lty=2, ...)
           lines(z$ages,z$upper01, lty=2, ...)
         }, type="l", ylim=c(0,.2), xlim=c(10,35), xlab="Age (years)",
         ylab="Hazard")
dev.off()
pdf("~/src/R/purged/test/coplot-3.pdf")
coplot(one ~ ages | factor(cohort),
         data=summ.1,
         subscripts=T,
         panel=function(x,y,subscripts,...){
           z <- summ.1[subscripts,]
           lines(z$ages,z$a12, ...)
           lines(z$ages,z$lower12, lty=2, ...)
           lines(z$ages,z$upper12, lty=2, ...)
         }, type="l", ylim=c(0,.15), xlim=c(10,70), xlab="Age (years)",
         ylab="Hazard")
dev.off()
pdf("~/src/R/purged/test/coplot-4.pdf")
coplot(one ~ ages | factor(cohort),
         data=summ.2,
         subscripts=T,
         panel=function(x,y,subscripts,...){
           z <- summ.2[subscripts,]
           lines(z$ages,z$a12, ...)
           lines(z$ages,z$lower12, lty=2, ...)
           lines(z$ages,z$upper12, lty=2, ...)
         }, type="l", ylim=c(0,.15), xlim=c(10,70), xlab="Age (years)",
         ylab="Hazard")
dev.off()

## scratch
display(out[[1]])
display(out[[1+21]])
plot.1 <- display(out[[1]])
plot.2 <- display(out[[1+21]])
plot.1 <- display(out[[11]])
plot.2 <- display(out[[11+21]])
with(plot.1, plot(ages,a01,type="l",xlim=c(0,40)))
with(plot.2, lines(ages,a01,lty=2))
with(plot.1, plot(ages,a12,type="l"))
with(plot.2, lines(ages,a12,lty=2))




## Identifiability constraint for APC model
require(mgcv)
require(rstpm2)
somedata <- within(expand.grid(age=seq(0,80,by=5.0),
                                  year=seq(1980,2010,by=1.0)),
                      { cohort <- year-age
                        y <- rnorm(length(age))
                    })
gam1 <- gam(y~s(age)+s(cohort)+s(year),
            data=somedata)
##
##
## bsplinepen
require(splines)
require(fda)
basisobj <- fda::create.bspline.basis(c(0,1),13)
x=seq(0,1,length=11)
bs(x,knots=x[2:10],int=T) - predict(basisobj)
##
##
ns.Q <- 
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x)) 
{
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
            Boundary.knots[2L])
    }
    else outside <- FALSE
    if (!is.null(df) && is.null(knots)) {
        nIknots <- df - 1L - intercept
        if (nIknots < 0L) {
            nIknots <- 0L
            warning(gettextf("'df' was too small; have used %d", 
                1L + intercept), domain = NA)
        }
        knots <- if (nIknots > 0L) {
            knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L, 
                nIknots + 2L)]
            stats::quantile(x[!outside], knots)
        }
    }
    else nIknots <- length(knots) 
    Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
    const <- splineDesign(Aknots, Boundary.knots, 4, c(2, 2))
    if (!intercept) {
        const <- const[, -1, drop = FALSE]
    }
    qr.const <- qr(t(const))
    qmat <- qr.Q(qr.const, complete=TRUE)
    qmat[, -(1L:2L), drop = FALSE]
}
(bs(x,knots=x[2:10],int=T) %*% ns.Q(x,knots=x[2:10],int=T)) - ns(x,knots=x[2:10],int=T)
##
##
zapsmall(ns.Q(x,knots=x[2:10],int=T) %*% t(ns.Q(x,knots=x[2:10],int=T)))
bsplinepen(basisobj)

## Reimplement the example in R
require(splines)
require(deSolve)
mu0Data <- data.frame(age=seq(0.5,105.5,by=1),
                  rate=c(0.00219, 0.000304, 5.2e-05, 0.000139, 0.000141, 3.6e-05, 7.3e-05, 
		  0.000129, 3.8e-05, 0.000137, 6e-05, 8.1e-05, 6.1e-05, 0.00012, 
		0.000117, 0.000183, 0.000185, 0.000397, 0.000394, 0.000585, 0.000448, 
		0.000696, 0.000611, 0.000708, 0.000659, 0.000643, 0.000654, 0.000651, 
		0.000687, 0.000637, 0.00063, 0.000892, 0.000543, 0.00058, 0.00077, 
		0.000702, 0.000768, 0.000664, 0.000787, 0.00081, 0.000991, 9e-04, 
		0.000933, 0.001229, 0.001633, 0.001396, 0.001673, 0.001926, 0.002217, 
		0.002562, 0.002648, 0.002949, 0.002729, 0.003415, 0.003694, 0.004491, 
		0.00506, 0.004568, 0.006163, 0.006988, 0.006744, 0.00765, 0.007914, 
		0.009153, 0.010231, 0.011971, 0.013092, 0.013839, 0.015995, 0.017693, 
		0.018548, 0.020708, 0.022404, 0.02572, 0.028039, 0.031564, 0.038182, 
		0.042057, 0.047361, 0.05315, 0.058238, 0.062619, 0.074934, 0.089776, 
		0.099887, 0.112347, 0.125351, 0.143077, 0.153189, 0.179702, 0.198436, 
		0.240339, 0.256215, 0.275103, 0.314157, 0.345252, 0.359275, 0.41768, 
		0.430279, 0.463636, 0.491275, 0.549738, 0.354545, 0.553846, 0.461538, 
		0.782609))
mu1Data <- transform(mu0Data,rate=rate*1.4)
mu2Data <- transform(mu0Data,rate=rate*1.2)
mu0 <- splinefun(mu0Data$age, mu0Data$rate)
mu1 <- splinefun(mu1Data$age, mu1Data$rate)
mu2 <- splinefun(mu2Data$age, mu2Data$rate)
nsfun <- function(knots,centre=knots[1],intercept=0) {
    nsobj <- ns(knots,df=length(knots)-1)
    nsref <- predict(nsobj,centre)
    fun <- function(t,beta) apply(predict(nsobj,t),1,
                                  function(row) exp(intercept+sum((row-nsref) * beta)))
    fun
}
alpha01 <- nsfun(c(10,20,30),centre=20,intercept=-3)
alpha12 <- nsfun(c(20,40,60),centre=40,intercept=-4)
##plot(10:30,alpha01(10:30,c(1,-1)))
##plot(10:80,alpha12(10:80,c(1,-1)))
bounds <- function(t,lower,upper) pmax(pmin(t,upper),lower)
parms <- c(beta20=log(0.01),beta01=c(1,-1),beta12=c(1,-1)) # GLOBAL
Never <- 1; Current <- 2; Former <- 3; Reclassified <- 4; Death <- 5
odefunc <- function(t,y,parms) {
    a20 <- exp(parms["beta20"])
    a01 <- alpha01(t,c(parms["beta011"],parms["beta012"]))
    a12 <- alpha12(t,c(parms["beta121"],parms["beta122"]))
    m0 <- mu0(bounds(t,0.5,105.5))
    m1 <- mu1(bounds(t,0.5,105.5))
    m2 <- mu2(bounds(t,0.5,105.5))
    list(c(-y[Never]*(m0+a01),
           -y[Current]*(m1+a12)+y[Never]*a01,
           -y[Former]*(m2+a20)+y[Current]*a12,
           -y[Reclassified]*m0+y[Former]*a20),
         c(a01=a01,a12=a12,a20=a20,m0=m0,m1=m1,m2=m2))
}     
muij <- function(i,j,t,parms) {
    if (i==Never & j==Current) return(alpha01(t,c(parms["beta011"],parms["beta012"])))
    if (i==Current & j==Former) return(alpha12(t,c(parms["beta121"],parms["beta122"])))
    if (i==Former & j==Reclassified) return(exp(parms["beta20"]))
    if (j==Death) {
        m0 <- mu0(bounds(t,0.5,105.5))
        m1 <- mu1(bounds(t,0.5,105.5))
        m2 <- mu2(bounds(t,0.5,105.5))
        if (i==Never) return(m0)
        if (i==Current) return(m1)
        if (i==Former) return(m2)
        if (i==Reclassified) return(m0)
    }
    return(0)
    }
Pij <- function(i,j,t0,t1,parms,size=4) {
    y <- rep(0,size)
    y[i] <- 1
    ode(y=y, t=c(t0,t1), func = odefunc, parms,
            rtol=1e-8, atol=1e-8)[2,j+1]
}
PiK <- function(i,t0,t1,parms,size=4) {
    y <- rep(0,size)
    y[i] <- 1
    sum(ode(y=y, t=c(t0,t1), func = odefunc, parms,
            rtol=1e-8, atol=1e-8)[2,2:(size+1)])
}
ll <- function(finalState,s,t,u) {
    if(finalState==Never) return(log((Pij(Never,Never,0,u,parms)+Pij(Never,Reclassified,0,u,parms))/PiK(Never,0,u,parms)))
    if(finalState==Current)
        return(log(Pij(Never,Never,0,s,parms)*muij(Never,Current,s,parms)*
                   Pij(Current,Current,s,u,parms)/PiK(Never,0,u,parms)))
    if(finalState==Former)
        return(log(Pij(Never,Never,0,s,parms)*muij(Never,Current,s,parms)*
                   Pij(Current,Current,s,t,parms)*muij(Current,Former,t,parms)*
                   Pij(Former,Former,t,u,parms)/PiK(Never,0,u,parms)))
}
ll(Never,u=70)
ll(Current,20,u=70)
ll(Former,20,50,u=70)

## Simulated data
require(msm)
require(purged)
mu0Data <- data.frame(age=seq(0.5,105.5,by=1),
                  rate=c(0.00219, 0.000304, 5.2e-05, 0.000139, 0.000141, 3.6e-05, 7.3e-05, 
		  0.000129, 3.8e-05, 0.000137, 6e-05, 8.1e-05, 6.1e-05, 0.00012, 
		0.000117, 0.000183, 0.000185, 0.000397, 0.000394, 0.000585, 0.000448, 
		0.000696, 0.000611, 0.000708, 0.000659, 0.000643, 0.000654, 0.000651, 
		0.000687, 0.000637, 0.00063, 0.000892, 0.000543, 0.00058, 0.00077, 
		0.000702, 0.000768, 0.000664, 0.000787, 0.00081, 0.000991, 9e-04, 
		0.000933, 0.001229, 0.001633, 0.001396, 0.001673, 0.001926, 0.002217, 
		0.002562, 0.002648, 0.002949, 0.002729, 0.003415, 0.003694, 0.004491, 
		0.00506, 0.004568, 0.006163, 0.006988, 0.006744, 0.00765, 0.007914, 
		0.009153, 0.010231, 0.011971, 0.013092, 0.013839, 0.015995, 0.017693, 
		0.018548, 0.020708, 0.022404, 0.02572, 0.028039, 0.031564, 0.038182, 
		0.042057, 0.047361, 0.05315, 0.058238, 0.062619, 0.074934, 0.089776, 
		0.099887, 0.112347, 0.125351, 0.143077, 0.153189, 0.179702, 0.198436, 
		0.240339, 0.256215, 0.275103, 0.314157, 0.345252, 0.359275, 0.41768, 
		0.430279, 0.463636, 0.491275, 0.549738, 0.354545, 0.553846, 0.461538, 
		0.782609))
mu1Data <- transform(mu0Data,rate=rate*1.4)
mu2Data <- transform(mu0Data,rate=rate*1.2)
rpexpSlow <- function(n, rate, t, t0 = 0) {
    stopifnot(length(rate)==length(t) & t[1]==0 & all(diff(t)>0) & all(rate>=0))
    eps <- .Machine$double.xmin
    if (length(t0)==1)
        t0 <- rep(t0, length=n)
    y <- vector("numeric",n)
    maxt <- t[length(t)]
    ratemaxt <- rate[length(t)]
    if (length(i <- which(t0 >= maxt))>0)
        y[i] <- t0[i] + rexp(length(i),ratemaxt)
    if (length(i <- which(t0 < maxt))>0) 
        for (t0i in unique(t0[i])) {
            j <- which(t0 == t0i)
            if (t0i==0) {
                y[j] <- msm::rpexp(length(j),rate,t)
            } else {
                k <- which(t>=t0i)
                if (any(t==t0i)) {
                    ti <- c(0,t[k])
                    ratei <- c(eps,rate[k])
                } else {
                    ti <- c(0,t0i,t[k])
                    ratei <- c(eps, rate[c(k[1]-1,k)])
                }
                y[j] <- msm::rpexp(length(j), ratei, ti)
            }
        }
    y
}
rpexp <- function(n, rate, t, t0 = 0) {
    if (length(t0)==1) t0 <- rep(t0, length=n)
    y <- msm::rpexp(n, rate, t)
    counter <- 1
    ## i <- n <- NULL
    while ((n <- length(i <- which(y<t0)))>0 && counter<10) {
        counter = counter + 1
        y[i] <- msm::rpexp(n,rate,t)
    }
    if (n>0) y[i] <- rpexpSlow(n,rate,t,t0[i])
    structure(y,counter=counter)
}
rpexpTest <- function(n, rate, t, t0=0) {
    stopifnot(length(rate)==length(t) & t[1]==0 & all(diff(t)>0) & all(rate>=0))
    if (length(t0)==1)
        t0 <- rep(t0, length=n)
    .Call("purged_rpexp",as.integer(n),t, rate, t0, package="purged")
}
if (test <- FALSE) {
    set.seed(12345)
    system.time(plot(density(y1 <- rpexp(10000,c(0.1,0.2),c(0,1),30),from=30)))
    system.time(plot(density(y2 <- rpexpTest(10000,c(0.1,0.2),c(0,1),30),from=30)))
    t.test(y1,y2)
    qqplot(y1,y2)
    abline(0,1)
}
## TODO: invert an arbitrary monotone function
set.seed(12345)
n <- 10000
test0 <- data.frame(finalState=0,t1=NA,t2=NA,t3=floor(runif(n,30,80)),
                    death1=NA,death2=NA,death3=NA)
test0 <- transform(test0, death1=rpexp(n,mu0Data$rate,mu0Data$age-0.5))
test0 <- subset(test0,t3<death1)
##
test1 <- data.frame(finalState=1,t1=NA,t2=NA,t3=floor(runif(n,30,80)),
                    death1=NA,death2=NA,death3=NA)
test1 <- transform(test1,t1 = floor(runif(length(t1),12,30)))
test1 <- transform(test1, death1=rpexp(n,mu0Data$rate,mu0Data$age-0.5))
test1 <- transform(test1, death2=rpexp(n,mu0Data$rate*1.4,mu0Data$age-0.5,t1))
test1 <- subset(test1,t1<death1 & t3<death2)
##
test2 <- data.frame(finalState=2,t1=NA,t2=NA,t3=floor(runif(n,30,80)),
                    death1=NA,death2=NA,death3=NA)
test2 <- transform(test2,t1 = floor(runif(length(t1),12,30)))
test2 <- transform(test2,t2=floor(runif(length(t3),t1,t3)))
test2 <- transform(test2, death1=rpexp(n,mu0Data$rate,mu0Data$age-0.5))
test2 <- transform(test2, death2=rpexp(n,mu0Data$rate*1.4,mu0Data$age-0.5,t1))
test2 <- transform(test2, death3=rpexp(n,mu0Data$rate*1.2,mu0Data$age-0.5,t2))
test2 <- within(test2, {
    reclassified <- rbinom(n,1,1-exp(-t2*0.01))
    finalState <- ifelse(reclassified==1,0,2)
    reclassified <- NULL
})
test2 <- subset(test2,t1<death1 & t2<death2 & t3<death3)
##
test <- rbind(test0,test1,test2)
test$freq <- 1.0
test$recall <- 1L

objective <- function(beta) {
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- beta[(i+1):(i <- i+2)]
    beta12 <- beta[(i+1):(i <- i+2)]
    ## beta20 <- log(0.01)
    beta20 <- beta[i <- i+1]
    print(beta)
    .Call("gsl_main2Reclassified",
          test$finalState, # finalState
          test$recall,
          test$t1, # t1 = age started smoking
          test$t2, # t2 = age quit
          test$t3, # t3 = age observed
          test$freq, # freq = frequency (as double)
          int01,
          int12,
          c(10,20,30), # knots01
          c(20,40,60), # knots12
          beta01,
          beta12,
          beta20,
          mu0Data$age,
          mu0Data$rate,
          mu1Data$rate,
          mu2Data$rate,
          FALSE, # debug
          package="purged")
}
##objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)))
system.time(print(1-100960.6/objective(init <- c(-3,-4,1,-1,1,-1,log(0.01))))) # should be small
##init <- c(-2.3633004, -3.7141737,  4.4504238, -0.1386734, -0.3805543, -0.4052293, log(0.01))
##options(width=120)
##optim1 <- optim(init,objectiveReclassified,control=list(trace=2)) # SLOW

## RcppExport SEXP 
## gsl_main2ReclassifiedPS(SEXP _finalState,
##       		  SEXP _recall,
##       		  SEXP _time1, SEXP _time2, SEXP _time3,
##       		  SEXP _freq,
##       		  SEXP _int01, SEXP _int12,
##       		  SEXP _lower01, SEXP _upper01, SEXP _nterm01, SEXP _pmatrix01, SEXP _sp01,
##       		  SEXP _lower12, SEXP _upper12, SEXP _nterm12, SEXP _pmatrix12, SEXP _sp12,
##       		  SEXP _beta01, SEXP _beta12,
##       		  SEXP _beta20,
##       		  SEXP _ages0,
##       		  SEXP _mu0,
##       		  SEXP _mu1,
##       		  SEXP _mu2,
##       		  SEXP _debug)
require(survival)
objectivePS <- function(beta,df01=10L,df12=10L) {
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- beta[(i+1):(i <- i+df01)]
    beta12 <- beta[(i+1):(i <- i+df12)]
    ## beta20 <- log(0.01)
    beta20 <- beta[i <- i+1]
    nterm01 <- as.integer(df01 - 2L)
    nterm12 <- as.integer(df12 - 2L)
    print(beta)
    .Call("gsl_main2ReclassifiedPS",
          test$finalState, # finalState
          test$recall,
          test$t1, # t1 = age started smoking
          test$t2, # t2 = age quit
          test$t3, # t3 = age observed
          test$freq, # freq = frequency (as double)
          int01,
          int12,
          10, 60, nterm01, attr(pspline(c(10,60),nterm=nterm01),"pparm"), 1,
          10, 30, nterm12, attr(pspline(c(10,30),nterm=nterm12),"pparm"), 1,
          beta01,
          beta12,
          beta20,
          mu0Data$age,
          mu0Data$rate,
          mu1Data$rate,
          mu2Data$rate,
          FALSE, # debug
          package="purged")
}
init <- c(-3,
          -4,
          rnorm(10),
          rep(0.1,8+2),
          log(0.01))
system.time(print(objectivePS(init)))


require(parallel)
## split into blocks for computations
objective <- function(beta,mc.cores=getOption("mc.cores",2),test) {
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- beta[(i+1):(i <- i+2)]
    beta12 <- beta[(i+1):(i <- i+2)]
    beta20 <- log(0.01)
    ## beta20 <- beta[i <- i+1]
    print(beta)
    n <- nrow(test)
    i <- sort((0:(n-1)) %% mc.cores)
    do.call("sum",
            mclapply(0:(mc.cores-1), function(core,test) {
                j <- which(i==core)
                data <- test[j,]
                .Call("gsl_main2Reclassified",
                      data$finalState, # finalState
                      data$t1, # t1 = age started smoking
                      data$t2, # t2 = age quit
                      data$t3, # t3 = age observed
                      int01,
                      int12,
                      c(10,20,30), # knots01
                      c(20,40,60), # knots12
                      beta01,
                      beta12,
                      beta20,
                      mu0Data$age,
                      mu0Data$rate,
                      package="purged")
            }, test))
}
system.time(objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)),mc.cores=1,test))
system.time(objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)),mc.cores=2,test))
system.time(objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)),mc.cores=3,test))
system.time(objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)),mc.cores=4,test))

require(parallel)
system.time(mclapply(1:10, function(i,data) objective(init), mc.cores=2))


.Call("gsl_main2Reclassified",
      c(2,0), # finalState
      c(20,NA), # t1 = age started smoking
      c(50,NA), # t2 = age quit
      c(70,70), # t3 = age observed
      c(1.0, 1.0), # freq
      -3, # int01
      -4, # int12
      c(10,20,30), # knots01
      c(20,40,60), # knots12
      c(1,-1),     # beta01
      c(1,-1),     # beta12
      log(0.01),   # beta20,
      mu0Data$age,
      mu0Data$rate,
      mu1Data$rate,
      mu2Data$rate,
      TRUE, # debug
      package="purged")




