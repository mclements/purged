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
    class(smokingList) <- c("lookupList","list")
    subset.lookupList <- function(obj,subset) {
        e <- substitute(subset)
        obj[sapply(obj, function(obji) eval(e,obji))]
    }
    save(smokingList,mortList,subset.lookupList,file="~/Documents/clients/ted/smokingList-20140528.RData")
}

length(subset(smokingList,sex==1 & cohort == 1948))
length(subset(smokingList,sex==1))



load(file="~/Documents/clients/ted/smokingList-20140528.RData")
require(purged)
test <- subset(smokingList,sex==1 & cohort<1950 & cohort>=1945)
objective <- function(beta) {
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- beta[(i+1):(i <- i+2)]
    beta12 <- beta[(i+1):(i <- i+2)]
    beta20 <- beta[i <- i+1]
    print(beta)
    sum(sapply(test, function(obj) {
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
##objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)))
system.time(print(objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)))))

init <- c(-1.93636374725587, -3.4001339744127, 3.01006761619528, -3.02733667354355, 
1.01152056888029, 0.394911534677803, -3.98158254291634)
options(width=120)
optim1 <- optim(init,objective,control=list(trace=2),hessian=TRUE) # SLOW

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





## Display the estimated transition intensities
require(splines)
display <- function(fit) {
    beta <- fit$par
    sigma <- solve(fit$hessian)
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- c(int01,beta[(i+1):(i <- i+2)])
    beta12 <- c(int12,beta[(i+1):(i <- i+2)])
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
    matplot(ages,a01*exp(cbind(0,-1.96*sd01,1.96*sd01)),type="l")
    a12 <- alpha12(ages,beta12)
    sd12 <- as.vector(sqrt(colSums(X12(ages)* (sigma12 %*% X12(ages)))))
    matplot(ages,a12*exp(cbind(0,-1.96*sd12,1.96*sd12)),type="l")
    return(c(alpha20=exp(beta20),lower=exp(beta20-1.96*sqrt(sigma20)),
             lower=exp(beta20+1.96*sqrt(sigma20))))
}
par(mfrow=c(1,2))
##debug(display)
display(optim1)


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
system.time(print(objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)))))


init <- c(-2.3633004, -3.7141737,  4.4504238, -0.1386734, -0.3805543, -0.4052293, log(0.01))
options(width=120)
optim1 <- optim(init,objectiveReclassified,control=list(trace=2)) # SLOW





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




