
## We need some simulated data
## We also need to use the GSL C/C++ code through Rcpp
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
test2 <- subset(test2,t1<death1 & t2<death2 & t3<death3)
##
test <- rbind(test0,test1,test2)


objective <- function(beta) {
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- beta[(i+1):(i <- i+2)]
    beta12 <- beta[(i+1):(i <- i+2)]
    beta20 <- log(0.01)
    ## beta20 <- beta[i <- i+1]
    print(beta)
    .Call("gsl_main2",
          test$finalState, # finalState
          test$t1, # t1 = age started smoking
          test$t2, # t2 = age quit
          test$t3, # t3 = age observed
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
}
##objective(init <- c(-3,-4,1,-1,1,-1,log(0.01)))
print(objective(init <- c(-3,-4,1,-1,1,-1)))

objectiveReclassified <- function(beta) {
    int01 <- beta[i <- 1]
    int12 <- beta[i <- i+1]
    beta01 <- beta[(i+1):(i <- i+2)]
    beta12 <- beta[(i+1):(i <- i+2)]
    beta20 <- log(0.01)
    ## beta20 <- beta[i <- i+1]
    print(beta)
    .Call("gsl_main2Reclassified",
          test$finalState, # finalState
          test$t1, # t1 = age started smoking
          test$t2, # t2 = age quit
          test$t3, # t3 = age observed
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
}
print(objectiveReclassified(init <- c(-3,-4,1,-1,1,-1)))
