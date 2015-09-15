# EVT
## for each possible state, comput expected max value
## 

## Supose exponential alternative benefits
## http://stats.stackexchange.com/questions/48496/conditional-expectation-of-exponential-random-variable
theta <- 0.25

m <- 40:50
ms <- 0.2*m
patt <- pexp(ms, theta)
EV <- patt*ms + (1- patt)*(ms + 1/theta)

EV5 <- function(m, theta){
  ms <- 0.2*m
  patt <- pexp(ms, theta)
  EV <- patt*ms + (1- patt)*(ms + 1/theta)
 # cat('in EV5, patt = ',patt,'m = ',m,'\n')
  return(EV)
}


EV4 <- function(m, theta){
  ms <- 0.2*m
  patt <- pexp(ms+EV5(m+5, theta) - EV5(m,theta), theta)
  EV <- ms + EV5(m +5, theta)+ (1- patt)*(1/theta)
  return(EV)
}

EV3 <- function(m, theta){
  ms <- 0.2*m
  patt <- pexp(ms+EV4(m+5, theta) - EV4(m,theta), theta)
  EV <- ms + EV4(m+5, theta)+ (1- patt)*(1/theta)
  return(EV)
}

EV2 <- function(m, theta){
  ms <- 0.2*m
  patt <- pexp(ms+EV3(m+5, theta) - EV3(m,theta), theta)
  EV <- ms + EV3(m+5, theta)+ (1- patt)*(1/theta)
  return(EV)
}

EV1 <- function(m, theta){
  ms <- 0.2*m
  patt <- pexp(ms+EV2(m+5, theta) - EV2(m,theta), theta)
  EV <- ms + EV2(m+5, theta)+ (1- patt)*(1/theta)
  return(EV)
}

EV1(0, theta)
plot(0:100,EV1(0:100,theta), type = 'l')
#

pAtt <- function(stage, m, theta, delta = 5){
  stage <- as.character(stage)
  ms <- 0.2*m
  switch(stage, 
         '5' =  pexp(ms, theta),
         '4' = pexp(ms+EV5(m+delta, theta) - EV5(m,theta), theta),
         '3' = pexp(ms+EV4(m+delta, theta) - EV4(m,theta), theta),
         '2' = pexp(ms+EV3(m+delta, theta) - EV3(m,theta), theta),
         '1' = pexp(ms+EV2(m+delta, theta) - EV2(m,theta), theta))
}

pAtt(5, 0,theta)
sapply(1:5, pAtt,80, theta)
plot(1:5, sapply(1:5, pAtt, 50, theta))
plot(0:20,EV4(0:20,theta), type = 'l', col = 'pink', lwd=5, ylim = c(4,10))
lines(0:20,EV5(0:20,theta), type = 'l', col = 'lightblue', lwd=5)


## DGP
delt <- 5
simMarks <- function(mark=50, theta = 0.25, delta = delt){
#theta <- 0.25
m <- mark
#m <- c(10,90)
ms <- 0.2*m
#stage 1
z <- rexp(m,theta)
att1 <- (z < ms+EV2(m+delta, theta) - EV2(m,theta))
m <- m + att1*delta
ms <- 0.2*m
##stage 2
z <- rexp(m,theta)
att2 <- (z < ms+EV3(m+delta, theta) - EV3(m,theta))
m <- m + att2*delta
ms <- 0.2*m
##stage 3
z <- rexp(m,theta)
att3 <- (z < ms+EV4(m+delta, theta) - EV4(m,theta))
m <- m + att3*delta
ms <- 0.2*m
##stage 4
z <- rexp(m,theta)
att4 <- (z < ms+EV5(m+delta, theta) - EV5(m,theta))
m <- m + att4*delta
ms <- 0.2*m
##stage 5
z <- rexp(m,theta)
att5 <- (z < ms)
return(cbind(att1, att2, att3, att4, att5))
}

simMarks(mark = c(10,50, 100), delta = delt, theta = theta)
set.seed(23)
mks <- c(sample(1:50,100, replace = TRUE), sample(50:70,200, replace = TRUE), 
         sample(71:80,200, replace = TRUE))
choices <- simMarks(mks, delta = delt, theta = theta)
mode(choices) <- 'numeric'
df <- data.frame(mks, choices)


ndf <- cbind(df, mks+delt*t(apply(df[,2:5],1,cumsum)))
mp <- ndf[,c(1,7:10)]
names(mp) <- paste0('m', 1:5)
dat <- cbind(mp, choices)



NLL <- function(theta, data){
  choices <- data[,c('att1', 'att2','att3','att4','att5')]
  mks <- data[,c('m1', 'm2','m3','m4','m5')]
  probs <- sapply(1:5, function(x) pAtt(x, mks[,x], theta = theta))
  -sum(log(choices*probs + (1-choices)*(1-probs)))
}

NLL(0.23, dat)

(result <- optim(0.8, NLL,  data = dat,  hessian = TRUE, method = 'BFGS'))
sqrt(1/result$hessian)
plot(x <- seq(0.1,0.35,by = 0.001), sapply(x, NLL, dat), type = 'l', lwd =4)






library(parallel)
cl <- makeCluster(detectCores()-1)  
#get library support needed to run the code
#clusterEvalQ(cl,library(repsych))
#put objects in place that might be needed for the code
#clusterExport(cl,c("myData"))
#... then parallel replicate...
parSapply(cl, 1:10000, function(i,...) { x <- rnorm(10); mean(x)/sd(x) } )
#stop the cluster
stopCluster(cl)


## updating priors
m.prior <- 50 # prior mark
m.new <- 50
#m.true <- 60
##' simas are variances
sig.prior <- 10 #prior sigma
sig.s <- 10 # signal sigma

m.post <- (sig.prior/(sig.prior + sig.s)) * m.new + (sig.s/(sig.prior+sig.s))*m.prior
sig.post <- 1/((1/sig.prior) + (1/sig.s))

dnp <- dnorm(x<-0:100, m.prior,sig.prior)
dnpost<- dnorm(x, m.post, sig.post)

plot(x , dnp, ylim = c(0,max(dnp,dnpost)) , type = 'l')
lines(x, dnpost )
abline(v = m.new)

update.prior <- function(prior.mean, prior.sigma, signal, signal.sigma){
  m.post <- (prior.sigma/(prior.sigma + signal.sigma)) * signal + (signal.sigma/(prior.sigma+signal.sigma))*prior.mean
  sig.post <- 1/((1/prior.sigma) + (1/signal.sigma))
return(c(posterior.mean = m.post, posterior.sigma = sig.post))  
}

u <- update.prior(50, 20, 50, 10)
u <- update.prior(u[1], u[2], 50, 10)
u <- update.prior(u[1], u[2], 45, 10)
u

library(emg)


## mark learning model
## use discrete grid wioth Bayesian updating
##stage 0
m <- 0:100
p <- 1/length(m)
pdf <- data.frame(m,p)
with(pdf, plot(m, p, type = 'h'))

## updated after observes ave

mn <- 45
pdf[pdf$m == mn,]$p
## prob data | prob is prob observe not 45
1 - pdf[pdf$m == mn,]$p
#so eg for prob 0
L <- pdf[pdf$m == 0,]$p
L* (1 - pdf[pdf$m == mn,]$p)

# but for 45
L <- pdf[pdf$m == 45,]$p
L * pdf[pdf$m == mn,]$p
