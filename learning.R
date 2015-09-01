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

EV1(0, 0.25)
plot(0:100,EV1(0:100,0.25), type = 'l')
#

pAtt <- function(stage, m, theta){
  stage <- as.character(stage)
  ms <- 0.2*m
  switch(stage, 
         '5' =  pexp(ms, theta),
         '4' = pexp(ms+EV5(m+5, theta) - EV5(m,theta), theta),
         '3' = pexp(ms+EV4(m+5, theta) - EV4(m,theta), theta),
         '2' = pexp(ms+EV3(m+5, theta) - EV3(m,theta), theta),
         '1' = pexp(ms+EV2(m+5, theta) - EV2(m,theta), theta))
}

pAtt(5, 0,0.25)
sapply(1:5, pAtt,80, 0.25)
plot(1:5, sapply(1:5, pAtt, 50, 0.25))
plot(0:20,EV4(0:20,0.25), type = 'l', col = 'pink', lwd=5, ylim = c(4,10))
lines(0:20,EV5(0:20,0.25), type = 'l', col = 'lightblue', lwd=5)



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
