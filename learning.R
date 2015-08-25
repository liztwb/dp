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

EV <- function(m, theta){
  ms <- 0.2*m
  patt <- pexp(ms, theta)
  EV <- patt*ms + (1- patt)*(ms + 1/theta)
  return(EV)
}

plot(0:20,EV(0:20,0.25), type = 'l')


## updating priors
m.prior <- 50 # prior mark
m.new <- 80
#m.true <- 60
##' simas are variances
sig.prior <- 10 #prior sigma
sig.s <- 5 # signal sigma

m.post <- (sig.prior/(sig.prior + sig.s)) * m.new + (sig.s/(sig.prior+sig.s))*m.prior
sig.post <- 1/((1/sig.prior) + (1/sig.s))





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
