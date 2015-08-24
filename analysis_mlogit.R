setwd("~/Desktop/IE/IELabs/Marks")
d <- readRDS('marks15_final.rds')
g <- readRDS('~/Desktop/IE/IELabs/Groups/sessdb.rds')
s <- readRDS('~/Desktop/IE/IELabs/sall.rds')
#library(MASS)
library(tidyr)
library(dplyr)
library(ggplot2)


labcols <- paste('Lab.',1:5,sep = '')
## make long form
dl <- d %>% dplyr::select(ID, starts_with('Lab') ) %>% gather(session, mark, Lab.1:Lab.5) #%>% head
dl %>% group_by(session) %>% dplyr::summarise(ave = mean(mark, na.rm = TRUE)) %>%   ggplot(aes(x = session, y = ave, group = 1)) + geom_line()
dl %>% ggplot(aes(x = session, y = mark)) + geom_boxplot(fill = 'lightblue') + geom_violin(alpha = 0.5)

## ID group and week
gmem <- gather(s,  wk, grp, gs1:gs5) %>% dplyr::select(ID, Course, title, grp,wk) %>% mutate(session = gsub('gs','',wk)) %>% dplyr::select(-wk)

## marks by week including group attended
mbywk <- merge(gmem,transmute(dl, ID, session = gsub('Lab.','',session),mark), by = c('ID','session'))

all <- g %>% dplyr::select(grp, sat, session = weeks) %>% merge(mbywk,by = c('grp','session'))
all <- mutate(all, rmark = ifelse(is.na(mark),0,mark))
all <- all %>% mutate(Gender = ifelse(title == 'Mr', 'M','F')) %>% dplyr::select(-title)

## Look at missed
dat <- all %>% group_by(ID, Gender) %>% dplyr::summarise( missed = sum(is.na(mark)), ave = mean(mark, na.rm = TRUE))

dat$att <- 5 - dat$missed

## Dyn Prog model

dat.dp <- dplyr::filter(dat, is.finite(ave)) %>% mutate(ave = ave +1)
m <- log(400)
sd <- log(8)

##exponentiial;
m <- 1
data <- dat.dp
ll <- data$att* log(pexp(0.2*data$ave, m)) + (5 - data$att) * log((1-pexp(0.2*data$ave,m)))

DPLL.gpd.fix <- function(theta, data){
library(evd)
    #p <- pnorm(theta[1])
p <- theta[1]
     ll <- data$att* log(p) + (5 - data$att) * log((1-p))
  return(-1*sum(ll))
}

DPLL.gpd.fix(c(.5), dat.dp)

res <- optim(c(.1), DPLL.gpd.fix, data = dat.dp, hessian = TRUE, method = 'L-BFGS-B', lower = c(0.00001), upper = c(.999999),
             control = list(maxit=500))
res
lnpar <- res$par
#pnorm(res$par)
lnhess <- res$hessian
ses <- diag(solve(lnhess))
ses
plot(x <- seq(0,20,by=.1), pgpd(x, lnpar[1],lnpar[2], lnpar[3]),  type = 'l')
plot(x <- seq(0,20,by=.1), dgpd(x, lnpar[1],lnpar[2], lnpar[3]),  type = 'l')
sqrt(diag(lnhess))
summary(rnorm(10000,lnpar[1], lnpar[2]))

DPLL.gpd <- function(theta, data){
  library(evd)
  loc <- theta[1]
  scale <- theta[2]
  shape <-  theta[3]
#  R <- theta[4]
  ll <- data$att* log(pgpd(0.2*data$ave, loc,scale,shape)) + (5 - data$att) * log((1-pgpd(0.2*data$ave,loc,scale,shape)))
  return(-1*sum(ll))
}

DPLL.gpd(c(.5,1,0), dat.dp)
nlminb(c(-.77,5,.1), DPLL.gpd, data = dat.dp, control = list(rel.tol = 1e-14, iter.max = 5000, trace = 1))
res.gpd <- optim(c(0,1,0), DPLL.gpd, data = dat.dp, hessian = TRUE, method = 'L-BFGS-B', lower = c(-Inf,0,-Inf), control = list(maxit=500))
res.gpd
## print results
myopt <- function(nll, start, data, ...){
  res <- optim(start, nll, data = data, ...)
  res$nobs <- dim(data)[1]
  return(res)
}



summary.myopt <- function(ob) {
  vcov <- solve(ob$hessian)
  se <- sqrt(diag(vcov))
  z <- ob$par/se
  df <- data.frame(par = ob$par, se = se, z = z)
  AIC <- 2*ob$value + 2*length(ob$par)
  BIC <- 2*ob$value + log(ob$nobs) * length(ob$par)
  print(df)
  cat('\n LL:', -ob$value,'\n')
  cat('\n AIC:', AIC,'\n')
  cat('\n BIC:', BIC,'\n')
  print(vcov)
  
}

res.gpd <- myopt(DPLL.gpd, c(0,1,0), data = dat.dp, hessian = TRUE, method = 'L-BFGS-B', lower = c(-Inf,0,-Inf), control = list(maxit=500))
res.gpd

res.gpd.fix <- myopt(DPLL.gpd.fix, c(.8), data = dat.dp, hessian = TRUE, method = 'L-BFGS-B', lower = c(0.00001), upper = c(.999999),
                     control = list(maxit=500))
res.gpd.fix




summary.myopt(res.gpd)
summary.myopt(res.gpd.fix)

res.gpd.m <- myopt(DPLL.gpd, c(0,1,0), data = filter(dat.dp, Gender == 'M'), hessian = TRUE, method = 'L-BFGS-B', lower = c(-Inf,0,-Inf), control = list(maxit=500))
res.gpd.m
summary.myopt(res.gpd.m)


res.gpd.f <- myopt(DPLL.gpd, c(0,1,0), data = filter(dat.dp, Gender == 'F'), hessian = TRUE, method = 'L-BFGS-B', lower = c(-Inf,0,-Inf), control = list(maxit=500))
res.gpd.f
summary.myopt(res.gpd.f)

res.gpd.fix.m <- myopt(DPLL.gpd.fix, c(.8), data = filter(dat.dp, Gender == 'M'), hessian = TRUE, method = 'L-BFGS-B', lower = c(0.00001), upper = c(.999999),
                       control = list(maxit=500))
summary.myopt(res.gpd.fix.m)
res.gpd.fix.f <- myopt(DPLL.gpd.fix, c(.8), data = filter(dat.dp, Gender == 'F'), hessian = TRUE, method = 'L-BFGS-B', lower = c(0.00001), upper = c(.999999),
                       control = list(maxit=500))
summary.myopt(res.gpd.fix.f)



lnpar <- res$par
DPLL.gpd(lnpar,dat.dp)
exp(res$par)
lnhess <- res$hessian
ses <- diag(solve(lnhess))

plot(x <- seq(0,20,by=.1), pgpd(x, res.gpd.f$par[1],res.gpd.f$par[2],res.gpd.f$par[3]),  type = 'l', ylim = c(0,1))
lines(x <- seq(0,20,by=.1), pgpd(x, res.gpd.m$par[1],res.gpd.m$par[2],res.gpd.m$par[3]),  col = 'lightblue', type = 'l')

plot(x <- seq(0,20,by=.1), dgpd(x, res.gpd.f$par[1],res.gpd.f$par[2],res.gpd.f$par[3]),  type = 'l', ylim = c(0,0.28))
lines(x <- seq(0,20,by=.1), dgpd(x, res.gpd.m$par[1],res.gpd.m$par[2],res.gpd.m$par[3]),  col = 'lightblue', type = 'l')

plot(x <- seq(0,20,by=.1), dgpd(x, lnpar[1],lnpar[2], lnpar[3]),  type = 'l')
sqrt(diag(lnhess))
summary(rnorm(10000,lnpar[1], lnpar[2]))



##############
library(bbmle)

DPLL.exp.bb <- function(lambda){
  ll <- att* log(pexp(0.2*ave, lambda)) + (5 - att) * log((1-pexp(0.2*ave,lambda)))
  return(-1*sum(ll))
}

att <- dat.dp$att; ave <- dat.dp$ave

DPLL.exp.bb(.2)
## with a fixed mark invariant trigger
spar <- list(lambda = .7)
mle.exp <- mle2(DPLL.exp.bb, start = spar,  method = 'BFGS', control = list(maxit=5000))
summary(mle.exp)
p0 <- profile(mle.exp, trace = TRUE, try_harder = TRUE)
confint(p0)
confint(p0,method="quad")
confint(p0,method="uniroot")
plot(p0,plot.confstr=TRUE)
AIC(mle.exp)
BIC(mle.exp) #npt works
## BIC
AIC(mle.exp, k = log(length(att)))
logLik(mle.exp)
length(att)

## for pareto
## http://stats.stackexchange.com/questions/78168/how-to-know-if-my-data-fits-pareto-distribution
ppareto <- function(q, xm, alpha) ifelse(q > xm , 1 - (xm/q)**alpha, 0 )
DPLL.pareto.bb <- function(xm, alpha){
  ll <- att* log(ppareto(0.2*ave, xm, alpha)) + (5 - att) * log((1-ppareto(0.2*ave,xm,alpha)))
  return(-1*sum(ll))
}

att <- dat.dp$att; ave <- dat.dp$ave

DPLL.pareto.bb(.00001,.2)
## with a fixed mark invariant trigger
spar <- list(xm = 0.001, alpha = .2)
mle.pareto <- mle2(DPLL.pareto.bb, start = spar,  method = 'L-BFGS-B', lower = c(0.00001,-Inf) ,control = list(maxit=5000))

summary(mle.pareto)
p0 <- profile(mle.pareto, trace = TRUE, try_harder = TRUE)
confint(p0)
confint(p0,method="quad")
confint(p0,method="uniroot")
plot(p0,plot.confstr=TRUE)
AIC(mle.pareto)


DPLL.gpd.bb <- function(loc, scale, shape){
  library(evd)
#  loc <- theta[1]
#  scale <- theta[2]
#  shape <-  theta[3]
#  R <- theta[4]
  ll <- att* log(pgpd(0.2*ave, loc,scale,shape)) + (5 - att) * log((1-pgpd(0.2*ave,loc,scale,shape)))
  return(-1*sum(ll))
}

att <- dat.dp$att
ave <- dat.dp$ave
DPLL.gpd.bb(-0.7717872,  3.4,  0.1136894)
spar <- list(loc = -.7, scale = 3 , shape = 0.1)
rmle <- mle2(DPLL.gpd.bb, start = spar,  method = 'L-BFGS-B', 
             lower = c(loc = -Inf,scale = 0,shape =-Inf), control = list(maxit=5000))
summary(rmle)
p0 <- profile(rmle, maxsteps = 5000, trace = FALSE, try_harder = TRUE)
confint(p0)
confint(p0,method="quad")
confint(p0,method="uniroot")
plot(p0,plot.confstr=TRUE)
AIC(rmle, mle.exp)
AIC(rmle, k = log(length(att)))
AIC(rmle, mle.exp, k = log(length(att)))
logLik(rmle)

ppexp <- pexp(0.2*dat.dp$ave,coef(mle.exp))^dat.dp$att 
ppgpd <- pgpd(0.2*dat.dp$ave, coef(rmle)[1], coef(rmle)[2], coef(rmle)[3] )^dat.dp$att 
summary(ppexp)
summary(ppgpd)
N <- length(ppexp)
wsq <- sum(log(ppexp/ppgpd)^2) /N - (sum(log(ppexp/ppgpd))/N)^2
wsq
w <- sqrt(wsq)
w

LR <- sum(log(ppexp/ppgpd)) - 0.5*(1-3)*log(N)
V <- LR/(sqrt(N) * w)
V

 
abline(coef=c(0,1))
predict(rmle)





res <- optim(c(0,1,0,80), DPLL.fm.gpd, data = dat.dp, hessian = TRUE, method = 'L-BFGS-B', lower = c(-Inf,0,-Inf,-Inf), 
             control = list(maxit=5000))
res
lnpar <- res$par
exp(res$par)
lnhess <- res$hessian
ses <- diag(solve(lnhess))
ses
plot(x <- seq(0,20,by=.1), pgpd(x, lnpar[1],lnpar[2], lnpar[3]),  type = 'l')
plot(x <- seq(0,20,by=.1), dgpd(x, lnpar[1],lnpar[2], lnpar[3]),  type = 'l')
sqrt(diag(lnhess))
summary(rnorm(10000,lnpar[1], lnpar[2]))






DPLL.exp <- function(theta, data){
  m <- theta[1]
  ll <- data$att* log(pexp(0.2*data$ave, m)) + (5 - data$att) * log((1-pexp(0.2*data$ave,m)))
  return(-1*sum(ll))
}
DPLL.exp(.2,dat.dp)

res <- optim(1, DPLL.exp, data = dat.dp, hessian = TRUE, method = 'BFGS', control = list(maxit=500))
res
lnpar <- res$par
exp(res$par)
lnhess <- res$hessian
plot(x <- seq(-20,50,by=.1), pnorm(x, lnpar[1],lnpar[2]),  type = 'l')
sqrt(diag(lnhess))
summary(rnorm(10000,lnpar[1], lnpar[2]))




DPLL.norm <- function(theta, data){
  m <- theta[1]
  sd <- theta[2]
  ll <- data$att* log(pnorm(0.2*data$ave,m,sd)) + (5 - data$att) * log((1-pnorm(0.2*data$ave,m,sd)))
  return(-1*sum(ll))
}



DPLL.norm(c(m,sd),dat.dp)
res <- optim(c(m,sd), DPLL.norm, data = dat.dp, hessian = TRUE, method = 'BFGS', control = list(maxit=500))
res
lnpar <- res$par
exp(res$par)
lnhess <- res$hessian
plot(x <- seq(-20,50,by=.1), pnorm(x, lnpar[1],lnpar[2]),  type = 'l')
sqrt(diag(lnhess))
summary(rnorm(10000,lnpar[1], lnpar[2]))



## median
exp(lnpar[1])

## fit lognor
DPLL.ln <- function(theta, data){
  m <- theta[1]
  sd <- theta[2]
  ll <- data$att* log(plnorm(0.2*data$ave,m,sd)) + (5 - data$att) * log((1-plnorm(0.2*data$ave,m,sd)))
  return(-1*sum(ll))
}

DPLL.ln(c(m,sd),dat.dp)
res <- optim(c(m,sd), DPLL.ln, data = dat.dp, hessian = TRUE, method = 'BFGS', control = list(maxit=500))
res
lnpar <- res$par
exp(res$par)
lnhess <- res$hessian
plot(x <- seq(0,2,by=.1), plnorm(x, lnpar[1],lnpar[2]),  type = 'l')
sqrt(diag(lnhess))


summary(rlnorm(1000,lnpar[1], lnpar[2]))

## median
exp(lnpar[1])


library(bbmle)

## needs a new func

DPLL2 <- function(theta, natt, save){
  m <- theta[1]
  sd <- theta[2]
  ll <- natt* log(pnorm(0.2*save,m,sd)) + (5 - natt) * log((1-pnorm(0.2*save,m,sd)))
  return(-1*sum(ll))
}


mle2(DPLL2, start = list(m,sd), data = list( natt = dat.dp$att, save = dat.dp$ave), vecpar = TRUE, parnames = list('m','sd'))
## ordered logits for latent cutoffs, treat as descriptive
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)

dat.opeg <- read.dta("http://www.ats.ucla.edu/stat/data/ologit.dta")
head(dat.opeg)
dat.op <- mutate(dat, att = factor(att)) %>% filter(!is.na(ave), ave > 0)
ol.out <- polr(att ~ log(ave), data = dat.op, Hess=TRUE)
summary(ol.out)
str(ol.out)
exp(coef(ol.out))
### refromat for mlogit
library(mlogit)
library(tidyr)
rdat <- dat[rep(seq(1,dim(dat)[1]), each =5),]
rdat <- cbind(rdat, alt = rep(1:5,dim(dat)[1]))
head(rdat)
rdat <- rdat %>%   mutate(Gender = factor(Gender), mark = ave*alt/5, logm = log(mark), sqrtm = sqrt(mark), sqm = mark^2, cubm = mark^3, tcost = alt, choice = ifelse(att == alt,TRUE,FALSE)) %>% 
  dplyr::select(-missed, -ave, -att) %>% dplyr::filter(!is.na(mark), mark > 0) 


rd <- mlogit.data(rdat, choice = 'choice', shape = 'long', alt.var = 'alt')   
ml.out <- mlogit(choice ~  logm | 0,  rpar = c(logm = 'ln'),  print.level = 2, halton = NULL, data = rd, R=1000)
summary(ml.out)
AIC(ml.out)
rpar(ml.out)
fitted(ml.out)
fitted(ml.out, outcome = FALSE)
?fitted
apply(fitted(ml.out, outcome=FALSE), 2, mean)
## standard mlogit

ml.lin.out <- mlogit(choice ~   mark | 0, print.level = 2, halton = NA,rpar = c(mark = 'n'), data = rd, R=100)
summary(ml.lin.out)


## mpn
mnp.out <- mlogit(choice ~  logm | 0 ,  probit = TRUE,  print.level = 2, data = rd, R=100)
## doesn't work
mnp.out <- mlogit(choice ~  mark | 0 ,  probit = TRUE,  print.level = 2, data = rd, R=100)
#doesn't work

summary(mnp.out)
AIC(ml.out)
rpar(ml.out)
fitted(ml.out)
fitted(ml.out, outcome = FALSE)
?fitted
apply(fitted(ml.out, outcome=FALSE), 2, mean)



head(model.matrix(choice ~  tcost | 0 | mark,  data = rd))



d.ml <- mutate(dat, ben0 = 0,  ben1 = ave*1/5, ben2 = ave*2/5, ben3 = ave*3/5, ben4 = ave*4/5,  ben5 = ave,
        c0 = 0, c1 = 1, c2 = 2, c3 = 3, c4 = 4, c5 = 5)
sd.ml <- d.ml %>%  dplyr::select(-missed, -ave) %>% dplyr::select( ID = as.numeric(ID), choice = att, ben0,  ben1,ben2, ben3,ben4,ben5, c0, c1,c2,c3,c4, c5)
#sd.ml %>% dplyr::select( ID, att,  ben1,ben2, ben3,ben4,ben5, c1,c2,c3,c4, c5)
#gather(sd.ml, attrib,val, ben1, ben2, ben3, ben4, ben5, c1,c2,c3,c4,c5) %>% arrange(ID)
mlogit.data(as.data.frame(sd.ml), shape = 'wide', id.var = 'ID', varying = list(3:14), sep = '', choice = 'choice')



with(dat,table(cut(ave,breaks=10),missed))

ggplot(dat, aes(x = factor(missed), y= ave)) + geom_boxplot()
ggplot(dat, aes(x = factor(att), y= ave)) + geom_boxplot()

ggplot(dat, aes(x = ave, y = missed)) + geom_point() +geom_smooth(method = 'lm')

## SAC for SD

ldat <- split(dat,dat$missed)
sapply(ldat,function(x) sd(x$ave, na.rm = TRUE))

## Look at linerar model, poisson model and neg binomial
lm.out <- lm(missed ~ ave ,data=dat)
summary(lm.out)
coef(lm.out)
## implictly the estimates may describe a quadratic
## for the effort cost function
### coeff of squared effort (n)
b <- -5/(2*coef(lm.out)[2])
a <- 2*b*(coef(lm.out)[1] - 5)


#m1 <- glm(missed ~ av,data=dat, family = 'poisson')

dat.df <- filter(dat, is.finite(ave) & ave > 0)

m1 <- glm(att ~ log(ave),data=dat.df, na.action = na.exclude, control=glm.control(trace=1,maxit=100), family = 'poisson')
summary(m1)
## gamma is the cost elasticity of attendance
gamm <- (1+coef(m1)[2])/coef(m1)[2]
## K is the multiplier for the effort cost
K <- exp(1/(gamm-1) - coef(m1)[1])/(5*gamm)
x<-seq(0,5,by=.1)

### recover cost function - units are marks given up
cost<- function(x,k,gamma) {k*x^gamma}
cost.df <- data.frame(n = seq(0,5,by=.1), cost = cost(x,k=K,gamma=gamm))

ggplot(cost.df, aes(x = n, y = cost)) + geom_line()


cost(5,K,gamm) - cost(4,K,gamm)
.3*.2*(cost(5,K,gamm) - cost(4,K,gamm))/12


## av = 20
lines(x,x*20/5)
lines(x,x*50/5)
## optimal n is beta1 * av^beta2
lnstar <- coef(m1)[1] + coef(m1)[2] * log(dat.df$ave)
nstar <- exp(lnstar)
plot(nstar,dat.df$att)


library(sandwich)
cov.m1<-vcovHC (m1, type="HC0")
std.err<-sqrt(diag(cov.m1))
r.est<-cbind(estimate=m1$coefficients, std.err,
       pvalue=round(2*(1-pnorm(abs(m1$coefficients/std.err))),5),
       lower=m1$coefficients-1.96*std.err,
       upper=m1$coefficients+1.96*std.err)

r.est
gof<-cbind(res.diviance=m1$deviance, df=m1$df.residual,
           p=1-pchisq(m1$deviance, m1$df.residual))
gof

library(msm)
estmean<-coef(m1)
var<-cov.m1
s2<-deltamethod (~ exp(x2), estmean, var)

irr<-cbind(irr=exp(estmean)[2], r.serr=s2,
           lower=exp(m1$coefficients-1.96*std.err)[2],
           upper=exp(m1$coefficients+1.96*std.err)[2])

irr

df <- data.frame(ave = seq(0,100,by=5))
 predict(m1,df,type = 'response')
plot(df$ave,predict(m1,df,type = 'response'))


######
m1.pois <- m1
library(car)
m1 <- glm.nb(missed ~ ave,data=dat.df)
summary(m1)
library(lmtest)

lrtest(m1.pois,m1)

exp(coef(m1))-1


df <- data.frame(ave = seq(0,100,by=5))
 predict(m1,df,type = 'response')
plot(df$ave,predict(m1,df,type = 'response'))

