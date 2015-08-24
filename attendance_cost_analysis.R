setwd("~/Desktop/IE/IELabs/Marks")
d <- readRDS('marks15_final.rds')
g <- readRDS('~/Desktop/IE/IELabs/Groups/sessdb.rds')
s <- readRDS('~/Desktop/IE/IELabs/sall.rds')

library(tidyr)
library(dplyr)
library(MASS)

labcols <- paste('Lab.',1:5,sep = '')
## make long form
dl <- d %>% select(ID, starts_with('Lab') ) %>% gather(session, mark, Lab.1:Lab.5) #%>% head
dl %>% group_by(session) %>% summarise(ave = mean(mark, na.rm = TRUE)) %>%   ggplot(aes(x = session, y = ave, group = 1)) + geom_line()
dl %>% ggplot(aes(x = session, y = mark)) + geom_boxplot(fill = 'lightblue') + geom_violin(alpha = 0.5)

## ID group and week
gmem <- gather(s,  wk, grp, gs1:gs5) %>% select(ID, Course, title, grp,wk) %>% mutate(session = gsub('gs','',wk)) %>% select(-wk)

## marks by week including group attended
mbywk <- merge(gmem,transmute(dl, ID, session = gsub('Lab.','',session),mark), by = c('ID','session'))

all <- g %>% select(grp, sat, session = weeks) %>% merge(mbywk,by = c('grp','session'))
all <- mutate(all, rmark = ifelse(is.na(mark),0,mark))
all <- all %>% mutate(Gender = ifelse(title == 'Mr', 'M','F')) %>% select(-title)
## tutor boxplot
ggplot(all, aes(x = sat, y = mark)) + geom_boxplot() + facet_grid( ~ session)

## group boxplot
ggplot(all, aes(x = factor(grp), y = mark)) + geom_boxplot() + facet_wrap( ~ session)

## week boxplot
ggplot(all, aes(x = factor(session), y = mark)) + geom_boxplot()

## Course boxplot
ggplot(all, aes(x = Course, y = mark)) + geom_boxplot()

## Gender boxplot
ggplot(all, aes(x = Gender, y = mark)) + geom_boxplot()

## regression
summary(lm(mark ~ factor(grp) + factor(session)  + factor(sat) + Gender + Course,data = all))


## Look at missed
dat <- all %>% group_by(ID) %>% summarise(missed = sum(is.na(mark)), ave = mean(mark, na.rm = TRUE))

dat$att <- 5 - dat$missed

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

