
library(msm)
library(Hmisc)

##################################### bath ##################################### 

packageDescription("msm", fields = "Version")
# ver "1.6.8" (16 December, 2019)
R.version$version.string
# "R version 4.0.4 (2021-02-15)"

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$bath3 <- dat$bathdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$bath3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(bath3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$bath3==3] <- dat$agem50[dat$bath3==3] + 0.1
dat$year[dat$bath3==3] <- dat$year[dat$bath3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(bath3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
bath.msm <- msm(bath3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
	qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines and gender: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(bath3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                 qmatrix = Q.crude,  fixedpars=TRUE)
loglik<-dif.loglik$minus2loglik 
bath.msm2 <- msm(bath3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                 qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(bath.msm, bath.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters
# -2 log LR df p
# bath.msm2  586.63 12 0

### Predictions
pmatrix.msm(bath.msm, t = 9)
pmatrix.msm(bath.msm2, t = 9)

pmatrix.msm(bath.msm, t = 20)
pmatrix.msm(bath.msm2, t = 20)

pnext.msm(bath.msm)
pnext.msm(bath.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(bath.msm , main = "Bath: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(bath.msm2 , main = "Bath: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))


##################################### bed ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$bed3 <- dat$beddif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$bed3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(bed3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$bed3==3] <- dat$agem50[dat$bed3==3] + 0.1 
dat$year[dat$bed3==3] <- dat$year[dat$bed3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(bed3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
bed.msm <- msm(bed3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
               qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(bed3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE) 
loglik<-dif.loglik$minus2loglik 
bed.msm2 <- msm(bed3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(bed.msm, bed.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(bed.msm, t = 9)
pmatrix.msm(bed.msm2, t = 9)

pmatrix.msm(bed.msm, t = 20)
pmatrix.msm(bed.msm2, t = 20)

pnext.msm(bed.msm)
pnext.msm(bed.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(bed.msm , main = "Bed: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(bed.msm2 , main = "Bed: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))


##################################### dress ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$dress3 <- dat$dressdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$dress3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(dress3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$dress3==3] <- dat$agem50[dat$dress3==3] + 0.1 
dat$year[dat$dress3==3] <- dat$year[dat$dress3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(dress3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
dress.msm <- msm(dress3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
                 qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(dress3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE)   
loglik<-dif.loglik$minus2loglik 
dress.msm2 <- msm(dress3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(dress.msm, dress.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(dress.msm, t = 9)
pmatrix.msm(dress.msm2, t = 9)

pmatrix.msm(dress.msm, t = 20)
pmatrix.msm(dress.msm2, t = 20)

pnext.msm(dress.msm)
pnext.msm(dress.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(dress.msm , main = "Dress: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(dress.msm2 , main = "Dress: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))

##################################### eat ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$eat3 <- dat$eatdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$eat3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(eat3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$eat3==3] <- dat$agem50[dat$eat3==3] + 0.1 
dat$year[dat$eat3==3] <- dat$year[dat$eat3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(eat3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
eat.msm <- msm(eat3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
               qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(eat3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE) 
loglik<-dif.loglik$minus2loglik 
eat.msm2 <- msm(eat3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(eat.msm, eat.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(eat.msm, t = 9)
pmatrix.msm(eat.msm2, t = 9)

pmatrix.msm(eat.msm, t = 20)
pmatrix.msm(eat.msm2, t = 20)

pnext.msm(eat.msm)
pnext.msm(eat.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(eat.msm , main = "Eat: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(eat.msm2 , main = "Eat: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))


##################################### toilet ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$toilet3 <- dat$toiletdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$toilet3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(toilet3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$toilet3==3] <- dat$agem50[dat$toilet3==3] + 0.1 
dat$year[dat$toilet3==3] <- dat$year[dat$toilet3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(toilet3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
toilet.msm <- msm(toilet3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
                  qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))


# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(toilet3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE)  
loglik<-dif.loglik$minus2loglik 
toilet.msm2 <- msm(toilet3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                   qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(toilet.msm, toilet.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(toilet.msm, t = 9)
pmatrix.msm(toilet.msm2, t = 9)

pmatrix.msm(toilet.msm, t = 20)
pmatrix.msm(toilet.msm2, t = 20)

pnext.msm(toilet.msm)
pnext.msm(toilet.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(toilet.msm , main = "Toilet: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(toilet.msm2 , main = "Toilet: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))

##################################### walk ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$walk3 <- dat$walkdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$walk3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(walk3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$walk3==3] <- dat$agem50[dat$walk3==3] + 0.1 
dat$year[dat$walk3==3] <- dat$year[dat$walk3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(walk3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
walk.msm <- msm(walk3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
                qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(walk3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE)  
loglik<-dif.loglik$minus2loglik 
walk.msm2 <- msm(walk3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                 qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(walk.msm, walk.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(walk.msm, t = 9)
pmatrix.msm(walk.msm2, t = 9)

pmatrix.msm(walk.msm, t = 20)
pmatrix.msm(walk.msm2, t = 20)

pnext.msm(walk.msm)
pnext.msm(walk.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(walk.msm , main = "Walk: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(walk.msm2 , main = "Walk: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))


##################################### phone ##################################### 

### once data create in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$phone3 <- dat$phonedif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$phone3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(phone3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$phone3==3] <- dat$agem50[dat$phone3==3] + 0.1 
dat$year[dat$phone3==3] <- dat$year[dat$phone3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(phone3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
phone.msm <- msm(phone3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
                 qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))
  
# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(phone3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE)  
loglik<-dif.loglik$minus2loglik 
phone.msm2 <- msm(phone3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(phone.msm, phone.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(phone.msm, t = 9)
pmatrix.msm(phone.msm2, t = 9)

pmatrix.msm(phone.msm, t = 20)
pmatrix.msm(phone.msm2, t = 20)

pnext.msm(phone.msm)
pnext.msm(phone.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(phone.msm , main = "Phone: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(phone.msm2 , main = "Phone: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))

##################################### money ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$money3 <- dat$moneydif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$money3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(money3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$money3==3] <- dat$agem50[dat$money3==3] + 0.1 
dat$year[dat$money3==3] <- dat$year[dat$money3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(money3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
money.msm <- msm(money3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
                 qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(money3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE) 
loglik<-dif.loglik$minus2loglik 
money.msm2 <- msm(money3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(money.msm, money.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(money.msm, t = 9)
pmatrix.msm(money.msm2, t = 9)

pmatrix.msm(money.msm, t = 20)
pmatrix.msm(money.msm2, t = 20)

pnext.msm(money.msm)
pnext.msm(money.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(money.msm , main = "Money: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(money.msm2 , main = "Money: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))

##################################### med ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$med3 <- dat$meddif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$med3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(med3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$med3==3] <- dat$agem50[dat$med3==3] + 0.1 
dat$year[dat$med3==3] <- dat$year[dat$med3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(med3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
med.msm <- msm(med3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
               qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(med3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE) 
loglik<-dif.loglik$minus2loglik 
med.msm2 <- msm(med3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(med.msm, med.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(med.msm, t = 9)
pmatrix.msm(med.msm2, t = 9)

pmatrix.msm(med.msm, t = 20)
pmatrix.msm(med.msm2, t = 20)

pnext.msm(med.msm)
pnext.msm(med.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(med.msm , main = "Med: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(med.msm2 , main = "Med: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))


##################################### shop ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$shop3 <- dat$shopdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$shop3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(shop3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$shop3==3] <- dat$agem50[dat$shop3==3] + 0.1 
dat$year[dat$shop3==3] <- dat$year[dat$shop3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(shop3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
shop.msm <- msm(shop3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
                qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(shop3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE) 
loglik<-dif.loglik$minus2loglik 
shop.msm2 <- msm(shop3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                 qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(shop.msm, shop.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(shop.msm, t = 9)
pmatrix.msm(shop.msm2, t = 9)

pmatrix.msm(shop.msm, t = 20)
pmatrix.msm(shop.msm2, t = 20)

pnext.msm(shop.msm)
pnext.msm(shop.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(shop.msm , main = "Shop: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(shop.msm2 , main = "Shop: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))

##################################### meal ##################################### 

### once data created in R can just read it in quickly
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$meal3 <- dat$mealdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$meal3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

dat$year <- dat$ageint-dat$ageintbl #create year variable

statetable.msm(meal3, newid, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix

### need to have death time later than last obs time
dat$agem50[dat$meal3==3] <- dat$agem50[dat$meal3==3] + 0.1 
dat$year[dat$meal3==3] <- dat$year[dat$meal3==3] + 0.1 

rownames(Q) <- colnames(Q) <- c("No Difficulty", "Difficulty", "Death")

Q.crude  <- crudeinits.msm(meal3 ~ year, newid, data=dat, qmatrix=Q)

# run model with year as time and age linear and gender as covariates
meal.msm <- msm(meal3 ~ year, subject=newid, data = dat, covariates=~ageint+ragender,
                qmatrix = Q.crude,  opt.method="fisher", control = list(trace=1, REPORT=1))
      

# run model using restricted cubic splines: https://rdrr.io/cran/Hmisc/man/rcspline.eval.html
splineVals<-rcspline.eval(dat$ageint, inclx=TRUE) 
colnames(splineVals)<-c('agelin', 'sp1', 'sp2', 'sp3')
# If knots location not given, knots will be estimated using default quantiles of x
# number of knots. Default is 5. The minimum value is 3.
# inclx=TRUE: to add x as the first column of the returned matrix
# number of splines=number of knots -1= 5-1=4splines (this includes linear time as the 1st column)

attributes(splineVals)
# $dim
#[1] 232174      4

#$knots
#[1]  53 60 67 76 88

class(splineVals)

dat$agelin<-splineVals[,'agelin']
colnames(dat)
summary(splineVals[,'agelin'])
summary(dat$agelin)

dat$sp1<-splineVals[,'sp1']
dat$sp2<-splineVals[,'sp2']
dat$sp3<-splineVals[,'sp3']

dif.loglik <- msm(meal3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                  qmatrix = Q.crude,  fixedpars=TRUE)  
loglik<-dif.loglik$minus2loglik 
meal.msm2 <- msm(meal3 ~ year, subject=newid, data = dat, covariates=~agelin+sp1+sp2+sp3+ragender,
                 qmatrix = Q.crude,  control = list(fnscale = loglik, maxit = 10000))

### Model comparison
options(digits=5)
lrtest.msm(meal.msm, meal.msm2) #Two or more fitted multi-state models, as returned by msm, ordered by increasing numbers of parameters

### Predictions
pmatrix.msm(meal.msm, t = 9)
pmatrix.msm(meal.msm2, t = 9)

pmatrix.msm(meal.msm, t = 20)
pmatrix.msm(meal.msm2, t = 20)

pnext.msm(meal.msm)
pnext.msm(meal.msm2)

# Survival plots
par(mfrow = c(1, 2))
plot(meal.msm , main = "Meal: Age linear model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
plot(meal.msm2 , main = "Meal: Age splines model", xlab = "Years since baseline interview", lwd = 1, legend.pos=c(10,10)) 
legend(10, 1 ,c("From No difficulty", "From Difficulty"), text.col=c("red", "turquoise1" ), bty = "n", cex=0.8)
par(mfrow = c(1, 1))

