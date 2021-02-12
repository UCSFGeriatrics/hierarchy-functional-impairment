# run individual adl/iadl using msm

#install.packages("msm")
library(msm)

##################################### bath ##################################### 

setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$bath3 <- dat$bathdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$bath3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(bath3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
### need to have death time later than last obs time
dat$agem50[dat$bath3==3] <- dat$agem50[dat$bath3==3] + 0.1

Q.crude  <- crudeinits.msm(bath3 ~ agem50, newid2, data=dat, qmatrix=Q)

## We decided to use opt.method="optim" with control = list(fnscale = initial value of -2loglik) since minimized -2loglik was smaller

# Get the initial value of the -2loglik
bath.msm <- msm(bath3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3218480.915094 

# run model with age as time and age and gender as covariate
bath.msm <- msm(bath3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
	qmatrix = Q.crude,  opt.method="optim", control = list(fnscale = 3218480.915094))

#### get the 1 year probability transition matrices for each age and gender
pmatoutbathG1 <- vector("list",length=51)
pmatoutbathG2 <- vector("list",length=51)

for (i in 1:51){
pmatoutbathG1[[i]] <- pmatrix.msm(bath.msm, t=1, covariates=list(agem50=i-1,ragender=1))
pmatoutbathG2[[i]] <- pmatrix.msm(bath.msm, t=1, covariates=list(agem50=i-1,ragender=2))

}

# put these prob trans matrices into right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
	ibot <- i+51
	ptdat[(2*i-1):(2*i),4:6] <- pmatoutbathG1[[i]][1:2,]
	ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutbathG2[[i]][1:2,]

}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
	outsim <- matrix(1, nsim, maxage-minage+1) 
	for (yidx in 1:(ncol(outsim)-1)){
		tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
		for (sidx in 1:(nstate-1)){
			whichidx <- (1:nsim)[outsim[,yidx]==sidx]
			numwhich <- length(whichidx)
			if (numwhich>0){
				simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
				outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
			}
		}
		whichdead <- (1:nsim)[outsim[,yidx]==nstate]
		numdead <- length(whichdead)
		if (numdead){
			outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
		}
	}
	
	return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.bath.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.bath.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)


save("out.bath.g1","out.bath.g2", file="bathSims.Rdata")

##################################### bed ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$bed3 <- dat$beddif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$bed3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(bed3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q

### need to have death time later than last obs time
dat$agem50[dat$bed3==3] <- dat$agem50[dat$bed3==3] + 0.1  

Q.crude  <- crudeinits.msm(bed3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
bed.msm <- msm(bed3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
               qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3310882.127655 

# run model with age as time and age and gender as covariate
bed.msm <- msm(bed3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
               qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =3310882.127655))

#### get the 1 year probability transition matrices for each age and gender
pmatoutbedG1 <- vector("list",length=51)
pmatoutbedG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutbedG1[[i]] <- pmatrix.msm(bed.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutbedG2[[i]] <- pmatrix.msm(bed.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutbedG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutbedG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.bed.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.bed.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.bed.g1","out.bed.g2", file="bedSims.Rdata")


##################################### dress ##################################### 
rm(list = ls(all.names = TRUE)) 
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$dress3 <- dat$dressdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$dress3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(dress3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$dress3==3] <- dat$agem50[dat$dress3==3] + 0.1  

Q.crude  <- crudeinits.msm(dress3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
dress.msm <- msm(dress3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3646163.739229 

# run model with age as time and age and gender as covariate
dress.msm <- msm(dress3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =3646163.739229))

#### get the 1 year probability transition matrices for each age and gender
pmatoutdressG1 <- vector("list",length=51)
pmatoutdressG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutdressG1[[i]] <- pmatrix.msm(dress.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutdressG2[[i]] <- pmatrix.msm(dress.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutdressG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutdressG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.dress.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.dress.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.dress.g1","out.dress.g2", file="dressSims.Rdata")

##################################### eat ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$eat3 <- dat$eatdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$eat3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(eat3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$eat3==3] <- dat$agem50[dat$eat3==3] + 0.1  

Q.crude  <- crudeinits.msm(eat3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
eat.msm <- msm(eat3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
               qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))

#initial  value 2712235.033516 
 
# run model with age as time and age and gender as covariate
eat.msm <- msm(eat3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
               qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =2712235.033516))

#### get the 1 year probability transition matrices for each age and gender
pmatouteatG1 <- vector("list",length=51)
pmatouteatG2 <- vector("list",length=51)

for (i in 1:51){
  pmatouteatG1[[i]] <- pmatrix.msm(eat.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatouteatG2[[i]] <- pmatrix.msm(eat.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatouteatG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatouteatG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.eat.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.eat.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.eat.g1","out.eat.g2", file="eatSims.Rdata")


##################################### toilet ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$toilet3 <- dat$toiletdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$toilet3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(toilet3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$toilet3==3] <- dat$agem50[dat$toilet3==3] + 0.1  

Q.crude  <- crudeinits.msm(toilet3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
toilet.msm <- msm(toilet3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                  qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3082156.471446 

# run model with age as time and age and gender as covariate
toilet.msm <- msm(toilet3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                  qmatrix = Q.crude,  opt.method="optim", control = list(fnscale = 3082156.471446))


#### get the 1 year probability transition matrices for each age and gender
pmatouttoiletG1 <- vector("list",length=51)
pmatouttoiletG2 <- vector("list",length=51)

for (i in 1:51){
  pmatouttoiletG1[[i]] <- pmatrix.msm(toilet.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatouttoiletG2[[i]] <- pmatrix.msm(toilet.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatouttoiletG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatouttoiletG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.toilet.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.toilet.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.toilet.g1","out.toilet.g2", file="toiletSims.Rdata")


##################################### walk ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$walk3 <- dat$walkdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$walk3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(walk3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$walk3==3] <- dat$agem50[dat$walk3==3] + 0.1  

Q.crude  <- crudeinits.msm(walk3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
walk.msm <- msm(walk3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3257989.194068 

# run model with age as time and age and gender as covariate
walk.msm <- msm(walk3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =3257989.194068))


#### get the 1 year probability transition matrices for each age and gender
pmatoutwalkG1 <- vector("list",length=51)
pmatoutwalkG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutwalkG1[[i]] <- pmatrix.msm(walk.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutwalkG2[[i]] <- pmatrix.msm(walk.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutwalkG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutwalkG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.walk.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.walk.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.walk.g1","out.walk.g2", file="walkSims.Rdata")


##################################### phone ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$phone3 <- dat$phonedif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$phone3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(phone3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$phone3==3] <- dat$agem50[dat$phone3==3] + 0.1  

Q.crude  <- crudeinits.msm(phone3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
phone.msm <- msm(phone3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 2823842.259317 
 
# run model with age as time and age and gender as covariate
phone.msm <- msm(phone3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =2823842.259317))


#### get the 1 year probability transition matrices for each age and gender
pmatoutphoneG1 <- vector("list",length=51)
pmatoutphoneG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutphoneG1[[i]] <- pmatrix.msm(phone.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutphoneG2[[i]] <- pmatrix.msm(phone.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutphoneG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutphoneG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.phone.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.phone.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.phone.g1","out.phone.g2", file="phoneSims.Rdata")


##################################### money ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$money3 <- dat$moneydif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$money3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(money3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time

dat$agem50[dat$money3==3] <- dat$agem50[dat$money3==3] + 0.1  

Q.crude  <- crudeinits.msm(money3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
money.msm <- msm(money3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3125314.768834 

# run model with age as time and age and gender as covariate
money.msm <- msm(money3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =3125314.768834))

#### get the 1 year probability transition matrices for each age and gender
pmatoutmoneyG1 <- vector("list",length=51)
pmatoutmoneyG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutmoneyG1[[i]] <- pmatrix.msm(money.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutmoneyG2[[i]] <- pmatrix.msm(money.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutmoneyG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutmoneyG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.money.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.money.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.money.g1","out.money.g2", file="moneySims.Rdata")

##################################### med ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$med3 <- dat$meddif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$med3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(med3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$med3==3] <- dat$agem50[dat$med3==3] + 0.1  

Q.crude  <- crudeinits.msm(med3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
med.msm <- msm(med3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
               qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 2810084.794593 

# run model with age as time and age and gender as covariate
med.msm <- msm(med3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
               qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =2810084.794593))


#### get the 1 year probability transition matrices for each age and gender
pmatoutmedG1 <- vector("list",length=51)
pmatoutmedG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutmedG1[[i]] <- pmatrix.msm(med.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutmedG2[[i]] <- pmatrix.msm(med.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutmedG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutmedG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.med.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.med.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.med.g1","out.med.g2", file="medSims.Rdata")


##################################### shop ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$shop3 <- dat$shopdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$shop3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(shop3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$shop3==3] <- dat$agem50[dat$shop3==3] + 0.1  

Q.crude  <- crudeinits.msm(shop3 ~ agem50, newid2, data=dat, qmatrix=Q)

# run model with age as time and age and gender as covariate
shop.msm <- msm(shop3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3348953.938219 

# run model with age as time and age and gender as covariate
shop.msm <- msm(shop3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =3348953.938219 ))


#### get the 1 year probability transition matrices for each age and gender
pmatoutshopG1 <- vector("list",length=51)
pmatoutshopG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutshopG1[[i]] <- pmatrix.msm(shop.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutshopG2[[i]] <- pmatrix.msm(shop.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutshopG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutshopG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.shop.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.shop.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.shop.g1","out.shop.g2", file="shopSims.Rdata")


##################################### meal ##################################### 
rm(list = ls(all.names = TRUE))
library(msm)
setwd("path")
load("originaldata.Rdata")

table(dat$ageint)
dat$agem50 <- dat$ageint-50 # start age at 0 not 50
dat$meal3 <- dat$mealdif+1 # make states 1/2/3 not 0/1/2
dat <- dat[!is.na(dat$meal3),] # drop missing ones
dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
dat$ragender <- dat$RAGENDER

statetable.msm(meal3, newid2, data=dat) # check transition counts

Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
Q
### need to have death time later than last obs time
dat$agem50[dat$meal3==3] <- dat$agem50[dat$meal3==3] + 0.1  

Q.crude  <- crudeinits.msm(meal3 ~ agem50, newid2, data=dat, qmatrix=Q)

# Get the initial value of the -2loglik
meal.msm <- msm(meal3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1))
#initial  value 3002337.250029 

# run model with age as time and age and gender as covariate
meal.msm <- msm(meal3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                qmatrix = Q.crude,  opt.method="optim", control = list(fnscale =3002337.250029))


#### get the 1 year probability transition matrices for each age and gender
pmatoutmealG1 <- vector("list",length=51)
pmatoutmealG2 <- vector("list",length=51)

for (i in 1:51){
  pmatoutmealG1[[i]] <- pmatrix.msm(meal.msm, t=1, covariates=list(agem50=i-1,ragender=1))
  pmatoutmealG2[[i]] <- pmatrix.msm(meal.msm, t=1, covariates=list(agem50=i-1,ragender=2))
  
}

# put these prob trans matrices into the right format
ptdat <- matrix(0,51*2*2,6) # one line per age per gender per non dead state
dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
ptdat[,1] <- rep(c(1,2),each=51*2)
ptdat[,2] <- rep(c(1,2),51*2)
ptdat[,3] <- rep(rep(50:100,each=2),2)
for (i in 1:51){
  ibot <- i+51
  ptdat[(2*i-1):(2*i),4:6] <- pmatoutmealG1[[i]][1:2,]
  ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutmealG2[[i]][1:2,]
  
}
ptdat <- as.data.frame(ptdat)

#### simulate cohort
runonedataset <- function(dat, minage=50, maxage=100, nstate=3, nsim=20000, dcode=3){
  outsim <- matrix(1, nsim, maxage-minage+1) 
  for (yidx in 1:(ncol(outsim)-1)){
    tptmp <- dat[(dat$AGE==yidx+minage-1),4:(4+nstate-1)] # read current tp matrix
    for (sidx in 1:(nstate-1)){
      whichidx <- (1:nsim)[outsim[,yidx]==sidx]
      numwhich <- length(whichidx)
      if (numwhich>0){
        simtmp <- sample(1:nstate, numwhich, replace=TRUE, prob=tptmp[sidx,])
        outsim[whichidx,yidx+1] <- simtmp # put next age of simulated outcomes
      }
    }
    whichdead <- (1:nsim)[outsim[,yidx]==nstate]
    numdead <- length(whichdead)
    if (numdead){
      outsim[whichdead,yidx+1] <- dcode # once dead stay dead!
    }
  }
  
  return(outsim)
}

# now run all the ADL/IADL for both gender 1 and 2 with 1M people each sim

out.meal.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=1000000)
out.meal.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=1000000)

save("out.meal.g1","out.meal.g2", file="mealSims.Rdata")

rm(list = ls(all.names = TRUE))













