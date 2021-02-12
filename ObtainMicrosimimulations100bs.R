# run individual adl/iadl using msm

#install.packages("msm")
#install.packages("haven")

library(msm)
library(haven)

#define function msm_dif that includes all the steps to generate microsimulations for each ADL/IADL

msm_dif<-function (dif=''){
  
  
  dat$agem50 <- dat$ageint-50 # start age at 0 not 50
  dat[,"dif3"]<-dat[,paste(dif)]+1 # make states 1/2/3 not 0/1/2
  dat <- dat[!is.na(dat$dif3),] # drop missing ones
  dat <- dat[!is.na(dat$RAGENDER),] # drop missing ones
  dat$ragender <- dat$RAGENDER
  
  Q <- rbind(c(0,0.25,0.25), c(0.25,0,0.25), c(0,0,0)) # setup initial instantaneous transition matrix
  dat$agem50[dat$dif3==3] <- dat$agem50[dat$dif3==3] + 0.1
  Q.crude  <- crudeinits.msm(dif3 ~ agem50, newid2, data=dat, qmatrix=Q)
  
  # get the initial value of the -2loglik
  dif.loglik <- msm(dif3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                    qmatrix = Q.crude,  opt.method="optim", control = list(trace=1, REPORT=1),fixedpars=TRUE)
  loglik<-dif.loglik$minus2loglik 
  
  # run model with age as time, age and gender as covariates, and initial value of the -2loglik
  dif.msm <- msm(dif3 ~ agem50, subject=newid2, data = dat, covariates=~agem50+ragender,
                 qmatrix = Q.crude,  opt.method="optim", control = list(fnscale = loglik))
  
  #### get the 1 year probability transition matrices for each age and gender
  pmatoutdifG1 <- vector("list",length=51)
  pmatoutdifG2 <- vector("list",length=51)
  
  for (i in 1:51){
    pmatoutdifG1[[i]] <- pmatrix.msm(dif.msm, t=1, covariates=list(agem50=i-1,ragender=1))
    pmatoutdifG2[[i]] <- pmatrix.msm(dif.msm, t=1, covariates=list(agem50=i-1,ragender=2))
    
  }
  
  # put these prob trans matrices into form used by SPACE macro so can use program for microsim based on SPACE output
  ptdat <- matrix(0,51*2*2,6) # one line per age per gender per nondead state
  dimnames(ptdat) <- list(NULL,c("ragender","BEGST","AGE","P1","P2","P3"))
  ptdat[,1] <- rep(c(1,2),each=51*2)
  ptdat[,2] <- rep(c(1,2),51*2)
  ptdat[,3] <- rep(rep(50:100,each=2),2)
  for (i in 1:51){
    ibot <- i+51
    ptdat[(2*i-1):(2*i),4:6] <- pmatoutdifG1[[i]][1:2,]
    ptdat[(2*ibot-1):(2*ibot),4:6] <- pmatoutdifG2[[i]][1:2,]
    
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
  
  # now run all the ADL/IADL for both gender 1 and 2 with 50000 people each sim
  out.dif.g1 <- runonedataset(ptdat[ptdat$ragender==1,],nsim=50000)
  out.dif.g2 <- runonedataset(ptdat[ptdat$ragender==2,],nsim=50000)
  
  # add gender and replicate variables
  out.dif.g1<-as.data.frame(out.dif.g1)
  out.dif.g1$replicate<-b
  out.dif.g1$gender<-1
  
  out.dif.g2<-as.data.frame(out.dif.g2)
  out.dif.g2$replicate<-b
  out.dif.g2$gender<-2
  
  # stack data sets for both genders
  out.dif<-rbind(out.dif.g1, out.dif.g2)
  
  # save data
  write.csv(out.dif, file = paste(dif, b, '.csv', sep=""),row.names=FALSE)
  
  rm(list=setdiff(ls(), c("dset", "msm_dif")))  #remove all objects except "dset" and "msm_dif" so they can be used for next ADLs/IADLs (dset) and replicates (msm_diff)
  
}

for (i in 1:100) {
  
  #read SAS data
  setwd("path")
  dset <- read_sas(paste('bs', i,".sas7bdat",sep="")) # dset will remain in memory for all the ADLS/IADLs
  dset <- dset[order(dset$newid2, dset$wavenum),]
  setwd("path")
  dat<-as.data.frame(dset) # copy dset to dat so dat can be modified for each ADL/IADL
  b<-i #save bootstrap sample number. I will use it later
  
  #call msm_diff function for each ADL/IADL
  msm_dif(dif='bathdif')
  msm_dif(dif='beddif')
  msm_dif(dif='dressdif')
  msm_dif(dif='eatdif')
  msm_dif(dif='toiletdif')
  msm_dif(dif='walkdif')
  msm_dif(dif='phonedif')
  msm_dif(dif='moneydif')
  msm_dif(dif='meddif')
  msm_dif(dif='shopdif')
  msm_dif(dif='mealdif')
  
  rm(list=setdiff(ls(), "msm_dif"))  #remove all objects except "msm_diff" function so it can be used for rest of the replicates
  
  
}



