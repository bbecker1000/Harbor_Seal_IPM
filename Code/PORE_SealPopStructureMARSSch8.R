
##%######################################################%##
#                                                          #
####               modified for PORE data               ####
#                                                          #
##%######################################################%##

## data processed in PORE_SealTrend.....
## change dataset name so no conflicts between SealTrend and SealPopStructure
dat2 ## get from SealTrend Script before transposed

dat3 <- t(dat2)
###################################################
### code chunk number 8: Cs01_set.up.data
###################################################
#years = dat3[,1] #first col is years
#n = nrow(dat3)-1 # for reference later
#sealData = t(dat3[,c(2:6)]) # transpose and name sealData and remove years

n = nrow(dat3)-1
years = dat3[1,]     # remove years
sealData = dat3[2:nrow(dat3),]
legendnames = (unlist(dimnames(sealData)[1]))
#legendnames = legendnames[-1] ## get rid of year

###################################################
### code chunk number 9: Cs02_fig1
###################################################

## use ggplot to plot seal data


###################################################
### code chunk number 10: Cs03_set.up.Z.models
###################################################
#H1 Site
Z1=factor(c("BL", "DE", "DP", "TB", "TP")) #site
#H2 coastal+estuary
Z2=factor(c("E","E","C","E","C")) 
#H3 N and S
Z3=factor(c("S","S","S","N","N")) 
#H4 N and S and BL
Z4=factor(c("BL","S","S","N","N"))
#H5 panmictic
Z5=factor(rep("pan",5)) 
#H6 coastal and estuary + BL
Z6=factor(c("BL","E","C","E","C"))

Z.models=list(Z1,Z2,Z3,Z4,Z5,Z6)
names(Z.models)=
  c("Site","Coast+Est", "North-South","North/South/BL", "Panmictic", "Coast+Est+BL")


###################################################
### code chunk number 11: Cs04_Q.models
###################################################
Q.models=c("diagonal and equal", "diagonal and unequal")

###################################################
### code chunk number 12: Cs04a_other.models
###################################################
##These all worked well for SealTrend Section.
U.model="unequal"
R.model="diagonal and equal"
A.model="scaling"
B.model="identity" #"unconstrained" #"
x0.model="unequal"
V0.model="zero"
model.constant=list(
  U=U.model, R=R.model, A=A.model, B=B.model,
  x0=x0.model, V0=V0.model, tinitx=0)


###################################################
### code chunk number 13: Cs05_run.the.models  (runs in 75 sec with H=5)
###################################################
out.tab=NULL
fits=list()
for(i in 1:length(Z.models)){
  for(Q.model in Q.models){
    fit.model = c(list(Z=Z.models[[i]], Q=Q.model), model.constant)
    fit = MARSS(sealData, model=fit.model,
                silent=TRUE, control=list(maxit=1000))
    out=data.frame(H=names(Z.models)[i], Q=Q.model, U=U.model, B=B.model,
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   m=length(unique(Z.models[[i]])),
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
    if(i==5) next #one m for panmictic so only run 1 Q
  }
}


###################################################
### code chunk number 14: Cs06_sort.results
###################################################
min.AICc=order(out.tab$AICc)
out.tab.1=out.tab[min.AICc,]


###################################################
### code chunk number 15: Cs07_add.delta.aicc
###################################################
out.tab.1=cbind(out.tab.1,
                delta.AICc=out.tab.1$AICc-out.tab.1$AICc[1])


###################################################
### code chunk number 16: Cs08_add.delta.aicc
###################################################
out.tab.1=cbind(out.tab.1, 
                rel.like=exp(-1*out.tab.1$delta.AICc/2))


###################################################
### code chunk number 17: Cs09_aic.weight
###################################################
out.tab.1=cbind(out.tab.1,
                AIC.weight = out.tab.1$rel.like/sum(out.tab.1$rel.like))


###################################################
### code chunk number 18: Cs10_print.table
###################################################
out.tab.1$delta.AICc = round(out.tab.1$delta.AICc, digits=2)
out.tab.1$AIC.weight = round(out.tab.1$AIC.weight, digits=3)
print(out.tab.1[,c("H","Q","B", "delta.AICc","AIC.weight")], row.names=FALSE)


###################################################
### code chunk number 19: Cs11_fignorthsouth
###################################################
best.fit=fits[min.AICc][[1]]
par(mfrow = c(1,1))
matplot(years, t(best.fit$states-best.fit$states[,1]), 
        xlab="abundance index", ylab="",
        type="l",lwd=2,col="black")
legend("topleft",c("Bolinas Lagoon","South (DE/DP)","North (TB/TP)"),lwd=2,lty=c(1:3),bty="n")

## why are these linear!!!  same problem with SealTrend




###################################################
### code chunk number 21: Cs12_new.Q.model  ## takes about 2-3 min to run with example data
###################################################
for(i in 1:length(Z.models)){
  if(i==5) next #don't rerun panmictic
  for(Q.model in c("equalvarcov","unconstrained")){
    fit.model = c(list(Z=Z.models[[i]], Q=Q.model), model.constant)
    fit = MARSS(sealData, model=fit.model,
                silent=TRUE, control=list(maxit=1000))
    out=data.frame(H=names(Z.models)[i], Q=Q.model, U=U.model,
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   m=length(unique(Z.models[[i]])),
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
  }
}


###################################################
### code chunk number 22: Cs13_out.tab.2
###################################################
min.AICc=order(out.tab$AICc)
out.tab.2=out.tab[min.AICc,]
fits=fits[min.AICc]
out.tab.2=cbind(out.tab.2,delta.AICc=out.tab.2$AICc-out.tab.2$AICc[1])
out.tab.2=cbind(out.tab.2,rel.like=exp(-1*out.tab.2$delta.AICc/2))
out.tab.2=cbind(out.tab.2,AIC.weight=out.tab.2$rel.like/sum(out.tab.2$rel.like))


###################################################
### code chunk number 23: Cs14_out.tab.2
###################################################
out.tab.2$AIC.weight = round(out.tab.2$AIC.weight, digits=3)
out.tab.2$delta.AICc = round(out.tab.2$delta.AICc, digits=2)
print(out.tab.2[1:10,c("H","Q","delta.AICc","AIC.weight")], row.names=FALSE)


###################################################
### code chunk number 24: Cs15_equalvarcov.weight
###################################################
c(
  sum(out.tab.2$AIC.weight[out.tab.2$Q=="equalvarcov"]),
  sum(out.tab.2$AIC.weight[out.tab.2$Q=="unconstrained"]),
  sum(out.tab.2$AIC.weight[out.tab.2$Q=="diagonal and equal"])
)


###################################################
### code chunk number 25: Cs16_Q.mat
###################################################
Q.unc=coef(fits[[3]],type="matrix")$Q


###################################################
### code chunk number 26: Cs17_Q.diag
###################################################
diag(Q.unc)


###################################################
### code chunk number 27: Cs18_Q.corr  
###################################################
h=diag(1/sqrt(diag(Q.unc)))
Q.corr=h%*%Q.unc%*%h
rownames(Q.corr)=unique(Z4)
colnames(Q.corr)=unique(Z4)

Q.corr

### let's ge the U's
###################################################
### code chunk number 25: Cs16_Q.mat
###################################################


###------------------
## setting up 3rd loop to vary the 
###################################################
### code chunk number 21: Cs12_new.Q.model  ## takes about 2-3 min to run with example data
###################################################
for(i in 1:length(Z.models)){
  if(i==5) next #don't rerun panmictic
  for(Q.model in c("equalvarcov","unconstrained")){
    fit.model = c(list(Z=Z.models[[i]], Q=Q.model), model.constant)
    fit = MARSS(sealData, model=fit.model,
                silent=TRUE, control=list(maxit=1000))
    out=data.frame(H=names(Z.models)[i], Q=Q.model, U=U.model,
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   m=length(unique(Z.models[[i]])),
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab=rbind(out.tab,out)
    fits=c(fits,list(fit))
  }
}


###################################################
### code chunk number 22: Cs13_out.tab.2
###################################################
min.AICc=order(out.tab$AICc)
out.tab.2=out.tab[min.AICc,]
fits=fits[min.AICc]
out.tab.2=cbind(out.tab.2,delta.AICc=out.tab.2$AICc-out.tab.2$AICc[1])
out.tab.2=cbind(out.tab.2,rel.like=exp(-1*out.tab.2$delta.AICc/2))
out.tab.2=cbind(out.tab.2,AIC.weight=out.tab.2$rel.like/sum(out.tab.2$rel.like))


###################################################
### code chunk number 23: Cs14_out.tab.2
###################################################
out.tab.2$AIC.weight = round(out.tab.2$AIC.weight, digits=3)
out.tab.2$delta.AICc = round(out.tab.2$delta.AICc, digits=2)
print(out.tab.2[1:10,c("H","Q","delta.AICc","AIC.weight")], row.names=FALSE)


###################################################
### code chunk number 24: Cs15_equalvarcov.weight
###################################################
c(
  sum(out.tab.2$AIC.weight[out.tab.2$Q=="equalvarcov"]),
  sum(out.tab.2$AIC.weight[out.tab.2$Q=="unconstrained"]),
  sum(out.tab.2$AIC.weight[out.tab.2$Q=="diagonal and equal"])
)

####---------------- end 3rd model loop 




## going to run best model Z4 standalone to easily get parameters for plotting, etc.
Z.model=Z4
Q.model="diagonal and equal"
U.model="unequal"
R.model="diagonal and equal"
A.model="scaling"
B.model="identity"
x0.model="unequal"
V0.model="zero"
model.constant=list(
  U=U.model, R=R.model, A=A.model, 
  x0=x0.model, V0=V0.model, tinitx=0)
ThreePopFinal = MARSS(sealData, model=model.constant, control=list(maxit=1000))
## get CIs
CIs <- MARSSparamCIs(ThreePopFinal)

## plot resids
plotdat = t(sealData)
matrix.of.biases = matrix(coef(ThreePopFinal,type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
par(mfrow=c(2,3))
for(i in 1:n){
  j=c(1,2,3,4,5)
  xs = ThreePopFinal$states[j[i],]
  resids = plotdat[,i]-matrix.of.biases[,i]-xs
  plot(resids[!is.na(resids)],ylab="residuals")
  title(paste("3 pops", legendnames[i]))
}

##Try plot again
#make figure
par(mfrow=c(1,1))
matplot(years, t(sealData),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,as.data.frame(t(ThreePopFinal$states-1.96*ThreePopFinal$states.se)),type="l",
      lwd=1,lty=2,col="red")
lines(years,ThreePopFinal$states+1.96*ThreePopFinal$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,ThreePopFinal$states,type="l",lwd=2)
title("Observations and total population estimate",cex.main=.9)


## get parameters
options(scipen = 999) ## turn off sci notation (use 1 to turn back on)
coef(ThreePopFinal, type="matrix")$R
coef(ThreePopFinal, type="matrix")$Q
coef(ThreePopFinal, type="matrix")$U
coef(ThreePopFinal, type="matrix")$B

## recover q-matrix
Q.unc=coef(ThreePopFinal,type="matrix")$Q
Q.unc
h=diag(1/sqrt(diag(Q.unc)))
Q.corr=h%*%Q.unc%*%h
print(Q.corr, figs = 3)

## recover R-matrix
R.unc=coef(ThreePopFinal,type="matrix")$R
R.unc
hR=diag(1/sqrt(diag(R.unc)))
R.corr=hR%*%R.unc%*%hR
print(R.corr, figs = 3)


########################################################################

## Parameter key ########
## R = observation errors          equal  (= observation variances)
## U = growth parameter            unequal  
## Z = design matrix 
## Q = hidden state process        diagonal and unequal (= independent year-to-year environmental variability)
## B = effect of column on row     unequal (these are the interactions)

#####################################################################


########----------
## let's get some growth rate summaries from best model


##---------- lets do a covariate model with MEI from Dec-Jan
MEI <- read.csv("C:/bbecker/Projects/phoca_trends/2011-/Phoca Long Term Analysis/Data 1997-2014 PROOFED/MEI.csv")
show(MEI)
##cut to seaData time series
MEI <- dplyr::filter(MEI, Year > 1995)
## make MEI a vector
MEI <- as.vector(MEI[,2])

Z.model=Z4
Q.model="diagonal and equal"
U.model="unequal"
R.model="diagonal and equal"
A.model="scaling"
B.model="identity"
x0.model="unequal"
V0.model="zero"
d.model="MEI"
model.constant=list(
  U=U.model, R=R.model, A=A.model, 
  x0=x0.model, V0=V0.model, d = d.model, tinitx=0)
ThreePopFinal = MARSS(sealData, model=model.constant, control=list(maxit=1000))

##----------------------------------------
## can't get covariate to work....go back to old school glmm/gamm
## set up dataframe
## get MEI
MEI <- read.csv("C:/bbecker/Projects/phoca_trends/2011-/Phoca Long Term Analysis/Data 1997-2014 PROOFED/MEI.csv")
show(MEI)
##cut to seaData time series
MEI <- dplyr::filter(MEI, Year > 1995)
## from SealTrend Script
pup.breed <- top1.pup.breed 
## add MEI
pup.breed <- dplyr::left_join(pup.breed, MEI, by = "Year")
all.breed <- dplyr::left_join(pup.breed, top1.adult.breed, by = c("Year", "Subsite"))
##rename columns
all.breed <- dplyr::rename(all.breed, Pup.Count = Count.x, Adult.Count = Count.y)

## from MARSS model, three pops, so create a nested field for mixed model

pup.breed$region <- ifelse(pup.breed$Subsite == "BL", "BL", 
                           ifelse(pup.breed$Subsite == "DE", "South", 
                                  ifelse(pup.breed$Subsite == "DP", "South", "North")))

## and the glmm
library(lme4)
m1.pup.p <- glmer(Count ~ I(Year-1995) + MEI_DEC_JAN + (1|region/Subsite), data = pup.breed, family = poisson)
summary(m1.pup.p)
par(mfrow = c(1,3))
plot(fitted(m1.pup.p), resid(m1.pup.p))

m1.pup.nb <- glmer(Count ~ I(Year-1995) + MEI_DEC_JAN + (1|region/Subsite), 
                   data = pup.breed, 
                   family = negative.binomial(1))
summary(m1.pup.nb)
plot(fitted(m1.pup.nb), resid(m1.pup.nb))
## nb binomial model best, with region/Subsite nesting better than Subsite
## so maybe better to use binomial # pups over adults?

m1.pup.bin <- glmer(cbind(Pup.Count, Adult.Count) ~ I(Year-1995) + MEI_DEC_JAN * (1|region/Subsite), 
                    data = all.breed, family = binomial)
summary(m1.pup.bin)
plot(fitted(m1.pup.bin), resid(m1.pup.bin))
AIC(m1.pup.p,m1.pup.nb, m1.pup.bin) 
library(effects)
plot(allEffects(m1.pup.bin))
plot(allEffects(m1.pup.nb))

##---------------------------------
## ok, rstanarm models
library(rstanarm)
library(parallel)
library(shinystan)
CORES <- 8
THIN <- 5
ITER = 1000
CHAINS = 3
m1.pup.stan.nb <- stan_glmer.nb(
  Count ~ I(Year - 1995) + MEI_DEC_JAN * Subsite + (1 | region / Subsite),
  data = pup.breed,
  cores = CORES, thin = THIN, iter = ITER, chains = CHAINS
)
launch_shinystan(m1.pup.stan.nb)
## see if correct random effects by removing region (runs way faster!)
m1.pup.stan.nb.2 <- stan_glmer.nb(
  Count ~ I(Year - 1995) + MEI_DEC_JAN + (1 | Subsite),
  data = pup.breed,
  cores = CORES, thin = THIN
)
launch_shinystan(m1.pup.stan.nb.2)

# whoops, binomial had better resids and AIC...
m1.pup.stan.bin <- stan_glmer(
  cbind(Pup.Count, Adult.Count) ~ I(Year - 1995) + MEI_DEC_JAN + (1 | region / Subsite),
  data = all.breed, family = binomial,
  cores = CORES, thin = THIN
)
launch_shinystan(m1.pup.stan.bin)

#####---------------------




#####--------------
## NOT USED YET
###################################################
### code chunk number 28: Cs19_add.hood.canal
###################################################
sealData.hc = rbind(sealData,harborSeal[,8])
rownames(sealData.hc)[12]="Hood.Canal"


###################################################
### code chunk number 29: Cs20_hood.z.models
###################################################
ZH1=factor(c("nc","nc","is","is","ps",
             "ps","sc","sc","nc","sc","is","ps")) 
ZH2=factor(c("nc","nc","is","is","ps",
             "ps","sc","sc","nc","sc","is","hc")) 
Z.models.hc=list(ZH1, ZH2)
names(Z.models.hc)=c("hood.in.ps","hood.separate")


###################################################
### code chunk number 30: Cs21_hood.uqr.models
###################################################
Q3=matrix(list("offdiag"),5,5)
diag(Q3)="q"
Q3[,5]=0; Q3[5,]=0; Q3[5,5]="q.hc"
Q.models=list("equalvarcov","unconstrained",Q3)
names(Q.models)=c("equalvarcov","unconstrained","hood.independent")


###################################################
### code chunk number 31: Cs22_hood-q3
###################################################
Q.models$hood.independent


###################################################
### code chunk number 32: Cs23_out.tab.hc
## runs in about 30 sec
###################################################
out.tab.hc=NULL
fits.hc=list()
for(i in 1:length(Z.models.hc)){
  for(j in 1:length(Q.models)){
    if(i==1 & j==3) next #Q3 is only for Hood Separate model
    Q.model=Q.models[[j]]
    fit.model = c(list(Z=Z.models.hc[[i]], Q=Q.model), model.constant)
    fit = MARSS(sealData.hc, model=fit.model,
                silent=TRUE, control=list(maxit=1000))
    out=data.frame(H=names(Z.models.hc)[i], Q=names(Q.models)[j], U=U.model,
                   logLik=fit$logLik, AICc=fit$AICc, num.param=fit$num.params,
                   m=length(unique(Z.models.hc[[i]])),
                   num.iter=fit$numIter, converged=!fit$convergence)
    out.tab.hc=rbind(out.tab.hc, out)
    fits.hc=c(fits.hc,list(fit))
  }
}


###################################################
### code chunk number 33: Cs24_sort.aicc.hc
###################################################
min.AICc=order(out.tab.hc$AICc)
out.tab.hc=out.tab.hc[min.AICc,]
out.tab.hc=cbind(out.tab.hc, delta.AICc=out.tab.hc$AICc-out.tab.hc$AICc[1])
out.tab.hc=cbind(out.tab.hc,rel.like=exp(-1*out.tab.hc$delta.AICc/2))
out.tab.hc=cbind(out.tab.hc,AIC.weight=out.tab.hc$rel.like/sum(out.tab.hc$rel.like))


###################################################
### code chunk number 34: Cs25_out.tab.2
###################################################
out.tab.hc$AIC.weight = round(out.tab.hc$AIC.weight, digits=3)
out.tab.hc$delta.AICc = round(out.tab.hc$delta.AICc, digits=2)
print(out.tab.hc[,c("H","Q","delta.AICc","AIC.weight")], row.names=FALSE)


