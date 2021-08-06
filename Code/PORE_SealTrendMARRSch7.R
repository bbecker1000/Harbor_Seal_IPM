##edited code from MARSS UserGuide chapter 7



library(MARSS)
library(ggplot2)


## B in
#  terms of the interaction strengths between species; bij
#  equals dfi /dXj, the change in the log population growth
#  rate of species i with respect to changes in the log
#  population abundance of species j. T

##%######################################################%##
#                                                          #
####                   data wrangling                   ####
#                                                          #
##%######################################################%##

## get PRNS data
Phoca <-
  read.csv("1997_2017_Phocadata.csv")

attach(Phoca) ## needed?
head(Phoca)

## need some year and julian date fields to plot by year and day of year
library(lubridate)
## tell lubridate the format
Phoca$Date2 <- mdy(Phoca$Date)
## add year column for faceting plots
Phoca$Year <- year(Phoca$Date2)
## get Julian Date (yday) for within year dates
Phoca$Julian <- yday(Phoca$Date2)

## Looks like need to remove Dead adult and deadpup
library(plyr)
library(dplyr)
Phoca <- dplyr::filter(Phoca, Age != "DEADADULT")
Phoca <- dplyr::filter(Phoca, Age != "DEADPUP")

## and need to convert HPUP (typo in database!) to PUP

Phoca$Age <- plyr::revalue(Phoca$Age, c("HPUP" = "PUP"))
Phoca$Yearf <- as.factor(Phoca$Year)

## now within year plots to look at breeding and molting season peaks
Phoca.Adult <- dplyr::filter(Phoca, Age == "ADULT")

seasonal.plot.adult <- plot.breeding <- ggplot(Phoca.Adult, aes(Julian, Count, colour = Year)) + geom_point() + geom_smooth(aes(colour = Year))
seasonal.plot.adult + facet_grid(. ~ Subsite)

## previous hard to see patterns due to molt and pups becoming big(= adults), so just look at pups?
Phoca.PUP <- dplyr::filter(Phoca, Age == "PUP")
## subset to the peak seal breeding season about March 20 - April 10
Phoca.breeding <- dplyr::filter(Phoca, Julian > 110 & Julian < 130)

d1 <- Phoca.PUP %>% 
  group_by(Subsite, Year) %>% 
  summarise_each(funs(MaxCount = max), Count)
d1

## julian is meaningless here, just for practice  
d2 <- Phoca.PUP %>% 
  group_by(Subsite, Year) %>% 
  summarise_each(funs(MaxJulian = max), Julian)
d2
## combine d1 and d2
d3 <- dplyr::bind_cols(d2, d1)
d3
## rename duplicate rows
## need to fix this... must have been a dplyr update
d3 <- d3[ -c(4:5)]
d3
## plot just the max
p1<-ggplot(d3, aes(Year, MaxCount)) +
  geom_point() +
  geom_smooth(se = TRUE)
p1 + facet_grid(. ~ Subsite)


Phoca.Adult.Breed <- dplyr::filter(Phoca.Adult, Julian > 110 & Julian < 130)

Phoca.breed.seas <- dplyr::filter(Phoca, Julian > 110 & Julian < 130)
top1 <- tbl_df(Phoca.breed.seas) %>% 
  group_by(Subsite, Yearf, Age) %>%
  top_n(n = 1, wt = Count)

plot.top1 <- ggplot(top1, aes(Year, Count, shape = Age, colour = Age)) + 
  geom_point() + 
  geom_smooth()
plot.top1 + facet_wrap(~ Subsite) + labs(title = "Top 1 Data Point")
ggsave("plot.top1.jpg", width = 8, height = 6, units = "in")

##remove any duplicate rows (some problem in the above queries!)
top1 <- dplyr::distinct(top1)

## now make 3 files of "top 1 for breeding season for pups, adults, and pups + adults (later)
top1.adult.breed <- dplyr::filter(top1, Age == "ADULT")
top1.pup.breed <- dplyr::filter(top1, Age == "PUP")
## MARSS only needs year, site and count
## dplyr:filter won't let me filter out fYear, so use base R
myvars <- c("Year", "Subsite", "Count")
top1.adult.breed <- top1.adult.breed[myvars]
top1.pup.breed <- top1.pup.breed[myvars]

## remove Point Bonita and Duxbury since no pups and some duplicates in the dataset!
top1.adult.breed <- dplyr::filter(top1.adult.breed, Subsite != "DR") ## what is "or" to do in 1-step?
top1.adult.breed <- dplyr::filter(top1.adult.breed, Subsite != "PB")
top1.pup.breed <- dplyr::filter(top1.pup.breed, Subsite != "DR")
top1.pup.breed <- dplyr::filter(top1.pup.breed, Subsite != "PB")
## and for the trial runs, lets get rid of PRH
top1.adult.breed <- dplyr::filter(top1.adult.breed, Subsite != "PR")
top1.pup.breed <- dplyr::filter(top1.pup.breed, Subsite != "PR")

# make sure no duplicate rows
top1.adult.breed <- dplyr::distinct(top1.adult.breed)
top1.pup.breed <- dplyr::distinct(top1.pup.breed)

## log all the counts, # log = ln in R
top1.adult.breed$Count <- log(top1.adult.breed$Count) 
top1.pup.breed$Count <- log(top1.pup.breed$Count)

## now need to make the data wide for MARSS
top1.adult.breed.spread <- tidyr::spread(top1.adult.breed, Subsite, Count)

## still a duplicate rows in the pup one???
top1.pup.breed <- dplyr::distinct(top1.pup.breed)
top1.pup.breed.spread <- tidyr::spread(top1.pup.breed, Subsite, Count)

dat2 <- top1.adult.breed.spread ## for use in SealPopStructure Script

## now use names in the basic MARSS code 
dat <- t(top1.adult.breed.spread)
## OR FOR PUPS
## dat <- t(top1.pup.breed.spread)



###################################################
### code chunk number 11: Cs2_Code1
###################################################
#Code to fit the single population model with i.i.d. errors
#Read in data
## dat=t(harborSealWA) #Transpose since MARSS needs time ACROSS columns
years = dat[1,]     # remove years
n = nrow(dat)-1
dat = dat[2:nrow(dat),]
legendnames = (unlist(dimnames(dat)[1]))


########################################################################

## Parameter key ########
## R = observation errors          equal
## U = growth parameter            unequal  
## Z = design matrix 
## Q = hidden state process        diagonal and unequal
## B = effect of column on row     unequal (these are the interactions)

#####################################################################


###NOW TO MODELS ##########
#estimate parameters
Z.model = factor(c(1,1,1,1,1))
R.model = "diagonal and unequal" 
kem_onepop = MARSS(dat, model=list(Z=Z.model, R=R.model))

#make figure
par(mfrow=c(1,1))
matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem_onepop$states-1.96*kem_onepop$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem_onepop$states+1.96*kem_onepop$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem_onepop$states,type="l",lwd=2)
title("Observations and total population estimate",cex.main=.9)

coef(kem_onepop, type="vector")  #show the estimated parameter elements as a vector
coef(kem_onepop, type="matrix")$R
#show estimated elements for each parameter matrix as a list
coef(kem_onepop) 

kem_onepop$logLik   #show the log-likelihood
kem_onepop$AIC  #show the AIC

#plot residuals
plotdat = t(dat)
matrix.of.biases = matrix(coef(kem_onepop, type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
xs = matrix(kem_onepop$states,
            nrow=dim(plotdat)[1],ncol=dim(plotdat)[2],byrow=F)
resids = plotdat-matrix.of.biases-xs
par(mfrow=c(2,3))
for(i in 1:n){
  plot(resids[!is.na(resids[,i]),i],ylab="residuals")
  title(paste("One Population", legendnames[i]))
}



###################################################
### code chunk number 18: Cs2_Code2
###################################################
#Code to fit the single population model with independent and unequal errors 

Z.model = factor(c(1,1,1,1,1))
R.model = "diagonal and equal" 
kem2 = MARSS(dat, model=list(Z=Z.model, R=R.model))

coef(kem2) #the estimated parameter elements
kem2$logLik #log likelihood
kem2$AIC  #AICs

#plot residuals
plotdat = t(dat)
matrix.of.biases = matrix(coef(kem2, type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
xs = matrix(kem2$states,
            nrow=dim(plotdat)[1],ncol=dim(plotdat)[2],byrow=F)
resids = plotdat-matrix.of.biases-xs
par(mfrow=c(2,3))
for(i in 1:n){
  plot(resids[!is.na(resids[,i]),i],ylab="residuals")
  title(paste("One Population", legendnames[i]))
}
 ## BL and DP residuals look a bit problematic

par(mfrow=c(1,1))

###################################################
### code chunk number 23: Cs2_Code4
###################################################
Z.model=factor(c(1,2,3,4,5))
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and unequal"
B.Model="identity"
kem4_ind=MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.Model),
                    control=list(maxit=1000, safe=TRUE)) ## still not getting convergence! (45 sec to run)

#plot residuals
plotdat = t(dat)
matrix.of.biases = matrix(coef(kem4_ind,type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
par(mfrow=c(2,3))
for(i in 1:n){
  j=c(1,2,3,4,5)
  xs = kem4_ind$states[j[i],]
  resids = plotdat[,i]-matrix.of.biases[,i]-xs
  plot(resids[!is.na(resids)],ylab="residuals")
  title(paste("Independent Populations", legendnames[i]))
}

####### THIS INDEPENDENT HAS BEST RESIDUALS SO FAR!! #######

coef(kem4_ind, type="vector")  #show the estimated parameter elements as a vector
coef(kem4_ind, type="matrix")$R
#show estimated elements for each parameter matrix as a list
coef(kem4_ind) 

## try opepop plot
par(mfrow=c(1,1))
#make figure
matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem4_ind$states-1.96*kem4_ind$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem4_ind$states+1.96*kem4_ind$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem4_ind$states,type="l",lwd=2)
title("Observations and total population estimate",cex.main=.9)

####################################
## let's try ocean vs bay model
###################################
Z.model=factor(c(1,1,2,1,2))
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and unequal"
B.model="unconstrained"
kemOB = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
              control=list(maxit=500, safe=TRUE))

#plot residuals
plotdat = t(dat)
matrix.of.biases = matrix(coef(kemOB,type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
par(mfrow=c(2,3))
for(i in 1:n){
  j=c(1,1,2,1,2)
  xs = kemOB$states[j[i],]
  resids = plotdat[,i]-matrix.of.biases[,i]-xs
  plot(resids[!is.na(resids)],ylab="residuals")
  title(paste("Bay vs Ocean", legendnames[i]))
}

## kem 4 still looks best residual and AIC wise.
coef(kemOB, type="vector")  #show the estimated parameter elements as a vector
coef(kemOB, type="matrix")$R
#show estimated elements for each parameter matrix as a list
coef(kemOB) 

##########################################################
#Two subpopulations with different population parameters  
## (North and South)
##########################################################
Z.model=factor(c(1,1,1,2,2))
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and equal"
B.model="diagonal and unequal"
kemNS = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model) )

#plot residuals
plotdat = t(dat)
matrix.of.biases = matrix(coef(kemNS,type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
par(mfrow=c(2,3))
for(i in 1:n){
  j=c(1,1,1,2,2)
  xs = kemNS$states[j[i],]
  resids = plotdat[,i]-matrix.of.biases[,i]-xs
  plot(resids[!is.na(resids)],ylab="residuals")
  title(paste("North vs South", legendnames[i]))
}



## N and S model not so good
coef(kemNS, type="vector")  #show the estimated parameter elements as a vector
coef(kemNS, type="matrix")$R
#show estimated elements for each parameter matrix as a list
coef(kemNS) 



##########################################################
#Three subpopulations with different population parameters  
## (BL + DE/DP + TP/DE)
##########################################################
Z.model=factor(c(1,2,2,3,3))
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and equal"
B.model= "identity"   # "diagonal and unequal"
kem3pop = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model) )

matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3", "4", "5"), ylim=c(3,9), bty="L")
lines(years,kem3pop$states-1.96*kem3pop$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states+1.96*kem3pop$states.se,type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states,type="l",lwd=2)
title("Observations and total population estimate",cex.main=.9)

plot(years,kem3pop$states[1,])


plotdat = t(dat)
matrix.of.biases = matrix(coef(kem3pop,type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
par(mfrow=c(2,3))
for(i in 1:n){
  j=c(1,2,2,3,3)
  xs = kem3pop$states[j[i],]
  resids = plotdat[,i]-matrix.of.biases[,i]-xs
  plot(resids[!is.na(resids)],ylab="residuals")
  title(paste("3 pops", legendnames[i]))
}

## get parameters
coef(kem3pop, type="matrix")$R
coef(kem3pop, type="matrix")$Q
coef(kem3pop, type="matrix")$U
coef(kem3pop, type="matrix")$B

## AICs
c(kem_onepop$AIC,kem2$AIC,kem4_ind$AIC,kemNS$AIC, kemOB$AIC, kem3pop$AIC)



############################
# leftover from vignette ## NOT USED
#Hood Canal covaries with the other regions
Z.model=factor(c(1,1,1,1,2))
U.model="unequal"
Q.model="equalvarcov"
R.model="diagonal and unequal"
kem = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model) )


##%######################################################%##
#                                                          #
####                 OK, now got models                 ####
####         working. On to hypothesis testing          ####
#                                                          #
##%######################################################%##





