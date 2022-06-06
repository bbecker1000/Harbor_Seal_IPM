##edited code from MARSS UserGuide chapter 7

library(MARSS)
library(ggplot2)
library("readxl")

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
Phoca <- read_excel("Data/1997_2019_Phocadata.xls")

Phoca <- Phoca[-c(3:4, 8:10)]

## need some year and julian date fields to plot by year and day of year
library(lubridate)
## tell lubridate the format
Phoca$Date2 <- ymd(Phoca$Date) #separates out month,day and year
## add year column for faceting plots
Phoca$Year <- year(Phoca$Date2)
## get Julian Date (yday) for within year dates
Phoca$Julian <- yday(Phoca$Date2)

## Looks like need to remove Dead adult and deadpup
library(plyr)
library(dplyr)

Phoca <- subset(Phoca, Age != "DEADPUP" | Age != "DEADADULT")

## and need to convert HPUP (typo in database!) to PUP
Phoca$Age[Phoca$Age == "HPUP"] <- "PUP"
Phoca$Yearf <- as.factor(Phoca$Year)

## now within year plots to look at breeding and molting season peaks
Phoca.Adult <- dplyr::filter(Phoca, Age == "ADULT")

seasonal.plot.adult <- plot.breeding <- ggplot(Phoca.Adult, aes(Julian, Count, colour = Year)) + geom_point() + geom_smooth(aes(colour = Year))
seasonal.plot.adult + facet_grid(. ~ Subsite)

## previous hard to see patterns due to molt and pups becoming big(= adults), so just look at pups?
Phoca.PUP <- dplyr::filter(Phoca, Age == "PUP")
## subset to the peak seal breeding season about April 20 - May 10
# (Silas) should this be April 15 to May 15 or April 20 to May 10??  right now using wider interval
Phoca.breeding <- dplyr::filter(Phoca, Julian > 105 & Julian < 135)

if (FALSE) { #remove this line and 80 to get this plot again
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
d3 <- rename(d3, Subsite = Subsite...1) 
d3 <- rename(d3, Year = Year...2)
## plot just the max
p1<-ggplot(d3, aes(Year, MaxCount)) +
  geom_point() +
  geom_smooth(se = TRUE)
p1 + facet_grid(. ~ Subsite) 
}

Phoca.Adult.Breed <- dplyr::filter(Phoca.Adult, Julian > 105 & Julian < 135) #April 15 to May 15

Phoca.breed.seas <- dplyr::filter(Phoca, Julian > 105 & Julian < 135)
top1 <- tbl_df(Phoca.breed.seas) %>% 
  group_by(Subsite, Yearf, Age) %>%
  top_n(n = 1, wt = Count)

##remove any duplicate rows (some problem in the above queries!)
top1 <- dplyr::distinct(top1)

#(silas) pr and pb missing lines connecting the dots 
plot.top1 <- ggplot(top1, aes(Year, Count, shape = Age, colour = Age)) + 
  geom_point() + 
  geom_smooth()
plot.top1 + facet_wrap(~ Subsite) + labs(title = "Top 1 Data Point") #PB and PR don't have line; remove deadpup??
#ggsave("plot.top1.jpg", width = 8, height = 6, units = "in")

## now make 3 files of "top 1 for breeding season for pups, adults, and pups + adults (later)
top1.adult.breed <- dplyr::filter(top1, Age == "ADULT")
top1.pup.breed <- dplyr::filter(top1, Age == "PUP")
## MARSS only needs year, site and count
## dplyr:filter won't let me filter out fYear, so use base R
myvars <- c("Year", "Subsite", "Count")
top1.adult.breed <- top1.adult.breed[myvars]
top1.pup.breed <- top1.pup.breed[myvars]

## remove Point Bonita and Duxbury since no pups and some duplicates in the dataset!
top1.adult.breed <- subset(top1.adult.breed, Subsite != "DR" & Subsite != "PB")
top1.pup.breed <- subset(top1.pup.breed, Subsite != "DR" & Subsite != "PB")
## and for the trial runs, lets get rid of PRH
#(Silas) should we add this back in ???
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
top1.pup.breed.spread <- tidyr::spread(top1.pup.breed, Subsite, Count)

dat2 <- top1.adult.breed.spread ## for use in SealPopStructure Script

## now use names in the basic MARSS code 
dat <- t(top1.adult.breed.spread)
## OR FOR PUPS
## dat <- t(top1.pup.breed.spread)


##########################################################################
### 8.3 Single population model with independent and non-identical errors
##########################################################################
#Code to fit the single population model with i.i.d. errors
#Read in data
#Transpose since MARSS needs time ACROSS columns
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
#coef(kem_onepop) 

kem_onepop$logLik   #show the log-likelihood
kem_onepop$AIC  #show the AIC

#plot residuals, from 8.3 (page 99)
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
#bad residuals: BL, DP
#good-ish residuals: DE, TB, TP(?)


##########################################################
### 8.2 A single well-mixed population with i.i.d. errors
##########################################################
#Code to fit the single population model with independent and equal errors
#(silas) model is linear and BAD if R. model = "diagonal and equal" changed to unequal 6/6
#(silas) now model follows similar patterns throughout, works better for some (TB + CI for example) than others (bad for BL + CI)

Z.model = factor(c(1,1,1,1,1))
R.model = "diagonal and unequal" 
kem2 = MARSS(dat, model=list(Z=Z.model, R=R.model))

coef(kem2) #the estimated parameter elements
kem2$logLik #log likelihood
kem2$AIC  #AICs

#plot residuals from page 99
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
 ## (Silas) I think most look somewhat problematic, maybe TB and TP are ok

par(mfrow=c(1,1))

###################################################
### code chunk number 23: Cs2_Code4
### Each population is independent

#(silas) BL and DE clearly nonlinear but DP, TB, and TP aren't
#(silas) is this because of the convergence problems??
#(Silas) tried running with maxit=10000 and it didn't change the plots (except the confidence estimates for TB??)

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
# (Silas) I don't like DP and TB too much but definite improvement

coef(kem4_ind, type="vector")  #show the estimated parameter elements as a vector
coef(kem4_ind, type="matrix")$R
#show estimated elements for each parameter matrix as a list
coef(kem4_ind) 

## try opepop plot
par(mfrow=c(2,3))
#make figure
matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem4_ind$states[1,]-1.96*kem4_ind$states.se[1,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[1,]+1.96*kem4_ind$states.se[1,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[1,],type="l",lwd=2)
title("Observations and total population estimate for site 1",cex.main=.9)

matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem4_ind$states[2,]-1.96*kem4_ind$states.se[2,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[2,]+1.96*kem4_ind$states.se[2,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[2,],type="l",lwd=2, col="red")
title("Observations and total population estimate for site 2",cex.main=.9)

matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem4_ind$states[3,]-1.96*kem4_ind$states.se[3,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[3,]+1.96*kem4_ind$states.se[3,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[3,],type="l",lwd=2, col="green")
title("Observations and total population estimate for site 3",cex.main=.9)

matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem4_ind$states[4,]-1.96*kem4_ind$states.se[4,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[4,]+1.96*kem4_ind$states.se[4,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[4,],type="l",lwd=2, col="blue")
title("Observations and total population estimate for site 4",cex.main=.9)
matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3","4","5"), ylim=c(3,9), bty="L")
lines(years,kem4_ind$states[5,]-1.96*kem4_ind$states.se[5,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[5,]+1.96*kem4_ind$states.se[5,],type="l",
      lwd=1,lty=2,col="purple")
lines(years,kem4_ind$states[5,],type="l",lwd=2, col="cadetblue1")
title("Observations and total population estimate for site 5",cex.main=.9)

#(silas)
#above shows the same thing as calling plot(kem4_ind) and looking just at the first plot
#still see that only BL and DE are nonlinear

####################################
## let's try ocean vs bay model
## Most similar to 8.4 Two Sub-populations
###################################
#(silas) ALL PLOTS ARE NONLINEAR YAY
#only looks good for DE and TB 
#why does the line look like it's following the same pattern for TP and DP but has different intercepts?

Z.model=factor(c(1,1,2,1,2))
U.model="unequal" #default is unequal, not necessary to include
Q.model="diagonal and unequal"
R.model="diagonal and unequal"
B.model="unconstrained" #what's going on with this?? only reference is in chapter 9
kemOB = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model),
              control=list(maxit=500, safe=TRUE)) #getting convergence warning

#plot residuals
plotdat = t(dat)
par(mfrow=c(2,3))
matrix.of.biases = matrix(coef(kemOB,type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)

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
#(silas) no plots are nonlinear

Z.model=factor(c(1,1,1,2,2))
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and equal"
B.model="diagonal and unequal"
kemNS = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model) ) #also convergence warning.. something to do with Q and B??

#plot residuals
par(mfrow=c(2,3))
plotdat = t(dat)
matrix.of.biases = matrix(coef(kemNS,type="matrix")$A,
                          nrow=nrow(plotdat),ncol=ncol(plotdat),byrow=T)
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
## (BL + DE/DP + TP/TB)
##########################################################
#(silas) i get convergence warnings
#only plot looks good for BL. all others are strictly linear with R unequal. all lin if R equal

Z.model=factor(c(1,2,2,3,3)) 
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and equal" #this makes the plots nonlinear, diagonal and equal makes them linear
B.model= "identity"   # "diagonal and unequal"
kem3pop = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model) )

par(mfrow=c(1,3))
matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3", "4", "5"), ylim=c(3,9), bty="L")
lines(years,kem3pop$states[1,]-1.96*kem3pop$states.se[1,],type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states[1,]+1.96*kem3pop$states.se[1,],type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states[1,],type="l",lwd=2)
title("Bolinas Lagoon",cex.main=.9)

matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3", "4", "5"), ylim=c(3,9), bty="L")
lines(years,kem3pop$states[2,]-1.96*kem3pop$states.se[2,],type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states[2,]+1.96*kem3pop$states.se[2,],type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states[2,],type="l",lwd=2)
title("Drakes Estero, Double Point",cex.main=.9)

matplot(years, t(dat),xlab="",ylab="index of log abundance",
        pch=c("1","2","3", "4", "5"), ylim=c(3,9), bty="L")
lines(years,kem3pop$states[3,]-1.96*kem3pop$states.se[3,],type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states[3,]+1.96*kem3pop$states.se[3,],type="l",
      lwd=1,lty=2,col="red")
lines(years,kem3pop$states[3,],type="l",lwd=2)
title("Tomales Bay and Tomales Point",cex.main=.9)

par(mfrow=c(1,1))
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

#Best AIC is 3 groups, then all independent 


##NEW PLOTS
df <- as.data.frame(t(dat))
mod <- as.data.frame(t(kem3pop$states))
names(mod) <- c("BL_m", "DEDP_m", "TBTP_m")
se <- as.data.frame(t(kem3pop$states.se))
df$BL_m <- mod$BL_m
df$BL_lower <- df$BL_m-1.96*se[,1]
df$BL_upper <- df$BL_m+1.96*se[,1]
df$year <- 1996:2019
ggplot(df) + geom_point(aes(years, BL)) + geom_line(aes(years, BL_m)) + geom_line(aes(years, BL_lower), linetype="dashed", color="red") + geom_line(aes(years, BL_upper), linetype="dashed", color="red") + labs(title="Bolinas Lagoon", y="index of log abundance")

df$DEDP_m <- mod$DEDP_m
df$DEDP_lower <- df$DEDP_m-1.96*se[,2]
df$DEDP_upper <- df$DEDP_m+1.96*se[,2]

ggplot(df) + geom_point(aes(years, DE), colour="salmon") + geom_point(aes(years, DP), colour="blue") + geom_line(aes(years, DEDP_m)) + geom_line(aes(years, DEDP_lower), linetype="dashed", color="red") + geom_line(aes(years, DEDP_upper), linetype="dashed", color="red") + labs(title="Drakes Estero and Double Point", y="index of log abundance") + theme(legend.title = element_text(face = "bold"))

df$TBTP_m <- mod$TBTP_m
#SE is all 0 for some reason 
df$TBTP_lower <- df$TBTP_m-1.96*se[,3]
df$TBTP_upper <- df$TBTP_m-1.96*se[,3]
ggplot(df) + geom_point(aes(years, TB), colour="salmon") + geom_point(aes(years, TP), colour="blue") + geom_line(aes(years, TBTP_m)) + labs(title="Tomales Bay and Tomales Point", y="index of log abundance") #+ geom_line(aes(years, TBTP_lower), linetype="dashed", color="red") + geom_line(aes(years, TBTP_upper), linetype="dashed", color="red") 

############################
# leftover from vignette ## NOT USED
#Hood Canal covaries with the other regions
if (FALSE) {
Z.model=factor(c(1,1,1,1,2))
U.model="unequal"
Q.model="equalvarcov"
R.model="diagonal and unequal"
kem = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model) )

}
##%######################################################%##
#                                                          #
####                 OK, now got models                 ####
####         working. On to hypothesis testing          ####
#                                                          #
##%######################################################%##


