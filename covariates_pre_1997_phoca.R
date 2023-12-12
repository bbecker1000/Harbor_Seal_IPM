#covariates and older seal data


#test comment
library("readxl")
library(dplyr)
library(tidyverse)
library(lubridate)


#Sarah Allen's 1977-1997 phoca counts for pupping and molt
hseals_1977_1997 <- read_excel("Data/Phoca_counts_1977-1997.xlsx")

#Add a month value

hseals_1977_1997$month <- month(hseals_1977_1997$DATE)

hseals_1977_1997_max_pupping <- hseals_1977_1997 %>%
  filter(month <=5) %>%
  group_by(YEAR, AGE) %>%
  slice(which.max(COUNT)) 

hseals_1977_1997_max_pupping <- hseals_1977_1997_max_pupping[,c(1, 3, 7, 8)]
hseals_1977_1997_max_pupping

unique(hseals_1977_1997_max_pupping$SUBSITE)



#merge with modern counts
#logged adult counts
pre_1997_adults_breeding <- hseals_1977_1997_max_pupping %>%
  filter(AGE != "PUP") 

pre_1997_adults_breeding <- pre_1997_adults_breeding[,-3] #remove age

pre_1997_adults_breeding %>%
  pivot_wider(names_from = SUBSITE, values_from = COUNT)


dat2






# MEI
MEI <- read_excel("Data/MEI.xlsx")
show(MEI)
#View(MEI)
##cut to sealData time series

MEI <- MEI %>% filter(Year < 2023)

#if NPGO
#NPGO <- MEI[,c(1,3)]
## make NPGO a vector
#NPGO <- NPGO[,2]
#if Coyote


# no coyote sightings at TB or TP so remove
Coyote <- MEI[,c(4:8)]
small_c <- as.matrix(t(Coyote))


Z.model=Z4 #3 populations
Q.model="diagonal and equal" #default is diagonal and unequal
U.model="unequal" #default
R.model="diagonal and equal" #default
#A.model="scaling" #default
B.model="identity" #default
x0.model="unequal" #???
V0.model="zero" #default
#d.model = "MEI"
#C.model = "unconstrained"
#c.model="MEI"

small_c <- as.matrix(t(NPGO))

model.constant = list(Z=Z.model, Q=Q.model, C="unconstrained", c=small_c)

model.constant = list(Z=Z.model, Q=Q.model, C="unconstrained", safe = TRUE)

ThreePopFinal = MARSS(sealData, model=model.constant, control=list(maxit=1000))

MARSSparamCIs(ThreePopFinal, method = "hessian", alpha = 0.05, nboot = 100, silent = FALSE, hessian.fun = "Harvey1989")




Z.model=factor(c(1,2,2,3,3)) 
U.model="unequal"
Q.model="diagonal and unequal"
R.model="diagonal and unequal" 
B.model= "identity"  
kem3pop2 = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, R=R.model, B=B.model), control=list(maxit=2000, safe=TRUE))
kem3pop2$AIC 

beepr::beep(0)

#add MEI covariate:
#small c in script chap 8...nbeed to fix up.  still has dummy values for 2021 and 2022
kem3pop2.NPGO = MARSS(dat, model=list(Z=Z.model, U=U.model, Q=Q.model, 
                                     R=R.model, B=B.model, c = small_c), 
                          control=list(maxit=5000, safe=TRUE))
kem3pop2.MEI$AIC  
beepr::beep(0)  

MARSSparamCIs(kem3pop2.NPGO, method = "hessian", alpha = 0.05, nboot = 100, silent = FALSE, hessian.fun = "Harvey1989")

  