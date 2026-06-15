#adults and molts all years models



## Parameter key ########
## Z = design matrix = Spatial population structure
## R = observation errors          equal
## U = growth parameter            (un)equal, or matrix [Z x 1]
## Q = hidden state process        diagonal and (un)equal
#                                  equalvarcov or unconstrained
## B = effect of column on row     unequal (these are the interactions)

#####################################################################

## B in -------------------------------------------
#  terms of the interaction strengths between species; bij
#  equals dfi /dXj, the change in the log population growth
#  rate of species i with respect to changes in the log
#  population abundance of species j. T

#Data structure----------------------------------

# Z models/Hypotheses ---------------------------
Z.model = factor(rep(c(1:3), times = 8))
#R.model = "diagonal and equal"
R.model=diag(0.025, 24)
Q.model=diag(0.025, 3)

AllYearAllClass.Z.class=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model))
  #x0 = x0.model, 
  #tinitx=1, 


Z.model = factor(rep(c(1:3), each = 8))
R.model=diag(0.025, 24)
Q.model=diag(0.025, 3)
AllYearAllClassZ.Site=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model))
#x0 = x0.model, 
#tinitx=1, 

Z.model = factor(1:24)
R.model=diag(0.025, 24)
Q.model=diag(0.025, 24)
AllYearAllClass=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model))

#timevarying u by site and class
#timevarying u by site and class
U1 <- matrix(c("t1_BL_A","t1_BL_M","t1_BL_P",
               "t1_DE_A","t1_DE_M","t1_DE_P",
               "t1_DP_A","t1_DP_M","t1_DP_P",
               "t1_DR_A","t1_DR_M","t1_DR_P",
               "t1_PB_A","t1_PB_M","t1_PB_P",
               "t1_PRH_A","t1_PRH_M","t1_PRH_P",
               "t1_TB_A","t1_TB_M","t1_TB_P",
               "t1_TP_A","t1_TP_M","t1_TP_P"),24, 1)
U2 <- matrix(c("t2_BL_A","t2_BL_M","t2_BL_P",
               "t2_DE_A","t2_DE_M","t2_DE_P",
               "t2_DP_A","t2_DP_M","t2_DP_P",
               "t2_DR_A","t2_DR_M","t2_DR_P",
               "t2_PB_A","t2_PB_M","t2_PB_P",
               "t2_PRH_A","t2_PRH_M","t2_PRH_P",
               "t2_TB_A","t2_TB_M","t2_TB_P",
               "t2_TP_A","t2_TP_M","t2_TP_P"),24, 1)
Ut.Site.Class <- array(U2, dim = c(dim(U1), dim(dat)[2]))
TT <- dim(dat)[2]
Ut.Site.Class[, , 1:ceiling(TT / 2)] <- U1


Z.model = factor(1:24)
R.model=diag(0.025, 24)
Q.model=diag(0.025, 24)
cntl.list <- list(maxit = 3000, allow.degen = FALSE)
AllYearAllClass.Ut=MARSS(dat, model=list(
  Z=Z.model, 
  R=R.model,
  Q=Q.model,
  U=Ut.Site.Class),
  control = cntl.list)

beepr::beep(0)

#look at AICs
AllYearAllClass$AIC
AllYearAllClass.Z.class$AIC
AllYearAllClassZ.Site$AIC
AllYearAllClass.Ut$AIC

autoplot(AllYearAllClass.Ut)

#use autoplot shortcut to gather data for plotting
p1<-autoplot(AllYearAllClass.Ut, plot.type = "fitted.ytT") +
  ylim(-2,8) +
  #scale_y_continuous(trans = exp(AllYearAllClass$ytT)) +
  ylab("estimate (log count)")# xtT


#now access the p1 object data
rownames <-as.character(p1[["data"]][[".rownames"]])
t <- p1[["data"]][["t"]]
y <- p1[["data"]][["y"]]
fitted <- p1[["data"]][[".fitted"]]
se <- p1[["data"]][[".se"]]
conf.low <- p1[["data"]][[".conf.low"]]
conf.up <- p1[["data"]][[".conf.up"]]

df <- data.frame(rownames, t, y, fitted, se, conf.low, conf.up)
head(df)

df$Site <- sub("_.*", "", df$rownames)
unique(df$Site)
df$Class <- sub(".*_", "", df$rownames)
head(df)
df$Year <- df$t+1974

#upper CIs for DR molt are unreasonable due to no pre 1997 data.
#remove predictions for DR and PB prior to 1997 since no data
df$y <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$y)
df$y <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$y)
df$fitted <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$fitted)
df$fitted <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$fitted)
df$conf.low <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$conf.low)
df$conf.low <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$conf.low)
df$conf.up <- ifelse(df$Site == "DR" & df$Year <= 1997, NA, df$conf.up)
df$conf.up <- ifelse(df$Site == "PB" & df$Year <= 1997, NA, df$conf.up)

#plot MARSS predictions
ggplot(df, aes(Year, exp(y), fill = Class, color = Class)) +
  geom_ribbon(aes(ymin = exp(conf.low), ymax = exp(conf.up), fill = Class), 
              alpha = 0.3, colour = NA) + 
  geom_point(alpha = 0.75) +
  geom_line(aes(x = Year, y = exp(fitted)), size = 1) +
  ylim(0,2000) +
  ylab("Seals") + 
  scale_x_continuous(breaks=seq(1980,2020,20)) +
  theme_gray(base_size = 18) +
  facet_wrap(.~Site, ncol = 4)


#total pop plots.
#exponentiate the estimates
df.total <- df
head(df.total)
df.total$y <- exp(df.total$y)
df.total$fitted <- exp(df.total$fitted)
df.total$conf.low <- exp(df.total$conf.low)
df.total$conf.up <- exp(df.total$conf.up)

#sum by year and class  
df.total.sum <- 
  df.total %>%
    dplyr::group_by(Class, Year) %>%
    dplyr::summarize(y_sum=sum(y, na.rm = T),
              fitted_sum=sum(fitted, na.rm = T),
              conf.low_sum=sum(conf.low, na.rm = T),
              conf.up_sum=sum(conf.up, na.rm = T)) 


#plot - recall that errors for DR and PB < 1997 have been removed
ggplot(df.total.sum, aes(Year, y_sum, color = Class)) +
  geom_ribbon(aes(ymin = conf.low_sum, ymax = conf.up_sum, fill = Class), 
              alpha = 0.3, colour = NA) +
  #geom_point() +
  geom_line(aes(x = Year, y = fitted_sum), size = 1) +
  ylim(0, 6000) +
  ylab("Seals") + 
  scale_x_continuous(breaks=seq(1980,2020,20)) +
  theme_gray(base_size = 18)
  

#Plot the U's from best model
AllYearAllClass.Ut
CIs <- MARSSparamCIs(AllYearAllClass.Ut, alpha = 0.89)
CIs

CIs <- tidy(CIs)
U <- CIs %>% filter(str_detect(term, "^U"))
U$time <- rep(c(1:2), each = 24)
#get the class from term
U$Class <- sub(".*_", "", U$term)
U$Class2 <- ifelse(U$Class == "A", "Adult",
                  ifelse(U$Class == "M", "Molt", "Pup" ))

#get the site from term
U$Site = sapply(strsplit(U$term, "_"), function(x) x[2])

#make a forest plot.
Time_Periods <- c(
  `1` = "1975-2004",
  `2` = "2005-2023")

#parameter plot of U'ss
ggplot(U, aes(Site, estimate, color = Class2)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), 
                  position = position_dodge(width = 0.9), size = 0.75) +
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Annual rate of change") + 
  theme_gray(base_size = 18) +
  ylim(-0.1, 0.3) +
  facet_wrap(.~time, labeller = as_labeller(Time_Periods))
