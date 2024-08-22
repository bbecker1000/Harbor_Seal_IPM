## all years adult molt PLOTS


#############################
# Model plots -------------------------------------
#############################
#BESTMODEL <- m.1997.2023.06.ut.Class.MOCI.dist
BESTMODEL <-   m.1997.2023.06.ut.Class.MOCI.ES.dist  

#BESTMODEL <-   m.1997.2023.06.ut.Class.All.MOCI.ES.dist


# m.1997.2023.06.ut.Class.MOCI.ES 
# m.1997.2023.06.ut.Site.Class.MOCI.ES 
# m.1997.2023.06.ut.Class.MOCI.ES
#  
#BESTMODEL <- m.1997.2023.06.ut.All.MOCI.ES # restricts trends... unrealistic
#m.1997.2022.06.ut.All.MOCI.ES # m.1997.2022.06.ut.All.MOCI #m.1997.2022.06.ut.All #m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U# m.1pop_Molt_Adult_Coyote_PDO_B_unc_2U #m.1pop_Coyote_PDO_Xo_fixed_B_unc_tinitx_1 #m.1pop_Coyote_PDO_Xo_fixed_B_unc# m.1pop_Coyote_PDO_B_unc #   m.1pop_Coyote_PDO_B_unc   m.5pop_Coyote_PDO_B_unc


autoplot(BESTMODEL)

autoplot(BESTMODEL, plot.type = "fitted.ytT") + # xtT
  ylim(3,8) +
  ylab("estimate (log count)")# xtT

autoplot(BESTMODEL, plot.type = "xtT") + # xtT
  ylim(3.1,7.6)

# plot temporal change from starting values

#BESTMODEL <- m.MOLT_IND_SITE_MOCI.equal_B_diag_uneq

CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11)  #crashes if in results # start 420

#CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11, hessian.fun = "fdHess")  #runs if few NAs, about 60 min

CIs

CIs$states.se  # change this to the CIs rather than SEs...  use 89%CIs for everything?



d <- as_tibble(t(BESTMODEL$states-BESTMODEL$states[,1]))
d.se <- as_tibble(t(BESTMODEL$states.se))
#add header names
names(d)[c(1:18)] <- c("BL_Breed", "BL_Molt", "BL_Pup",
                       "DE_Breed", "DE_Molt", "DE_Pup",
                       "DP_Breed", "DP_Molt", "DP_Pup",
                       "PRH_Breed", "PRH_Molt", "PRH_Pup",
                       "TB_Breed", "TB_Molt", "TB_Pup",
                       "TP_Breed", "TP_Molt", "TP_Pup"
                       # "DR_Breed", "DR_Molt", #"DR_Pup",
                       #  "PB_Breed", "PB_Molt" #"PB_Pup"
)

names(d.se)[c(1:18)] <- c("BL_Breed", "BL_Molt", "BL_Pup",
                          "DE_Breed", "DE_Molt", "DE_Pup",
                          "DP_Breed", "DP_Molt", "DP_Pup",
                          "PRH_Breed", "PRH_Molt", "PRH_Pup",
                          "TB_Breed", "TB_Molt", "TB_Pup",
                          "TP_Breed", "TP_Molt", "TP_Pup"
                          # "DR_Breed", "DR_Molt", #"DR_Pup",
                          #  "PB_Breed", "PB_Molt" #"PB_Pup"
)


d <- as_tibble(cbind(years,d))

d2 <- d %>% pivot_longer(cols = c(2:19), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d.se %>% pivot_longer(cols = c(1:18), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d2.se[,2]
colnames(d2.se)[1] = "log_se"

d4 <- as_tibble(cbind(d2,d2.se))

#add groupings for site and molt
d5 <- d4 %>% separate_wider_delim(Subsite_Season, "_", names = c("Subsite", "Season"))
d5$Subsite_Season <- paste0(d5$Subsite, "_", d5$Season)
d5

ggplot(d5, aes(x = years, y = log_est, group = Subsite_Season, color = Subsite)) +
  #geom_line(aes(linetype = Season), size = 1.25) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = log_est-log_se, ymax = log_est+log_se),  
              alpha = 0.2, colour = NA) +
  geom_hline(yintercept = c(-1,0,1), lty = 2) +
  xlim(1997, 2023) +
  ylim(-1.5, 2) + 
  theme_minimal(base_size = 20) +
  ylab("index of log abundance (N-No)") +
  xlab("Year") +
  facet_wrap(.~Season, ncol = 1)



## plot pup:adult ratios using predicted values

d <- as_tibble(t(BESTMODEL$states)) #-BESTMODEL$states[,1]))
d.se <- as_tibble(t(BESTMODEL$states.se))
#add header names
names(d)[c(1:18)] <- c("BL_Breed", "BL_Molt", "BL_Pup",
                       "DE_Breed", "DE_Molt", "DE_Pup",
                       "DP_Breed", "DP_Molt", "DP_Pup",
                       "PRH_Breed", "PRH_Molt", "PRH_Pup",
                       "TB_Breed", "TB_Molt", "TB_Pup",
                       "TP_Breed", "TP_Molt", "TP_Pup"
                       # "DR_Breed", "DR_Molt", #"DR_Pup",
                       #  "PB_Breed", "PB_Molt" #"PB_Pup"
)

names(d.se)[c(1:18)] <- c("BL_Breed", "BL_Molt", "BL_Pup",
                          "DE_Breed", "DE_Molt", "DE_Pup",
                          "DP_Breed", "DP_Molt", "DP_Pup",
                          "PRH_Breed", "PRH_Molt", "PRH_Pup",
                          "TB_Breed", "TB_Molt", "TB_Pup",
                          "TP_Breed", "TP_Molt", "TP_Pup"
                          # "DR_Breed", "DR_Molt", #"DR_Pup",
                          #  "PB_Breed", "PB_Molt" #"PB_Pup"
)


d <- as_tibble(cbind(years,d))

d2 <- d %>% pivot_longer(cols = c(2:19), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d.se %>% pivot_longer(cols = c(1:18), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d2.se[,2]
colnames(d2.se)[1] = "log_se"

d4 <- as_tibble(cbind(d2,d2.se))

#add groupings for site and molt
d5 <- d4 %>% separate_wider_delim(Subsite_Season, "_", names = c("Subsite", "Season"))
d5$Subsite_Season <- paste0(d5$Subsite, "_", d5$Season)
d5

# remove molt
d6 <- d5 %>% filter(Season != "Molt")
d6 <- d6 %>% select(-Subsite_Season)

# Create hi and lo ratios
d7 <- d6 %>% pivot_wider(names_from = Season, values_from = c(log_est, log_se))
d7$ratio <- d7$log_est_Pup / d7$log_est_Breed
d7$ratio.lo <- (d7$log_est_Pup - (d7$log_est_Pup * d7$log_se_Pup * 1.68)) / 
               (d7$log_est_Breed - (d7$log_est_Breed * d7$log_se_Breed * 1.68) )

d7$ratio.hi <- (d7$log_est_Pup + (d7$log_est_Pup * d7$log_se_Pup * 1.68)) / 
  (d7$log_est_Breed + (d7$log_est_Breed * d7$log_se_Breed * 1.68) )

mean(d7$ratio)


ggplot(d7, aes(x = years, y = ratio)) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = ratio.lo, ymax = ratio.hi), alpha = 0.2, colour = NA) +
  geom_hline(yintercept = 0.845, lty = 2) +
  xlim(1997, 2023) +
  ylim(0.6, 1.1) + 
  theme_grey(base_size = 20) +
  ylab("pups:adults") +
  xlab("Year") +
  facet_wrap(.~Subsite, ncol = 3)




## plot disturbance rates through time
 #from 4b
HumanDisturbance.A <- HumanDisturbance %>% filter(SiteCode != "DR" & SiteCode != "PB")

ggplot(HumanDisturbance.A, aes(x = Year, y = DistRate)) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = mean(HumanDisturbance$DistRate), lty = 2) +
  xlim(1997, 2023) +
  #ylim(0, 1) + 
  theme_grey(base_size = 20) +
  ylab("Anthropogenic distrubance rate") +
  xlab("Year") +
  facet_wrap(.~SiteCode, ncol = 3)



####----plot overall pop size for each age class summing sites

d.tot <- as_tibble(t(exp(BESTMODEL$states)))  #exp to original scale
d.tot.se <- as_tibble(t(BESTMODEL$states.se)) #exp to original scale
#add header names

names(d.tot)[c(1:18)] <- c("BL_Breed", "BL_Molt", "BL_Pup",
                           "DE_Breed", "DE_Molt", "DE_Pup",
                           "DP_Breed", "DP_Molt", "DP_Pup",
                           "PRH_Breed", "PRH_Molt", "PRH_Pup",
                           "TB_Breed", "TB_Molt", "TB_Pup",
                           "TP_Breed", "TP_Molt", "TP_Pup"
                           # "DR_Breed", "DR_Molt", #"DR_Pup",
                           #  "PB_Breed", "PB_Molt" #"PB_Pup"
)
names(d.tot.se)[c(1:18)] <- c("BL_Breed", "BL_Molt", "BL_Pup",
                              "DE_Breed", "DE_Molt", "DE_Pup",
                              "DP_Breed", "DP_Molt", "DP_Pup",
                              "PRH_Breed", "PRH_Molt", "PRH_Pup",
                              "TB_Breed", "TB_Molt", "TB_Pup",
                              "TP_Breed", "TP_Molt", "TP_Pup"
                              # "DR_Breed", "DR_Molt", #"DR_Pup",
                              #  "PB_Breed", "PB_Molt" #"PB_Pup"
)


#sum all molt and breeding sites by year
Breed.tot <- d.tot %>%
  dplyr::reframe(Breed = rowSums(across(c(BL_Breed, DE_Breed, DP_Breed, PRH_Breed, TB_Breed, TP_Breed))))
Molt.tot <- d.tot %>%
  dplyr::reframe(Molt = rowSums(across(c(BL_Molt, DE_Molt, DP_Molt, PRH_Molt, TB_Molt, TP_Molt))))
Pup.tot <- d.tot %>%
  dplyr::reframe(Pup = rowSums(across(c(BL_Pup, DE_Pup, DP_Pup, PRH_Pup, TB_Pup, TP_Pup))))

Breed.tot.se <- d.tot.se %>%
  dplyr::reframe(Breed.se = rowSums(across(c(BL_Breed, DE_Breed, DP_Breed, PRH_Breed, TB_Breed, TP_Breed)))/6) #get mean
Molt.tot.se <- d.tot.se %>%
  dplyr::reframe(Molt.se = rowSums(across(c(BL_Molt, DE_Molt, DP_Molt, PRH_Molt, TB_Molt, TP_Molt)))/6) #get mean
Pup.tot.se <- d.tot.se %>%
  dplyr::reframe(Pup.se = rowSums(across(c(BL_Pup, DE_Pup, DP_Pup, PRH_Pup, TB_Pup, TP_Pup)))/6) 


#put in a single table
d.Breed.Molt.Pup.tot <- as_tibble(cbind(years,Breed.tot, Molt.tot, Pup.tot))


d.Breed.Molt.Pup.tot <- d.Breed.Molt.Pup.tot %>% pivot_longer(cols = c(2:4), names_to = "Season", values_to = "estimate")

d.Breed.Molt.Pup.se <- as_tibble(cbind(years,Breed.tot.se, Molt.tot.se, Pup.tot.se))



d.Breed.Molt.Pup.se <- d.Breed.Molt.Pup.se %>% pivot_longer(cols = c(2:4), names_to = "Season", values_to = "SE")
d.Breed.Molt.Pup.se <- d.Breed.Molt.Pup.se[,3] #remove season and year

d4 <- as_tibble(cbind(d.Breed.Molt.Pup.tot,d.Breed.Molt.Pup.se))


#total pop plot with ses
ggplot(d4, aes(x = years, y = estimate, color = Season)) +
  #geom_line(aes(linetype = Season), size = 1.25) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1.1) +
  # CI at 90% use SE = 1.645, 80% = 1.28
  geom_ribbon(aes(ymin = estimate-1.28*SE*estimate, ymax = estimate+1.28*SE*estimate, fill = Season),  
              alpha = 0.2, color = NA) +
  #geom_hline(yintercept = c(-1,0,1), lty = 2) +
  xlim(1997, 2023) +
  ylim(0, 5500) + 
  theme_classic(base_size = 20) +
  ylab("Estimated abundance") +
  xlab("Year") 



## change over time ---

# get data
## calculate percentage change in total pop from 1982-2003 and 2004-2022

#column names
new_names <- c("STARTPOP", "ENDPOP", "STARTSE", "ENDSE")

Change_1997_2004_Breed <- d4 %>% filter(years == 1997 | years == 2004) %>%
  filter(Season == "Breed") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE)) %>%
  set_names(new_names)
  
## calculate percentage change in total pop from 1982-2003 and 2004-2023
Change_2004_2023_Breed <- d4 %>% filter(years == 2004 | years == 2023) %>%
  filter(Season == "Breed") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE)) %>%
  set_names(new_names)

Change_1997_2004_Molt <- d4 %>% filter(years == 1997 | years == 2004) %>%
  filter(Season == "Molt") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2023
Change_2004_2023_Molt <- d4 %>% filter(years == 2004 | years == 2023) %>%
  filter(Season == "Molt") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

Change_1997_2004_Pup <- d4 %>% filter(years == 1997 | years == 2004) %>%
  filter(Season == "Pup") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2023
Change_2004_2023_Pup <- d4 %>% filter(years == 2004 | years == 2023) %>%
  filter(Season == "Pup") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2023
Change_1997_2023_Breed <- d4 %>% filter(years == 1997 | years == 2023) %>%
  filter(Season == "Breed") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2023
Change_1997_2023_Molt <- d4 %>% filter(years == 1997 | years == 2023) %>%
  filter(Season == "Molt") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2023
Change_1997_2023_Pup <- d4 %>% filter(years == 1997 | years == 2023) %>%
  filter(Season == "Pup") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)


change.df <- bind_rows(Change_1997_2004_Breed, Change_2004_2023_Breed,
                       Change_1997_2004_Molt,  Change_2004_2023_Molt, 
                       Change_1997_2004_Pup,  Change_2004_2023_Pup, 
                       Change_1997_2023_Breed, Change_1997_2023_Molt, Change_1997_2023_Pup)


change.df$Season <- c("Breed", "Breed", "Molt", "Molt", "Pup", "Pup", "Breed", "Molt", "Pup")
change.df$Range <- c("1997-2004", "2004-2023", "1997-2004", "2004-2023", "1997-2004", "2004-2023","1997-2023", "1997-2023", "1997-2023" )
change.df$Duration <- c(7, 19, 
                        7, 19, 
                        7, 19, 
                        26,26, 26)
change.df

# Plot for percentage changes

change.df$Change <- (change.df$ENDPOP-change.df$STARTPOP)/change.df$STARTPOP
change.df$Change.lo <-  ((change.df$ENDPOP-change.df$ENDSE*change.df$ENDPOP)-
                           (change.df$STARTPOP-change.df$STARTSE*change.df$STARTPOP))/
                              (change.df$STARTPOP-change.df$STARTSE*change.df$STARTPOP)
change.df$Change.hi  <-  ((change.df$ENDPOP+change.df$ENDSE*change.df$ENDPOP)-
                            (change.df$STARTPOP+change.df$STARTSE*change.df$STARTPOP))/
                               (change.df$STARTPOP+change.df$STARTSE*change.df$STARTPOP)
change.df$Label <- as.character(paste(change.df$Season, change.df$Range, sep = " "))


dodge <- position_dodge(width=0.5)  

ggplot(change.df, aes(x=factor(Range, levels = c("1997-2004", "2004-2023", "1997-2023")),
                               y = Change, color = Season)) +
  geom_pointrange(aes(ymin=Change.lo, ymax=Change.hi), size = .5, position = dodge) + 
  geom_hline(yintercept=0, linetype=2) +
  xlab("Year Range") + 
  ylab("Percent change") +
  theme_classic(base_size = 18) +
  theme(legend.title=element_blank(),
        legend.position = c(0.8, 0.8))




####-----plot covariate effects

coef.data <- tidy(CIs)

#get just the coefficient terms
coef.data <- coef.data %>%
  filter(str_detect(term, "^C."))

coef.data$Season <- ifelse(str_detect(coef.data$term, "_A"), "Breeding Adults", 
                           ifelse(str_detect(coef.data$term, "_M"), "Molting", 
                              ifelse(str_detect(coef.data$term, "_P"), "Pups", "All")))

coef.data$Significant <- ifelse(coef.data$conf.up < 0, "Negative", 
                           ifelse(coef.data$conf.low > 0, "Positive", "Neutral"))

coef.data$Labels <- c("Coyote BL breed","Coyote BL molt","Coyote BL pup",
                        "Coyote DE breed", "Coyote DE molt","Coyote DE pup",
                        "Coyote DP breed","Coyote DP molt", "Coyote DP pup",
                        "MOCI breed",  
                        "MOCI molt", 
                        "MOCI pup", 
                      "MOCI AMJ lag breed",
                      "MOCI AMJ lag molt",
                      "MOCI AMJ lag pup",
                      "MOCI OND lag breed",
                      "MOCI OND lag molt",
                      "MOCI OND lag pup",
                      "eSeal DP-DE-PRH",
                      "Dist_BL",
                      "Dist_DE",
                      "Dist_DP",
                      "Dist_PRH",
                      "Dist_TB",
                      "Dist_TP")

coef.data$Labels <- c("COYOTE BL Adult","COYOTE BL molt","COYOTE BL pup",
                      "COYOTE DE Adult", "COYOTE DE molt","COYOTE DE pup",
                      "COYOTE DP Adult","COYOTE DP molt", "COYOTE DP pup",
                      "MOCI breed",  
                      "MOCI molt", 
                      "MOCI pup", 
                      "MOCI AMJ lag Adult",
                      "MOCI AMJ lag molt",
                      "MOCI AMJ lag pup",
                      "MOCI OND lag Adult",
                      "MOCI OND lag molt",
                      "MOCI OND lag pup",
                      "eSeal DP-DE-PRH",
                      "DISTURB BL",
                      "DISTURB DE",
                      "DISTURB DP",
                      "DISTURB PRH",
                      "DISTURB TB",
                      "DISTURB TP")





coef.data %>% 
ggplot(aes(fct_reorder(Labels, estimate), estimate, shape = Season, color = Significant)) +
#ggplot(aes(term, estimate, shape = Season, color = Significant)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), size = 0.8) +
  coord_flip() +
  theme(legend.position = c(0.2, 0.2)) + 
  scale_color_manual(values=c("red", "#999999", "#56B4E9")) +
  scale_shape_manual(values=c(15, 16, 17, 1)) +
  geom_hline(yintercept = 0, lty = 2) +
  #scale_x_discrete( labels = coef.data$Labels) + 
  xlab(NULL) +
  theme_sjplot2(base_size = 18)


########################################################################

## Parameter key ########
## Z = design matrix = Spatial population structure
## R = observation errors          equal
## U = growth parameter            (un)equal, or matrix [Z x 1]
## Q = hidden state process        diagonal and (un)equal
#                                  equalvarcov or unconstrained
## B = effect of column on row     unequal (these are the interactions)

#####################################################################
## some coefs
coef(BESTMODEL, type="matrix")$R  # observation errors  
coef(BESTMODEL, type="matrix")$Q  # hidden state process
(coef(BESTMODEL, type="matrix")$U)  # growth | su


B <- coef(BESTMODEL, type="matrix")$B  # effect of column on row

library("corrplot")

#Scale to get to -1 to +1
max(B)
min(B)
B <- B/abs(max(B))

#or scale with floor ceiling.

#B <- ifelse(B > 3, 3, B)
#B <- ifelse(B < -3, -3, B)

#replace o with NA
B[B==0]<-NA



dimnames(B) <- list(c("BL Adult", "BL Molt", "BL Pup",
                      "DE Adult", "DE Molt", "DE Pup",
                      "DP Adult", "DP Molt", "DP Pup",
                      #"DR Adult", "DR Molt", #"DP Pup",
                      #"DP Adult", "PB Molt", #"DP Pup",
                      "PRH Adult", "PRH Molt", "PRH Pup",
                      "TB Adult", "TB Molt", "TB Pup",
                      "TP Adult", "TP Molt", "TP Pup"),
                    c("BL Adult", "BL Molt", "BL Pup",
                      "DE Adult", "DE Molt", "DE Pup",
                      "DP Adult", "DP Molt", "DP Pup",
                      # "DR Adult", "DR Molt", #"DP Pup",
                      #"DP Adult", "PB Molt", #"DP Pup",
                      "PRH Adult", "PRH Molt", "PRH Pup",
                      "TB Adult", "TB Molt", "TB Pup",
                      "TP Adult", "TP Molt", "TP Pup"))
par(mfrow = c(1,1))
corrplot(B, method="color", is.corr = FALSE, rect.col = "gray", 
         outline = TRUE,
         addgrid.col = "gray", tl.col = "black",
         na.label = ".", na.label.col = "gray")

## Process correlation between sites --------
#recover q.matrix p 115-116
# use m.MOLT_IND_SITE_B_diag_uneq_q_unc to show all correlations for Q



Q.unc <- coef(BESTMODEL, type="matrix")$Q
h <- diag(1 / sqrt(diag(Q.unc)))
Q.corr <- h %*% Q.unc %*% h

## get Z.model from appropriate model
Z.model=factor(1:18) 

rownames(Q.corr) <- unique(Z.model)
colnames(Q.corr) <- unique(Z.model)
# Process Correlations between sites
Q.corr

### PREDICT -----------

##foreward 5 years
# must add c data
# lets say coyote at DE and DP and warm PDO

summary(t(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag))

# try 4 scenarios
#  A  Coyote increase, normal MOCI
#  A  Coyote increase, Poor MOCI
#  A  Coyote removal, normal MOCI
#  A  Coyote removal, Poor MOCI

# example A
#no coyotes, good upwelling
c_forecast_GGG=matrix(c(rep(0, times = 10),
                    rep(0, times = 10), #1,1,1,1,1,  Contrast with Y/N coyote
                    rep(0, times = 10), #1,1,1,1,1,
                    rep(0, times = 10), #1,1,1,1,1,
                    rep(0, times = 10), #1,1,1,1,1,
                    
                    rep(0, times = 10), #1,1,1,1,1,
                    -1,-1,-2,0,-2,  0,1,1,-1,-2, #MOCI
                    0,-2,1,0,2, -2,-1,0,1,-1,  #MOCI AMJ lag
                    0,2,-1,0,-2, -1,-1,0,1,-2,  #MOCI OMD lag
                    rep(0, times = 10),  #  eSeal
                    
                    rep(-1, times = 10), # dist BL
                    rep(-1, times = 10), # dist DE
                    rep(-1, times = 10), # dist DP
                    rep(-1, times = 10), # dist PRH
                    rep(-1, times = 10), # dist TB
                    
                    rep(-1, times = 10)), # dist TP
                    
                  nrow = 16, ncol = 10,
                  byrow = TRUE)


#TODO:
c_new_GGG <- cbind(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag, c_forecast_GGG)

forecast_GGG <- predict(BESTMODEL, type="ytT", n.ahead = 10, interval = "prediction", 
                     nsim = 100,
                     newdata = list(c = c_forecast_GGG))


forecast_GGG_plot_data <- forecast_GGG$pred

d.1 <- forecast_GGG_plot_data[,1:2]
d.2 <- exp(forecast_GGG_plot_data[,3:9])
forecast_GGG_plot_data <- tibble(d.1, d.2)

summary(forecast_GGG_plot_data$t) #t max = 37
#just keep the forecasted data...also a function in predict , include = 0, that might work
forecast_GGG_plot_data <- filter(forecast_GGG_plot_data, t > 37-10)

#END TODO:


# example B
# more coyotes, less upwelling
c_forecast_BBB=matrix(c(-0.8,-0.8,-0.8,-0.8,-0.8,1,1,1,1,1,  #BL
                    rep(3.02, times = 10),    #DE
                    rep(2.62, times = 10),    #DP
                    rep(0, times = 10),   #PRH
                    rep(0, times = 10),      #TB
                    
                    rep(0, times = 10),      #TP 
                    0,2,2,0,0,1,1,-2,0,2, # Warm MOCI
                    1,0,2,-1,0,0,1,-1,1,0, #MOCI AMJ lag
                    1,0,2,-1,0,0,1,-1,1,0, #MOCI OMD lag
                    rep(0, times = 10),  #  eSeal
                    
                    rep(1, times = 10), # dist BL
                    rep(1, times = 10), # dist DE
                    rep(1, times = 10), # dist DP
                    rep(1, times = 10), # dist PRH
                    rep(1, times = 10), # dist TB
                    
                    rep(1, times = 10)), # dist TP),  #eSeal), 
                  nrow = 16, ncol = 10,
                  byrow = TRUE)


c_new <- cbind(small_c_Coyote_3yr_MOCI_MOCI_Dist_lag, c_forecast_BBB)

forecast1_BBB <- predict(BESTMODEL, type="ytT", n.ahead = 10, interval = "prediction", 
                     nsim = 100,
                     newdata = list(c = c_forecast_BBB))


forecast1_BBB_plot_data <- forecast1_BBB$pred

d.1 <- forecast1_BBB_plot_data[,1:2]
d.2 <- exp(forecast1_BBB_plot_data[,3:9])
forecast1_BBB_plot_data <- tibble(d.1, d.2)

names(forecast1_BBB_plot_data)

ggplot(forecast1_BBB_plot_data, aes(t+1996, estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = `Lo 80`, ymax = `Hi 80`), alpha =0.2) + 
  geom_point(aes(t+1996, y), color = "blue4") + 
  geom_vline(xintercept = 2023.5, linetype = 2) + 
  geom_vline(xintercept = 2004, linetype = 2, color = "red3") + #show TV timepoint
  geom_line(data = forecast_GGG_plot_data, 
            aes(x = t+1996, y = estimate)) +
  geom_ribbon(data = forecast_GGG_plot_data, 
              aes(ymin = `Lo 80`, ymax = `Hi 80`), alpha =0.2, color = "blue4") + 
  xlab("Year") + 
  ylab("Seals") +
  facet_wrap(.~.rownames, ncol = 3)


# autoplot(forecast1) +
#   theme_grey(base_size = 18) +
#   #ylim(3,7.6) #+
#   scale_y_continuous(trans="exp")















