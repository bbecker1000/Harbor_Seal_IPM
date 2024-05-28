## all years adult molt PLOTS


#############################
# Model plots -------------------------------------
#############################
BESTMODEL <- m06.ut.All #m.Ind_Molt_Adult_Coyote_PDO_B_equalcov_1U_TV_Site # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U_R0 # # m.1pop_Molt_Adult_Coyote_PDO_B_custom_1U# m.1pop_Molt_Adult_Coyote_PDO_B_unc_2U #m.1pop_Coyote_PDO_Xo_fixed_B_unc_tinitx_1 #m.1pop_Coyote_PDO_Xo_fixed_B_unc# m.1pop_Coyote_PDO_B_unc #   m.1pop_Coyote_PDO_B_unc   m.5pop_Coyote_PDO_B_unc



autoplot(BESTMODEL)

autoplot(BESTMODEL, plot.type = "fitted.ytT") + # xtT
  ylim(4,8) # xtT

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
names(d)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                       "DE_Breed", "DE_Molt", #"DE_Pup",
                       "DP_Breed", "DP_Molt", #"DP_Pup",
                       "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                       "TB_Breed", "TB_Molt", #"TB_Pup",
                       "TP_Breed", "TP_Molt" #"TP_Pup",
                       # "DR_Breed", "DR_Molt", #"DR_Pup",
                       # "PB_Breed", "PB_Molt" #"PB_Pup"
)
names(d.se)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                          "DE_Breed", "DE_Molt", #"DE_Pup",
                          "DP_Breed", "DP_Molt", #"DP_Pup",
                          "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                          "TB_Breed", "TB_Molt", #"TB_Pup",
                          "TP_Breed", "TP_Molt" #"TP_Pup",
                          # "DR_Breed", "DR_Molt", #"DR_Pup",
                          #  "PB_Breed", "PB_Molt" #"PB_Pup"
)


d <- as_tibble(cbind(years,d))

d2 <- d %>% pivot_longer(cols = c(2:13), names_to = "Subsite_Season", values_to = "log_est")
d2.se <- d.se %>% pivot_longer(cols = c(1:12), names_to = "Subsite_Season", values_to = "log_est")
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
  xlim(1982, 2022) +
  ylim(-1.5, 2) + 
  theme_minimal(base_size = 20) +
  ylab("index of log abundance (N-No)") +
  xlab("Year") +
  facet_wrap(.~Season, ncol = 1)


####----plot overall pop size for each age class summing sites

d.tot <- as_tibble(t(exp(BESTMODEL$states)))  #exp to original scale
d.tot.se <- as_tibble(t(BESTMODEL$states.se)) #exp to original scale
#add header names

names(d.tot)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                           "DE_Breed", "DE_Molt", #"DE_Pup",
                           "DP_Breed", "DP_Molt", #"DP_Pup",
                           "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                           "TB_Breed", "TB_Molt", #"TB_Pup",
                           "TP_Breed", "TP_Molt" #"TP_Pup",
                           # "DR_Breed", "DR_Molt", #"DR_Pup",
                           # "PB_Breed", "PB_Molt" #"PB_Pup"
)
names(d.tot.se)[c(1:12)] <- c("BL_Breed", "BL_Molt", #"BL_Pup",
                              "DE_Breed", "DE_Molt", #"DE_Pup",
                              "DP_Breed", "DP_Molt", #"DP_Pup",
                              "PRH_Breed", "PRH_Molt",# "PRH_Pup",
                              "TB_Breed", "TB_Molt", #"TB_Pup",
                              "TP_Breed", "TP_Molt" #"TP_Pup",
                              # "DR_Breed", "DR_Molt", #"DR_Pup",
                              #  "PB_Breed", "PB_Molt" #"PB_Pup"
)

#sum all molt and breeding sites by year
Breed.tot <- d.tot %>%
  dplyr::reframe(Breed = rowSums(across(c(BL_Breed, DE_Breed, DP_Breed, PRH_Breed, TB_Breed, TP_Breed))))
Molt.tot <- d.tot %>%
  dplyr::reframe(Molt = rowSums(across(c(BL_Molt, DE_Molt, DP_Molt, PRH_Molt, TB_Molt, TP_Molt))))

Breed.tot.se <- d.tot.se %>%
  dplyr::reframe(Breed.se = rowSums(across(c(BL_Breed, DE_Breed, DP_Breed, PRH_Breed, TB_Breed, TP_Breed)))/6) #get mean
Molt.tot.se <- d.tot.se %>%
  dplyr::reframe(Molt.se = rowSums(across(c(BL_Molt, DE_Molt, DP_Molt, PRH_Molt, TB_Molt, TP_Molt)))/6) #get mean


#put in a single table
d.Breed.Molt.tot <- as_tibble(cbind(years,Breed.tot, Molt.tot))


d.Breed.Molt.tot <- d.Breed.Molt.tot %>% pivot_longer(cols = c(2:3), names_to = "Season", values_to = "estimate")

d.Breed.Molt.se <- as_tibble(cbind(years,Breed.tot.se, Molt.tot.se))



d.Breed.Molt.se <- d.Breed.Molt.se %>% pivot_longer(cols = c(2:3), names_to = "Season", values_to = "SE")
d.Breed.Molt.se <- d.Breed.Molt.se[,3] #remove season and year

d4 <- as_tibble(cbind(d.Breed.Molt.tot,d.Breed.Molt.se))


#total pop plot with ses
ggplot(d4, aes(x = years, y = estimate, color = Season)) +
  #geom_line(aes(linetype = Season), size = 1.25) +
  geom_point(size = 2) + 
  geom_line(linewidth = 1.1) +
  # CI at 90% use SE = 1.645, 80% = 1.28
  geom_ribbon(aes(ymin = estimate-1.28*SE*estimate, ymax = estimate+1.28*SE*estimate, fill = Season),  
              alpha = 0.2, color = NA) +
  #geom_hline(yintercept = c(-1,0,1), lty = 2) +
  xlim(1982, 2022) +
  ylim(0, 5500) + 
  theme_minimal(base_size = 20) +
  ylab("Estimated abundance") +
  xlab("Year") 



## change over time ---

# get data
## calculate percentage change in total pop from 1982-2003 and 2004-2022

#column names
new_names <- c("STARTPOP", "ENDPOP", "STARTSE", "ENDSE")

Change_1983_2003_Breed <- d4 %>% filter(years == 1983 | years == 2003) %>%
  filter(Season == "Breed") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE)) %>%
  set_names(new_names)
  
## calculate percentage change in total pop from 1982-2003 and 2004-2022
Change_2004_2022_Breed <- d4 %>% filter(years == 2004 | years == 2022) %>%
  filter(Season == "Breed") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE)) %>%
  set_names(new_names)

Change_1983_2003_Molt <- d4 %>% filter(years == 1983 | years == 2003) %>%
  filter(Season == "Molt") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2022
Change_2004_2022_Molt <- d4 %>% filter(years == 2004 | years == 2022) %>%
  filter(Season == "Molt") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2022
Change_1983_2022_Breed <- d4 %>% filter(years == 1983 | years == 2022) %>%
  filter(Season == "Breed") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

## calculate percentage change in total pop from 1982-2003 and 2004-2022
Change_1983_2022_Molt <- d4 %>% filter(years == 1983 | years == 2022) %>%
  filter(Season == "Molt") %>%
  pivot_wider(names_from = c(years,Season), values_from = c(estimate, SE))%>%
  set_names(new_names)

change.df <- bind_rows(Change_1983_2003_Breed, Change_2004_2022_Breed,Change_1983_2003_Molt,
          Change_2004_2022_Molt, Change_1983_2022_Breed, Change_1983_2022_Molt)


change.df$Season <- c("Breed", "Breed", "Molt", "Molt", "Breed", "Molt")
change.df$Range <- c("1983-2003", "2004-2022", "1983-2003", "2004-2022", "1983-2022", "1983-2022")
change.df$Duration <- c(20, 18, 20, 18, 39,39)
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

ggplot(change.df, aes(x=factor(Range, levels = c("1983-2003", "2004-2022", "1983-2022")),
                               y = Change, color = Season)) +
  geom_pointrange(aes(ymin=Change.lo, ymax=Change.hi), size = .3, position = dodge) + 
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
                           ifelse(str_detect(coef.data$term, "_M"), "Molting", "Pups"))

ggplot(coef.data, aes(term, estimate, color = Season, shape = Season)) +
  
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), size = 0.8) +
  coord_flip() +
  theme(legend.position = c(0.2, 0.2)) + 
  #ylim(-1.25, 0.25) +
  
  scale_x_discrete(name ="Covariate",
                   labels=c("Coyote BL breed","Coyote BL molt","Coyote DE breed", "Coyote DE molt",
                            "Coyote DP breed","Coyote DP molt", 
                            "UI breed", "UI lag breed", "UI lag molt", "UI molt"), 
                   limits=rev(levels(coef.data$term))) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_sjplot(base_size = 18)


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

exp(0.2644765)
exp(0.237)

B <- coef(BESTMODEL, type="matrix")$B  # effect of column on row

library("corrplot")

#Scale to get to -1 to +1
max(B)
min(B)
B <- B/abs(max(B))

#or scale with floor ceiling.

#B <- ifelse(B > 3, 3, B)
#B <- ifelse(B < -3, -3, B)



dimnames(B) <- list(c("BL Adult", "BL Molt", #"BL Pup",
                      "DE Adult", "DE Molt", #"DE Pup",
                      "DP Adult", "DP Molt", #"DP Pup",
                      #"DR Adult", "DR Molt", #"DP Pup",
                      #"DP Adult", "PB Molt", #"DP Pup",
                      "PRH Adult", "PRH Molt",
                      "TB Adult", "TB Molt", #"TB Pup",
                      "TP Adult", "TP Molt"), #"TP Pup"),
                    c("BL Adult", "BL Molt", #"BL Pup",
                      "DE Adult", "DE Molt", #"DE Pup",
                      "DP Adult", "DP Molt", #"DP Pup",
                      # "DR Adult", "DR Molt", #"DP Pup",
                      #"DP Adult", "PB Molt", #"DP Pup",
                      "PRH Adult", "PRH Molt",
                      "TB Adult", "TB Molt", #"TB Pup",
                      "TP Adult", "TP Molt")) #"TP Pup"),

par(mfrow = c(1,1))
corrplot(B, method="color", is.corr = FALSE)

## Process correlation between sites --------
#recover q.matrix p 115-116
# use m.MOLT_IND_SITE_B_diag_uneq_q_unc to show all correlations for Q



Q.unc <- coef(BESTMODEL, type="matrix")$Q
h <- diag(1 / sqrt(diag(Q.unc)))
Q.corr <- h %*% Q.unc %*% h

## get Z.model from appropriate model
Z.model=factor(1:12) 

rownames(Q.corr) <- unique(Z.model)
colnames(Q.corr) <- unique(Z.model)
# Process Correlations between sites
Q.corr

### PREDICT -----------

##foreward 5 years
# must add c data
# lets say coyote at DE and DP and warm PDO

summary(t(small_c_Coyote_3yr_UI_UI_lag))



# example A
#no coyotes, good upwelling
c_forecast=matrix(c(0,0,0,0,0,  0,0,0,0,0,
                    0,0,0,0,0,  0,0,0,0,0, #1,1,1,1,1,  Contrast with Y/N coyote
                    0,0,0,0,0,  0,0,0,0,0, #1,1,1,1,1,
                    0,0,0,0,0,  0,0,0,0,0, #1,1,1,1,1,
                    0,0,0,0,0,  0,0,0,0,0, #1,1,1,1,1,
                    0,0,0,0,0,  0,0,0,0,0, #1,1,1,1,1,
                    2,1,2,2,2,  2,2,3,3,0, #UI
                    0,2,1,0,2,-2,  -1,0,2,2), #UI lag
                  
                  nrow = 8, ncol = 10,
                  byrow = TRUE)

# example B
# more coyotes, less upwelling
c_forecast=matrix(c(-0.8,-0.8,-0.8,-0.8,-0.8,1,1,1,1,1,  #BL
                    2,2,2,2,2,2,2,2,2,2,  #DE
                    2,2,2,2,2,2,2,2,2,2,   #DP
                    -0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,-0.8,  #PRH
                    1,1,1,1,1,1,1,1,1,1,  #TB
                    -0.8,-0.8,-0.8,-0.8,-0.8,1,1,1,1,1,  #TP 
                    0,-2,-2,0,0,-1,-1,-1,-1,0, #UI
                    -1,0,-2,-2,0,0,-1,-1,-1,-1), #UI lag
                  nrow = 8, ncol = 10,
                  byrow = TRUE)


c_new <- cbind(small_c_Coyote_3yr_UI_UI_lag, c_forecast)

forecast1 <- predict(BESTMODEL, n.ahead = 10, interval = "prediction", 
                     nsim = 100,
                     newdata = list(c = c_forecast))


forecast1_plot_data <- forecast1$pred

d.1 <- forecast1_plot_data[,1:2]
d.2 <- exp(forecast1_plot_data[,3:9])
forecast1_plot_data <- tibble(d.1, d.2)

names(forecast1_plot_data)

ggplot(forecast1_plot_data, aes(t+1981, estimate)) +
  geom_line() +
  geom_ribbon(aes(ymin = `Lo 80`, ymax = `Hi 80`), alpha =0.2) + 
  geom_point(aes(t+1981, y), color = "blue4") + 
  geom_vline(xintercept = 2022.5, linetype = 2, color = "blue4") + 
  geom_vline(xintercept = 2003.5, linetype = 2) + #show TV timepoint
  xlab("Year") + 
  ylab("Seals") +
  facet_wrap(.~.rownames, ncol = 2)


autoplot(forecast1) +
  theme_grey(base_size = 18) +
  #ylim(3,7.6) #+
  scale_y_continuous(trans="exp")


