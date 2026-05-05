# ============================================================================
# PLOTTING CODE FOR BEST MODEL (INDEPENDENT - MODEL A)
# Updated for 2025 models with new C matrix structure
# ============================================================================

library(MARSS)
library(tidyverse)
library(ggplot2)

# Set best model
BESTMODEL <- m.A.indep

summary(BESTMODEL)

# Define years vector
years <- 1997:(1997 + ncol(dat) - 1)

# ============================================================================
# 1. BASIC MODEL DIAGNOSTICS
# ============================================================================

autoplot(BESTMODEL)

autoplot(BESTMODEL, plot.type = "fitted.ytT") +
  ylim(3, 8) +
  ylab("Estimate (log count)")

autoplot(BESTMODEL, plot.type = "xtT") +
  ylim(3.1, 7.6)

# ============================================================================
# 2. CONFIDENCE INTERVALS
# ============================================================================

CIs <- MARSSparamCIs(BESTMODEL, alpha = 0.11)
CIs

# ============================================================================
# 3. CHANGE BY SITE AND CLASS
# ============================================================================

# Extract states (change from starting values)
d <- as_tibble(t(BESTMODEL$states - BESTMODEL$states[, 1]))
d.se <- as_tibble(t(BESTMODEL$states.se))

# Column names - updated to use Adult instead of Breed
state_names <- c(
  "BL_Adult", "BL_Molt", "BL_Pup",
  "DE_Adult", "DE_Molt", "DE_Pup",
  "DP_Adult", "DP_Molt", "DP_Pup",
  "PRH_Adult", "PRH_Molt", "PRH_Pup",
  "TB_Adult", "TB_Molt", "TB_Pup",
  "TP_Adult", "TP_Molt", "TP_Pup"
)

names(d) <- state_names
names(d.se) <- state_names

# Add years and reshape
d <- cbind(years, d)

d_long <- d %>%
  pivot_longer(cols = -years, names_to = "Site_Class", values_to = "log_est") %>%
  separate(Site_Class, into = c("Site", "Class"), sep = "_", remove = FALSE)

d_se_long <- d.se %>%
  pivot_longer(cols = everything(), names_to = "Site_Class", values_to = "log_se")

d_plot <- bind_cols(d_long, d_se_long %>% select(log_se))

# Plot log abundance change
ggplot(d_plot, aes(x = years, y = log_est, color = Site, group = Site)) +
  geom_ribbon(aes(ymin = log_est - log_se, ymax = log_est + log_se, fill = Site),
              alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(-1, 0, 1), lty = 2) +
  xlim(min(years), max(years)) +
  ylim(-2, 2) +
  theme_minimal(base_size = 20) +
  ylab("Index of log abundance (N - N₀)") +
  xlab("Year") +
  facet_wrap(~Class, ncol = 1)

ggsave("Output/Plots/logAbundance.jpeg", width = 20, height = 30, units = "cm")



# ============================================================================
# 3. CHANGE BY SITE AND CLASS - PERCENTAGE CHANGE FROM YEAR 0
# ============================================================================

# Extract states (change from starting values on log scale)
d <- as_tibble(t(BESTMODEL$states - BESTMODEL$states[, 1]))
d.se <- as_tibble(t(BESTMODEL$states.se))

# Column names
state_names <- c(
  "BL_Adult", "BL_Molt", "BL_Pup",
  "DE_Adult", "DE_Molt", "DE_Pup",
  "DP_Adult", "DP_Molt", "DP_Pup",
  "PRH_Adult", "PRH_Molt", "PRH_Pup",
  "TB_Adult", "TB_Molt", "TB_Pup",
  "TP_Adult", "TP_Molt", "TP_Pup"
)

names(d) <- state_names
names(d.se) <- state_names

# Add years and reshape
d <- cbind(years, d)

d_long <- d %>%
  pivot_longer(cols = -years, names_to = "Site_Class", values_to = "log_diff") %>%
  separate(Site_Class, into = c("Site", "Class"), sep = "_", remove = FALSE)

d_se_long <- d.se %>%
  pivot_longer(cols = everything(), names_to = "Site_Class", values_to = "log_se")

d_plot <- bind_cols(d_long, d_se_long %>% select(log_se))

# Convert log difference to percentage change
# If log(N_t) - log(N_0) = x, then N_t/N_0 = exp(x)
# Percentage change = (exp(x) - 1) * 100
d_plot <- d_plot %>%
  mutate(
    pct_change = (exp(log_diff) - 1) * 100,
    pct_change_lo = (exp(log_diff - log_se) - 1) * 100,
    pct_change_hi = (exp(log_diff + log_se) - 1) * 100
  )

# Plot percentage change free scales
ggplot(d_plot, aes(x = years, y = pct_change, color = Site, group = Site)) +
  geom_ribbon(aes(ymin = pct_change_lo, ymax = pct_change_hi, fill = Site),
              alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(-50, 0, 100, 200), lty = 2) +
  xlim(min(years), max(years)) +
  theme_minimal(base_size = 20) +
  ylab("% Change from 1997") +
  xlab("Year") +
  facet_wrap(~Class, ncol = 1, scales = "free_y")

ggsave("Output/Plots/PctChangeAbundance.jpeg", width = 20, height = 30, units = "cm")

# Alternative with fixed y-axis for comparison across classes
ggplot(d_plot, aes(x = years, y = pct_change, color = Site, group = Site)) +
  geom_ribbon(aes(ymin = pct_change_lo, ymax = pct_change_hi, fill = Site),
              alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = c(-50, 0, 100, 200, 300), lty = 2) +
  xlim(min(years), max(years)) +
  ylim(-75, 400) +
  theme_minimal(base_size = 20) +
  ylab("% Change from 1997") +
  xlab("Year") +
  facet_wrap(~Class, ncol = 1)

ggsave("Output/Plots/PctChangeAbundance_fixed.jpeg", width = 20, height = 30, units = "cm")








# ============================================================================
# 4. PUP:ADULT RATIOS
# ============================================================================

# Extract raw states (not differenced)
d_raw <- as_tibble(t(BESTMODEL$states))
d_raw_se <- as_tibble(t(BESTMODEL$states.se))

names(d_raw) <- state_names
names(d_raw_se) <- state_names

d_raw <- cbind(years, d_raw)

d_raw_long <- d_raw %>%
  pivot_longer(cols = -years, names_to = "Site_Class", values_to = "log_est") %>%
  separate(Site_Class, into = c("Site", "Class"), sep = "_", remove = FALSE)

d_raw_se_long <- d_raw_se %>%
  pivot_longer(cols = everything(), names_to = "Site_Class", values_to = "log_se")

d_raw_plot <- bind_cols(d_raw_long, d_raw_se_long %>% select(log_se))

# Filter to Adults and Pups only
d_ratio <- d_raw_plot %>%
  filter(Class != "Molt") %>%
  select(-Site_Class) %>%
  pivot_wider(names_from = Class, values_from = c(log_est, log_se))

# Calculate ratios with CIs (using simplified SE approach)
d_ratio <- d_ratio %>%
  mutate(
    ratio = log_est_Pup / log_est_Adult,
    avg_se = (log_se_Pup + log_se_Adult) / 2,
    ratio_lo = ratio - 1.68 * avg_se * ratio,
    ratio_hi = ratio + 1.68 * avg_se * ratio
  )

mean_ratio <- mean(d_ratio$ratio, na.rm = TRUE)

ggplot(d_ratio, aes(x = years, y = ratio, color = Site, group = Site)) +
  geom_ribbon(aes(ymin = ratio_lo, ymax = ratio_hi, fill = Site),
              alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = mean_ratio, lty = 2) +
  xlim(min(years), max(years)) +
  theme_grey(base_size = 20) +
  ylab("Pups:Adults") +
  xlab("Year") +
  theme(legend.position = "none") +
  facet_wrap(~Site, ncol = 3)

ggsave("Output/Plots/PupAdultRatio.jpeg", width = 30, height = 20, units = "cm")



# ============================================================================
# 4B TEST ESTUARY VS OUTER COAST PUP:ADULT RATIOS
# ============================================================================

library(lme4)
library(glmmTMB)  # Better for beta regression

# Add location type to the ratio data
d_ratio <- d_ratio %>%
  mutate(
    Location = case_when(
      Site %in% c("BL", "TB", "DE") ~ "Estuary",
      Site %in% c("PRH", "TP", "DP") ~ "Outer Coast",
      TRUE ~ NA_character_
    ),
    Location = factor(Location)
  )

# Check the data
d_ratio %>%
  group_by(Location) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    sd_ratio = sd(ratio, na.rm = TRUE),
    n = n()
  )

# Beta regression requires values strictly between 0 and 1
# Your ratios look OK (0.6-1.1), but some might exceed 1
# Check range
range(d_ratio$ratio, na.rm = TRUE)

# If any ratios >= 1, need to transform or use different approach
# Option 1: Use Gaussian if ratios are well-behaved
# Option 2: Transform to keep in (0,1)

# ----------------------------------------------------------------------------
# MODEL 1: Gaussian mixed model (simpler, often fine for ratios near 0.8)
# ----------------------------------------------------------------------------

m1_gaussian <- lmer(ratio ~ Location + (1|years) + (1|Site), 
                    data = d_ratio)
summary(m1_gaussian)

# Check assumptions
plot(m1_gaussian)
qqnorm(resid(m1_gaussian))
qqline(resid(m1_gaussian))

# ----------------------------------------------------------------------------
# MODEL 2: Beta regression (if ratios strictly in 0-1)
# ----------------------------------------------------------------------------

# Transform ratios > 1 slightly if needed
d_ratio_beta <- d_ratio %>%
  filter(!is.na(ratio)) %>%
  mutate(
    # Squeeze into (0,1) if any exceed bounds
    ratio_beta = case_when(
      ratio >= 1 ~ 0.999,
      ratio <= 0 ~ 0.001,
      TRUE ~ ratio
    )
  )

# Check transformation
range(d_ratio_beta$ratio_beta)

# Fit beta regression with glmmTMB
m2_beta <- glmmTMB(ratio_beta ~ Location + (1|years) + (1|Site),
                   family = beta_family(link = "logit"),
                   data = d_ratio_beta)
summary(m2_beta)

# ----------------------------------------------------------------------------
# MODEL COMPARISON AND INFERENCE
# ----------------------------------------------------------------------------

# Test significance of Location effect
# For Gaussian model:
library(lmerTest)  # Adds p-values to lmer
m1_gaussian <- lmer(ratio ~ Location + (1|years) + (1|Site), 
                    data = d_ratio)
summary(m1_gaussian)
anova(m1_gaussian)

# For beta model - likelihood ratio test
m2_null <- glmmTMB(ratio_beta ~ 1 + (1|years) + (1|Site),
                   family = beta_family(link = "logit"),
                   data = d_ratio_beta)
anova(m2_null, m2_beta)

# ----------------------------------------------------------------------------
# VISUALIZE DIFFERENCES
# ----------------------------------------------------------------------------

# Summary by location
location_summary <- d_ratio %>%
  group_by(Location, Site) %>%
  summarise(
    mean_ratio = mean(ratio, na.rm = TRUE),
    se_ratio = sd(ratio, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

ggplot(d_ratio, aes(x = Location, y = ratio, fill = Location)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Site), width = 0.2, alpha = 0.5) +
  geom_hline(yintercept = mean_ratio, lty = 2) +
  theme_classic(base_size = 18) +
  ylab("Pup:Adult Ratio") +
  xlab("Location Type") +
  scale_fill_manual(values = c("Estuary" = "skyblue", "Outer Coast" = "coral"))

ggsave("Output/Plots/PupAdultRatio_by_Location.jpeg", width = 20, height = 15, units = "cm")

# Time series by location type
ggplot(d_ratio, aes(x = years, y = ratio, color = Location, group = Site)) +
  geom_line(alpha = 0.5) +
  geom_smooth(aes(group = Location), method = "loess", se = TRUE) +
  geom_hline(yintercept = mean_ratio, lty = 2) +
  theme_minimal(base_size = 18) +
  ylab("Pup:Adult Ratio") +
  xlab("Year") +
  scale_color_manual(values = c("Estuary" = "blue", "Outer Coast" = "red"))

ggsave("Output/Plots/PupAdultRatio_TimeSeries_by_Location.jpeg", width = 25, height = 15, units = "cm")


# Time series with site colors grouped by location
ggplot(d_ratio, aes(x = years, y = ratio, group = Site)) +
  geom_line(aes(color = Site), alpha = 0.7, linewidth = 0.8) +
  geom_smooth(aes(group = Location, linetype = Location), 
              method = "loess", se = TRUE, color = "blue", fill = "gray50", alpha = 0.3) +
  geom_hline(yintercept = mean_ratio, lty = 2, alpha = 0.5) +
  theme_minimal(base_size = 18) +
  ylab("Pup:Adult Ratio") +
  xlab("Year") +
  scale_color_manual(values = c("BL" = "#E41A1C", "DE" = "#FF7F00", "TB" = "#984EA3",
                                "DP" = "#377EB8", "PRH" = "#4DAF4A", "TP" = "#00CED1")) +
  guides(color = guide_legend(title = "Site"),
         linetype = guide_legend(title = "Location"))

ggsave("Output/Plots/PupAdultRatio_TimeSeries_by_Location.jpeg", width = 25, height = 15, units = "cm")

# Time series by location type with unique colors per site
ggplot(d_ratio, aes(x = years, y = ratio, color = Site, group = Site)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +
  geom_smooth(aes(group = Location, color = NULL, fill = Location), 
              method = "loess", se = TRUE, linewidth = 1.5) +
  geom_hline(yintercept = mean_ratio, lty = 2) +
  theme_minimal(base_size = 18) +
  ylab("Pup:Adult Ratio") +
  xlab("Year") +
  scale_fill_manual(values = c("Estuary" = "blue", "Outer Coast" = "red"))

ggsave("Output/Plots/PupAdultRatio_TimeSeries_by_Location.jpeg", width = 25, height = 15, units = "cm")



# ============================================================================
# 5. TOTAL POPULATION BY CLASS
# ============================================================================

d_tot <- as_tibble(t(exp(BESTMODEL$states)))
d_tot_se <- as_tibble(t(BESTMODEL$states.se))

names(d_tot) <- state_names
names(d_tot_se) <- state_names

# Sum across sites by class
totals <- tibble(
  years = years,
  Adult = rowSums(d_tot[, grep("Adult", state_names)]),
  Molt = rowSums(d_tot[, grep("Molt", state_names)]),
  Pup = rowSums(d_tot[, grep("Pup", state_names)])
)

totals_se <- tibble(
  Adult_se = rowMeans(d_tot_se[, grep("Adult", state_names)]),
  Molt_se = rowMeans(d_tot_se[, grep("Molt", state_names)]),
  Pup_se = rowMeans(d_tot_se[, grep("Pup", state_names)])
)

totals_long <- totals %>%
  pivot_longer(cols = -years, names_to = "Class", values_to = "estimate")

totals_se_long <- totals_se %>%
  mutate(years = years) %>%
  pivot_longer(cols = -years, names_to = "Class", values_to = "SE") %>%
  mutate(Class = gsub("_se", "", Class))

totals_plot <- left_join(totals_long, totals_se_long, by = c("years", "Class"))

ggplot(totals_plot, aes(x = years, y = estimate, color = Class, group = Class)) +
  geom_ribbon(aes(ymin = estimate - 1.28 * SE * estimate,
                  ymax = estimate + 1.28 * SE * estimate,
                  fill = Class),
              alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.1) +
  xlim(min(years), max(years)) +
  ylim(0, 5500) +
  theme_minimal(base_size = 20) +
  ylab("Estimated abundance") +
  xlab("Year") +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank())

ggsave("Output/Plots/TotalPop.jpeg", width = 22, height = 15, units = "cm")

# ============================================================================
# 6. POPULATION CHANGE RATES
# ============================================================================

# Get data for key years
key_years <- c(1997, 2004, max(years))

change_data <- totals_plot %>%
  filter(years %in% key_years) %>%
  select(years, Class, estimate, SE)

# Function to calculate change
calc_period_change <- function(data, start_yr, end_yr) {
  start_data <- data %>% filter(years == start_yr) %>% select(Class, estimate, SE)
  end_data <- data %>% filter(years == end_yr) %>% select(Class, estimate, SE)
  
  names(start_data) <- c("Class", "STARTPOP", "STARTSE")
  names(end_data) <- c("Class", "ENDPOP", "ENDSE")
  
  left_join(start_data, end_data, by = "Class") %>%
    mutate(
      Range = paste0(start_yr, "-", end_yr),
      Duration = end_yr - start_yr,
      Change = (ENDPOP - STARTPOP) / STARTPOP,
      Change_lo = ((ENDPOP - ENDSE * ENDPOP) - (STARTPOP + STARTSE * STARTPOP)) /
        (STARTPOP + STARTSE * STARTPOP),
      Change_hi = ((ENDPOP + ENDSE * ENDPOP) - (STARTPOP - STARTSE * STARTPOP)) /
        (STARTPOP - STARTSE * STARTPOP)
    )
}

change_df <- bind_rows(
  calc_period_change(change_data, 1997, 2004),
  calc_period_change(change_data, 2004, max(years)),
  calc_period_change(change_data, 1997, max(years))
)

ggplot(change_df, aes(x = factor(Range, levels = c("1997-2004", 
                                                   paste0("2004-", max(years)), 
                                                   paste0("1997-", max(years)))),
                      y = Change, shape = Class, color = Class)) +
  geom_pointrange(aes(ymin = Change_lo, ymax = Change_hi),
                  size = 0.5, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = 2) +
  xlab("Year Range") +
  ylab("Change") +
  theme_classic(base_size = 18) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8),
        legend.background = element_rect(linewidth = 0.5, linetype = "solid", color = "black"))

ggsave("Output/Plots/PopChangeRate.jpeg", width = 22, height = 15, units = "cm")

# ============================================================================
# 7. COVARIATE EFFECTS
# ============================================================================

# ============================================================================
# 7. COVARIATE EFFECTS
# ============================================================================
# ============================================================================
# 7. COVARIATE EFFECTS
# ============================================================================

# ============================================================================
# 7. COVARIATE EFFECTS
# ============================================================================

coef_data <- tidy(CIs) %>%
  filter(str_detect(term, "^C\\."))

coef_data <- coef_data %>%
  mutate(
    Class = case_when(
      str_detect(term, "_A$|_A,") ~ "Adult",
      str_detect(term, "_M$|_M,") ~ "Molting",
      str_detect(term, "_P$|_P,") ~ "Pup",
      TRUE ~ "All"
    ),
    Significant = case_when(
      conf.up < 0 ~ "Negative",
      conf.low > 0 ~ "Positive",
      TRUE ~ "Neutral"
    ),
    # Create readable labels from term names
    Label = term %>%
      str_remove("^C\\.") %>%
      str_replace_all("_", " ") %>%
      # Apply label replacements
      str_replace_all("AMJ", "Spring") %>%
      str_replace_all("JFM", "Winter") %>%
      str_replace_all("OND", "Prior Fall") %>%
      str_replace_all("\\bCoy\\b", "Coyote") %>%
      # Handle new ES-site format
      str_replace_all("ES-DE", "Elephant Seals DE") %>%
      str_replace_all("ES-PRH", "Elephant Seals PRH") %>%
      str_replace_all("\\bES\\b", "Elephant Seals") %>%
      str_replace_all("\\bDist\\b", "Disturbance") %>%
      # Replace class suffixes
      str_replace_all("\\bP$", "(Pup)") %>%
      str_replace_all("\\bA$", "(Adult)") %>%
      str_replace_all("\\bM$", "(Molting)") %>%
      str_replace_all(" P ", " (Pup) ") %>%
      str_replace_all(" A ", " (Adult) ") %>%
      str_replace_all(" M ", " (Molting) ") %>%
      # Title case FIRST
      str_to_title() %>%
      # THEN fix site codes to uppercase
      str_replace_all("\\bBl\\b", "BL") %>%
      str_replace_all("\\bDe\\b", "DE") %>%
      str_replace_all("\\bDp\\b", "DP") %>%
      str_replace_all("\\bPrh\\b", "PRH") %>%
      str_replace_all("\\bTb\\b", "TB") %>%
      str_replace_all("\\bTp\\b", "TP") %>%
      str_replace_all("\\bMoci\\b", "MOCI")
  ) %>%
  # Sort by Class then alphabetically within class (DESCENDING for reverse order)
  arrange(factor(Class, levels = c("Adult", "Molting", "Pup", "All")), desc(Label)) %>%
  mutate(Label = fct_inorder(Label))

ggplot(coef_data, aes(x = Label, y = estimate,
                      shape = Class, color = Significant)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), size = 0.8) +
  coord_flip() +
  scale_color_manual(values = c("Negative" = "red", "Neutral" = "#999999", "Positive" = "#56B4E9")) +
  scale_shape_manual(values = c(15, 16, 17, 1)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab(NULL) +
  ylab("Coefficient estimate") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.85))

ggsave("Output/Plots/Covariates.jpeg", width = 18, height = 30, units = "cm")


#-----------
#OR Sort alphabetical#
coef_data <- coef_data %>%
  mutate(
    Class = case_when(
      str_detect(term, "_A$|_A,") ~ "Adult",
      str_detect(term, "_M$|_M,") ~ "Molting",
      str_detect(term, "_P$|_P,") ~ "Pup",
      TRUE ~ "All"
    ),
    Significant = case_when(
      conf.up < 0 ~ "Negative",
      conf.low > 0 ~ "Positive",
      TRUE ~ "Neutral"
    ),
    # Create readable labels from term names
    Label = term %>%
      str_remove("^C\\.") %>%
      str_replace_all("_", " ") %>%
      # Apply label replacements
      str_replace_all("AMJ", "Spring") %>%
      str_replace_all("JFM", "Winter") %>%
      str_replace_all("OND", "Prior Fall") %>%
      str_replace_all("\\bCoy\\b", "Coyote") %>%
      str_replace_all("\\bES\\b", "Elephant Seals") %>%
      str_replace_all("\\bDist\\b", "Disturbance") %>%
      # Replace class suffixes
      str_replace_all("\\bP$", "Pup") %>%
      str_replace_all("\\bA$", "Adult") %>%
      str_replace_all("\\bM$", "Molting") %>%
      str_replace_all(" P ", " Pup ") %>%
      str_replace_all(" A ", " Adult ") %>%
      str_replace_all(" M ", " Molting ") %>%
      # Title case FIRST
      str_to_title() %>%
      # Add parentheses around class names
      str_replace_all("\\bPup\\b", "(Pup)") %>%
      str_replace_all("\\bMolting\\b", "(Molting)") %>%
      str_replace_all("\\bAdult\\b", "(Adult)") %>%
      # THEN fix site codes to uppercase
      str_replace_all("\\bBl\\b", "BL") %>%
      str_replace_all("\\bDe\\b", "DE") %>%
      str_replace_all("\\bDp\\b", "DP") %>%
      str_replace_all("\\bPrh\\b", "PRH") %>%
      str_replace_all("\\bTb\\b", "TB") %>%
      str_replace_all("\\bTp\\b", "TP") %>%
      str_replace_all("\\bMoci\\b", "MOCI")
  ) %>%
  # Sort by Class then alphabetically within class
  arrange(factor(Class, levels = c("Adult", "Molting", "Pup", "All")), Label) %>%
  mutate(Label = fct_inorder(Label))

ggplot(coef_data, aes(x = Label, y = estimate,
                      shape = Class, color = Significant)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.up), size = 0.8) +
  coord_flip() +
  scale_color_manual(values = c("Negative" = "red", "Neutral" = "#999999", "Positive" = "#56B4E9")) +
  scale_shape_manual(values = c(15, 16, 17, 1)) +
  geom_hline(yintercept = 0, lty = 2) +
  xlab(NULL) +
  ylab("Coefficient estimate") +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_blank(),
        legend.position = c(0.85, 0.5))

ggsave("Output/Plots/Covariates.jpeg", width = 18, height = 30, units = "cm")





# ============================================================================
# 8. MODEL COEFFICIENTS
# ============================================================================

cat("\n--- Model Coefficients ---\n")
cat("R (observation error):\n")
print(coef(BESTMODEL, type = "matrix")$R[1, 1])

cat("\nQ (process error):\n")
print(coef(BESTMODEL, type = "matrix")$Q[1, 1])

cat("\nU (growth rates):\n")
print(coef(BESTMODEL, type = "matrix")$U)

# ============================================================================
# 9. TEN-YEAR FORECAST
# ============================================================================

# ============================================================================
# FORECAST WITH SCALED COVARIATES
# ============================================================================

n_forecast <- 10

# Good scenario: all covariates at -1 SD (favorable conditions)
c_forecast_good <- matrix(
  c(
    rep(-1, n_forecast),   # MOCI_JFM: cool
    rep(-1, n_forecast),   # MOCI_AMJ: cool
    rep(-1, n_forecast),   # MOCI_OND: cool
    rep(-0.5, n_forecast), # Dist_BL: low
    rep(-0.5, n_forecast), # Dist_DE: low
    rep(-0.5, n_forecast), # Dist_DP: low
    rep(0, n_forecast),    # Dist_PRH: mean (no variance)
    rep(-0.5, n_forecast), # Dist_TB: low
    rep(-0.5, n_forecast), # Dist_TP: low
    rep(-0.5, n_forecast), # Coyote_BL: low
    rep(-0.5, n_forecast), # Coyote_DE: low
    rep(-0.5, n_forecast), # Coyote_DP: low
    rep(0, n_forecast),    # Coyote_PRH: none
    rep(0, n_forecast),    # Coyote_TB: none
    rep(0, n_forecast),    # Coyote_TP: none
    rep(-0.5, n_forecast)  # eSeal: low
  ),
  nrow = 16, ncol = n_forecast, byrow = TRUE
)

# Poor scenario: all covariates at +1 SD (unfavorable conditions)
c_forecast_poor <- matrix(
  c(
    rep(1.5, n_forecast),  # MOCI_JFM: warm
    rep(1.5, n_forecast),  # MOCI_AMJ: warm
    rep(1.5, n_forecast),  # MOCI_OND: warm
    rep(1, n_forecast),    # Dist_BL: high
    rep(1, n_forecast),    # Dist_DE: high
    rep(1, n_forecast),    # Dist_DP: high
    rep(0, n_forecast),    # Dist_PRH: mean
    rep(1, n_forecast),    # Dist_TB: high
    rep(1, n_forecast),    # Dist_TP: high
    rep(1.5, n_forecast),  # Coyote_BL: high
    rep(1.5, n_forecast),  # Coyote_DE: high
    rep(1.5, n_forecast),  # Coyote_DP: high
    rep(0, n_forecast),    # Coyote_PRH: none
    rep(0, n_forecast),    # Coyote_TB: none
    rep(0, n_forecast),    # Coyote_TP: none
    rep(1, n_forecast)     # eSeal: high
  ),
  nrow = 16, ncol = n_forecast, byrow = TRUE
)

# Run forecasts
forecast_good <- predict(BESTMODEL, 
                         type = "ytT", 
                         n.ahead = n_forecast, 
                         interval = "prediction",
                         nsim = 100,
                         newdata = list(c = c_forecast_good))

forecast_poor <- predict(BESTMODEL, 
                         type = "ytT", 
                         n.ahead = n_forecast, 
                         interval = "prediction",
                         nsim = 100,
                         newdata = list(c = c_forecast_poor))

# ----------------------------------------------------------------------------
# RUN FORECASTS
# ----------------------------------------------------------------------------

forecast_good <- predict(BESTMODEL, 
                         type = "ytT", 
                         n.ahead = n_forecast, 
                         interval = "prediction",
                         nsim = 100,
                         newdata = list(c = c_forecast_good))

forecast_poor <- predict(BESTMODEL, 
                         type = "ytT", 
                         n.ahead = n_forecast, 
                         interval = "prediction",
                         nsim = 100,
                         newdata = list(c = c_forecast_poor))

# ----------------------------------------------------------------------------
# PREPARE DATA FOR PLOTTING
# ----------------------------------------------------------------------------

prep_forecast <- function(forecast_obj, scenario_name) {
  pred <- forecast_obj$pred
  
  tibble(
    t = pred$t,
    Site_Class = pred$.rownames,
    y = pred$y,
    estimate = exp(pred$estimate),
    Lo_80 = exp(pred$`Lo 80`),
    Hi_80 = exp(pred$`Hi 80`),
    Scenario = scenario_name
  ) %>%
    mutate(
      Year = t + min(years) - 1,
      # Fix: split on underscore properly
      Site = sub("_.*", "", Site_Class),
      Class = sub(".*_", "", Site_Class),
      # Standardize class names
      Class = case_when(
        Class == "ADULT" ~ "Adult",
        Class == "MOLTING" ~ "Molt",
        Class == "PUP" ~ "Pup",
        TRUE ~ Class
      )
    )
}

# Recreate forecast data
forecast_good_df <- prep_forecast(forecast_good, "Good")
forecast_poor_df <- prep_forecast(forecast_poor, "Poor")

# Verify
forecast_good_df %>% 
  filter(Year > max(years)) %>%
  select(Site, Class, Year, estimate) %>%
  distinct() %>%
  print(n = 30)

# Combine for plotting
forecast_combined <- bind_rows(forecast_good_df, forecast_poor_df)

# Breakpoint year
breakpoint_year <- years[ceiling(length(years) / 4)]

# ----------------------------------------------------------------------------
# PLOT: FACETED BY SITE AND CLASS
# ----------------------------------------------------------------------------

ggplot() +
  # Good scenario - historical + forecast
  geom_line(data = forecast_good_df, 
            aes(x = Year, y = estimate), color = "blue4", linewidth = 0.8) +
  geom_ribbon(data = filter(forecast_good_df, Year > max(years)),
              aes(x = Year, ymin = Lo_80, ymax = Hi_80),
              alpha = 0.2, fill = "blue4") +
  # Poor scenario - forecast only (line diverges after max(years))
  geom_line(data = filter(forecast_poor_df, Year > max(years)), 
            aes(x = Year, y = estimate), color = "red3", linewidth = 0.8) +
  geom_ribbon(data = filter(forecast_poor_df, Year > max(years)),
              aes(x = Year, ymin = Lo_80, ymax = Hi_80),
              alpha = 0.2, fill = "red3") +
  # Observed data points
  geom_point(data = forecast_good_df, 
             aes(x = Year, y = exp(y)), color = "black", size = 1) +
  # Reference lines
  geom_vline(xintercept = max(years) + 0.5, linetype = 3, linewidth = 0.5) +
  geom_vline(xintercept = breakpoint_year, linetype = 3, color = "gray50") +
  # Labels
  xlab("Year") +
  ylab("Estimated Abundance") +
  ylim(0, 2500) +
  theme_minimal(base_size = 16) +
  theme(panel.spacing = unit(1.5, "lines")) +
  facet_grid(Site ~ Class)

ggsave("Output/Plots/Ten_Year_Predictions.jpeg", width = 28, height = 40, units = "cm")

# ----------------------------------------------------------------------------
# PLOT: TOTAL POPULATION BY CLASS (SUMMED ACROSS SITES)
# ----------------------------------------------------------------------------

# Sum across sites for each scenario
totals_forecast <- forecast_combined %>%
  group_by(Year, Class, Scenario) %>%
  summarise(
    Total = sum(estimate),
    Total_Lo = sum(Lo_80),
    Total_Hi = sum(Hi_80),
    .groups = "drop"
  )

ggplot(totals_forecast, aes(x = Year, y = Total, color = Scenario, fill = Scenario)) +
  geom_line(linewidth = 1) +
  geom_ribbon(data = filter(totals_forecast, Year > max(years)),
              aes(ymin = Total_Lo, ymax = Total_Hi),
              alpha = 0.2, color = NA) +
  geom_vline(xintercept = max(years) + 0.5, linetype = 3) +
  geom_vline(xintercept = breakpoint_year, linetype = 3, color = "gray50") +
  scale_color_manual(values = c("Good" = "blue4", "Poor" = "red3")) +
  scale_fill_manual(values = c("Good" = "blue4", "Poor" = "red3")) +
  xlab("Year") +
  ylab("Total Estimated Abundance") +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.85, 0.85),
        legend.title = element_blank()) +
  facet_wrap(~Class, ncol = 1, scales = "free_y")

ggsave("Output/Plots/Ten_Year_Totals_by_Class.jpeg", width = 22, height = 30, units = "cm")

# ============================================================================
# 10. SAVE RESULTS
# ============================================================================

save(m.A.indep, file = "Output/m.A.indep.RData")
save(m.B.estuary, file = "Output/m.B.estuary.RData")
save(m.C.one, file = "Output/m.C.one.RData")
save(m.D.molt, file = "Output/m.D.molt.RData")

# Save AIC comparison
write.csv(df_aic, "Output/model_comparison_AIC.csv", row.names = FALSE)
