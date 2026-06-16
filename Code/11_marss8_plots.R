# ============================================================================
# 11_marss8_plots.R
# ----------------------------------------------------------------------------
# 8-SITE MARSS — PLOTS (Appendix A). Style matches IPM v3.2 (89% CI, palette).
#
# Prereqs (either):
#   source("Code/08_data_prep_8site.R"); source("Code/09_covariates_8site.R")
#   source("Code/10_marss8_models.R")
#   OR load saved objects (auto-loaded below if missing):
#   load("Output/m.A_8site.RData"); load("Output/CIs_8site.RData")
# ============================================================================

library(MARSS)
library(tidyverse)
library(patchwork)
library(broom)   # tidy() for marssMLE / marssParamCIs

dir.create("Output/Plots", showWarnings = FALSE)

# ── Constants matching IPM v3.2 ─────────────────────────────────────────────
CI_LO <- 0.055; CI_HI <- 0.945; CI_LABEL <- "89% CI"

SITE_TYPE <- c(BL="Breeding", DE="Breeding", DP="Breeding",
               DR="Haul-out only", PB="Haul-out only",
               PRH="Breeding", TB="Breeding", TP="Breeding")
SITE_COLS <- c(BL="#E41A1C", DE="#FF7F00", DP="#4DAF4A",
               DR="#984EA3", PB="#00CED1",
               PRH="#377EB8", TB="#A65628", TP="#F781BF")
SITE_LTY  <- c(BL=1, DE=1, DP=1, DR=2, PB=2, PRH=1, TB=1, TP=1)

theme_seal <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.grid.major = element_line(colour="grey88", linewidth=0.4),
      panel.border     = element_rect(colour="grey70", fill=NA, linewidth=0.5),
      axis.title       = element_text(size=rel(0.95)),
      legend.position  = "bottom",
      strip.text       = element_text(size=rel(0.90), face="bold"),
      strip.background = element_rect(fill="grey94", colour="grey80"),
      plot.title       = element_text(size=rel(1.05), face="bold"),
      plot.caption     = element_text(size=rel(0.78), colour="grey50", hjust=1),
      plot.margin      = margin(10,14,10,10))
}

# ── Load model + CIs if not in memory ───────────────────────────────────────
if (!exists("m.A_8site"))  { load("Output/m.A_8site.RData");  cat("Loaded m.A_8site.RData\n") }
if (!exists("CIs_8site"))  { load("Output/CIs_8site.RData");  cat("Loaded CIs_8site.RData\n") }
if (!exists("dat_8site"))  stop("dat_8site not found — run 08_data_prep_8site.R.")
if (!exists("years_8site")) stop("years_8site not found — run 08_data_prep_8site.R.")

BESTMODEL <- m.A_8site
years     <- years_8site
breakpoint_year <- 2004

state_names_8site <- as.vector(t(outer(
  c("BL","DE","DP","DR","PB","PRH","TB","TP"),
  c("Adult","Molt","Pup"), paste, sep = "_")))

# ── Helper: long-form smoothed state data ───────────────────────────────────
build_state_df <- function(model, state_names, yrs) {
  d    <- as_tibble(t(model$states))
  d.se <- as_tibble(t(model$states.se))
  names(d) <- names(d.se) <- state_names
  d$Year <- yrs
  d_long <- d %>%
    pivot_longer(-Year, names_to="Site_Class", values_to="log_est") %>%
    separate(Site_Class, into=c("Site","Class"), sep="_", remove=FALSE)
  d_se_long <- d.se %>%
    pivot_longer(everything(), names_to="Site_Class", values_to="log_se")
  bind_cols(d_long, d_se_long %>% select(log_se)) %>%
    mutate(SiteType = SITE_TYPE[Site],
           lo89 = log_est + qnorm(CI_LO)*log_se,
           hi89 = log_est + qnorm(CI_HI)*log_se)
}
d_plot <- build_state_df(BESTMODEL, state_names_8site, years)

# ── PLOT 1: log-abundance index relative to 1997 ────────────────────────────
d_diff <- d_plot %>%
  group_by(Site_Class) %>%
  mutate(log_diff = log_est - first(log_est),
         lo89_d   = log_diff + qnorm(CI_LO)*log_se,
         hi89_d   = log_diff + qnorm(CI_HI)*log_se) %>%
  ungroup()

p1 <- ggplot(d_diff, aes(Year, log_diff, colour=Site, linetype=Site, group=Site)) +
  geom_hline(yintercept=0, linetype=2, colour="grey40") +
  geom_ribbon(aes(ymin=lo89_d, ymax=hi89_d, fill=Site), alpha=0.12, colour=NA) +
  geom_line(linewidth=1) +
  scale_colour_manual(values=SITE_COLS) + scale_fill_manual(values=SITE_COLS) +
  scale_linetype_manual(values=SITE_LTY) +
  facet_wrap(~Class, ncol=1) +
  labs(x="Year", y="Log-abundance index (rel. 1997)",
       title="8-Site Population Trends: Breeding vs Haul-out",
       subtitle=paste0("Bands = ", CI_LABEL, "; dashed = DR & PB (haul-out only)"),
       caption="DR/PB pup states unobserved (NA) — excluded from pup panel") +
  theme_seal() + theme(legend.position="right")
ggsave("Output/Plots/8site_log_abundance.jpeg", p1, width=22, height=30, units="cm")

# ── PLOT 1b: smoothed trajectories (real scale) ─────────────────────────────
d_smoothed <- d_plot %>%
  mutate(Class=factor(Class, levels=c("Pup","Adult","Molt")),
         estimate=exp(log_est), lo89_n=exp(lo89), hi89_n=exp(hi89)) %>%
  dplyr::filter(!(Class=="Pup" & Site %in% c("DR","PB")))
p1b <- ggplot(d_smoothed, aes(Year, estimate, colour=Site, linetype=Site, group=Site)) +
  geom_vline(xintercept=breakpoint_year, linetype=2, colour="grey40", linewidth=0.6) +
  geom_ribbon(aes(ymin=lo89_n, ymax=hi89_n, fill=Site), alpha=0.12, colour=NA) +
  geom_line(linewidth=1) +
  scale_colour_manual(values=SITE_COLS) + scale_fill_manual(values=SITE_COLS) +
  scale_linetype_manual(values=SITE_LTY) +
  scale_y_continuous(labels=scales::comma) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="Estimated abundance",
       title="MARSS Smoothed State Trajectories — All 8 Sites",
       subtitle=paste0("Bands = ", CI_LABEL, "; dashed = DR (purple) & PB (cyan)"),
       caption=paste0("Vertical dashed = ", breakpoint_year,
                      " breakpoint. Pup panel excludes DR/PB.")) +
  theme_seal() + theme(legend.position="right")
ggsave("Output/Plots/8site_smoothed_states.jpeg", p1b, width=22, height=30, units="cm")

# ── PLOT 1c: total abundance summed across sites ────────────────────────────
d_total <- d_plot %>%
  mutate(Class=factor(Class, levels=c("Pup","Adult","Molt")),
         estimate=exp(log_est), lo89_n=exp(lo89), hi89_n=exp(hi89)) %>%
  group_by(Year, Class) %>%
  summarise(total=sum(estimate, na.rm=TRUE), total_lo=sum(lo89_n, na.rm=TRUE),
            total_hi=sum(hi89_n, na.rm=TRUE), .groups="drop")
p1c <- ggplot(d_total, aes(Year, total)) +
  geom_ribbon(aes(ymin=total_lo, ymax=total_hi), alpha=0.20, fill="#2166AC", colour=NA) +
  geom_line(colour="#2166AC", linewidth=1.2) +
  scale_y_continuous(limits=c(0,NA), labels=scales::comma, expand=c(0,0)) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="Total estimated abundance",
       title="MARSS Total Estimated Abundance — Sum Across All 8 Sites",
       subtitle=paste0("Bands = ", CI_LABEL),
       caption="Pup total = 6 breeding sites only (DR/PB pups unobserved)") +
  theme_seal()
ggsave("Output/Plots/8site_total_abundance.jpeg", p1c, width=22, height=24, units="cm")

# ── PLOT 2: haul-out sites vs breeding-site mean ────────────────────────────
breed_mean <- d_diff %>%
  dplyr::filter(Class %in% c("Adult","Molt"), !Site %in% c("DR","PB")) %>%
  group_by(Year, Class) %>%
  summarise(log_diff=mean(log_diff), lo89_d=mean(lo89_d), hi89_d=mean(hi89_d),
            .groups="drop") %>%
  mutate(Group="Breeding site mean")
haulout_sites <- d_diff %>%
  dplyr::filter(Site %in% c("DR","PB"), Class %in% c("Adult","Molt")) %>%
  mutate(Group=paste0(Site, " (haul-out)"))
compare_df <- bind_rows(
  breed_mean %>% select(Year, Class, log_diff, lo89_d, hi89_d, Group),
  haulout_sites %>% select(Year, Class, log_diff, lo89_d, hi89_d, Group))
p2 <- ggplot(compare_df, aes(Year, log_diff, colour=Group, fill=Group, linetype=Group)) +
  geom_ribbon(aes(ymin=lo89_d, ymax=hi89_d), alpha=0.15, colour=NA) +
  geom_line(linewidth=1.1) +
  geom_hline(yintercept=0, linetype=2, colour="grey40") +
  scale_colour_manual(values=c("Breeding site mean"="grey30",
                               "DR (haul-out)"=SITE_COLS["DR"], "PB (haul-out)"=SITE_COLS["PB"])) +
  scale_fill_manual(values=c("Breeding site mean"="grey50",
                             "DR (haul-out)"=SITE_COLS["DR"], "PB (haul-out)"=SITE_COLS["PB"])) +
  scale_linetype_manual(values=c("Breeding site mean"=1, "DR (haul-out)"=2, "PB (haul-out)"=2)) +
  facet_wrap(~Class, ncol=1) +
  labs(x="Year", y="Log-abundance index rel. 1997",
       title="Haul-out Sites vs Breeding Site Mean Trend",
       subtitle=paste0("Do DR & PB track the breeding population? Bands = ", CI_LABEL),
       colour=NULL, fill=NULL, linetype=NULL) +
  theme_seal()
ggsave("Output/Plots/8site_haulout_vs_breeding.jpeg", p2, width=22, height=24, units="cm")

# ── PLOT 3: between-site correlation (adult class) ──────────────────────────
adult_states <- d_plot %>%
  dplyr::filter(Class=="Adult") %>%
  select(Year, Site, log_est) %>%
  pivot_wider(names_from=Site, values_from=log_est) %>%
  select(-Year)
cor_mat_8site <- cor(adult_states, use="pairwise.complete.obs")
site_order <- c("BL","DE","DP","PRH","TB","TP","DR","PB")
cor_long_8site <- as.data.frame(cor_mat_8site) %>%
  tibble::rownames_to_column("Site1") %>%
  pivot_longer(-Site1, names_to="Site2", values_to="r") %>%
  mutate(Site1=factor(Site1, levels=site_order), Site2=factor(Site2, levels=site_order))
p3 <- ggplot(cor_long_8site, aes(Site1, fct_rev(Site2), fill=r)) +
  geom_tile(colour="white", linewidth=0.5) +
  geom_text(aes(label=sprintf("%.2f", r), colour=abs(r)>0.5), size=3.5) +
  scale_colour_manual(values=c("TRUE"="white","FALSE"="grey30"), guide="none") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B",
                       midpoint=0, limits=c(-1,1), name="r") +
  geom_hline(yintercept=2.5, linewidth=1.2, colour="grey40") +
  geom_vline(xintercept=6.5, linewidth=1.2, colour="grey40") +
  labs(x=NULL, y=NULL,
       title="Between-site Synchrony: Adult Correlation (All 8 Sites)",
       subtitle="Grey lines separate breeding (left/bottom) from haul-out (right/top)") +
  coord_fixed() + theme_seal() +
  theme(axis.text=element_text(size=11), panel.grid=element_blank(), legend.position="right")
ggsave("Output/Plots/8site_correlation_matrix.jpeg", p3, width=22, height=20, units="cm")

# ── PLOT 4: % of abundance at haul-out sites (adult+molt) ────────────────────
d_abund <- d_plot %>%
  mutate(est_count=exp(log_est)) %>%
  dplyr::filter(Class != "Pup") %>%
  group_by(Year) %>% mutate(total_all=sum(est_count, na.rm=TRUE)) %>%
  group_by(Year, Site) %>%
  summarise(site_total=sum(est_count, na.rm=TRUE), total_all=first(total_all),
            .groups="drop") %>%
  mutate(SiteType=SITE_TYPE[Site])
haulout_prop <- d_abund %>%
  dplyr::filter(SiteType=="Haul-out only") %>%
  group_by(Year) %>%
  summarise(haulout_total=sum(site_total), total_all=first(total_all), .groups="drop") %>%
  mutate(prop=haulout_total/total_all)
p4 <- ggplot(haulout_prop, aes(Year, prop*100)) +
  geom_area(fill="#984EA3", alpha=0.4) +
  geom_line(colour="#984EA3", linewidth=1.2) +
  scale_y_continuous(limits=c(0,NA), labels=function(x) paste0(x,"%")) +
  labs(x="Year", y="% of total adult + molt abundance",
       title="Proportion of Estimated Abundance at Haul-out-Only Sites (DR + PB)",
       subtitle="Adult and molting counts only (no pup counts at DR/PB)") +
  theme_seal()
ggsave("Output/Plots/8site_haulout_proportion.jpeg", p4, width=22, height=14, units="cm")

# ── PLOT 5: covariate effects forest plot ───────────────────────────────────
coef_df <- tidy(CIs_8site) %>%
  tibble::as_tibble() %>%                          # tidy returns data.frame
  dplyr::filter(str_detect(term, "^C\\.")) %>%
  mutate(
    Class = case_when(str_detect(term,"_A$|_A,")~"Adult",
                      str_detect(term,"_M$|_M,")~"Molting",
                      str_detect(term,"_P$|_P,")~"Pup", TRUE~"All"),
    Significant = case_when(conf.up<0~"Negative", conf.low>0~"Positive", TRUE~"Neutral"),
    Group = case_when(str_detect(term,"MOCI")~"MOCI", str_detect(term,"Dist")~"Disturbance",
                      str_detect(term,"Coy")~"Coyote", str_detect(term,"ES|eSeal")~"Elephant seal",
                      TRUE~"Other"),
    Site = str_extract(term, "BL|DE|DP|DR|PB|PRH|TB|TP"),
    SiteType = case_when(Site %in% c("DR","PB")~"Haul-out", !is.na(Site)~"Breeding", TRUE~"Shared"),
    Label = term %>% str_remove("^C\\.") %>% str_replace_all("_"," ")) %>%
  arrange(Group, Class, Label) %>%
  mutate(Label = fct_inorder(Label))
grp_cols <- c(MOCI="#2166AC", Disturbance="#8C510A", Coyote="#B2182B", "Elephant seal"="#762A83")
p5 <- ggplot(coef_df, aes(estimate, Label)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey40") +
  geom_linerange(aes(xmin=conf.low, xmax=conf.up, colour=Group, alpha=SiteType), linewidth=0.8) +
  geom_point(aes(colour=Group, shape=Significant, size=SiteType=="Haul-out")) +
  scale_colour_manual(values=grp_cols) +
  scale_shape_manual(values=c("Negative"=25,"Neutral"=16,"Positive"=24), name=NULL) +
  scale_alpha_manual(values=c("Breeding"=0.8,"Haul-out"=1.0,"Shared"=0.7), name="Site type") +
  scale_size_manual(values=c("TRUE"=3.5,"FALSE"=2.5), guide="none") +
  facet_grid(Group ~ Class, scales="free_y", space="free_y") +
  labs(x="Coefficient estimate", y=NULL,
       title="MARSS Covariate Effects — 8-Site Analysis",
       subtitle=paste0(CI_LABEL, "; triangles = significant; larger = haul-out (DR/PB)"),
       colour="Covariate group") +
  theme_seal(base_size=12) +
  theme(strip.text.y=element_text(angle=0, face="bold"), panel.grid.major.y=element_blank())
ggsave("Output/Plots/8site_covariate_forest.jpeg", p5, width=32, height=36, units="cm", dpi=200)

# ── PLOT 7: total abundance by site type ────────────────────────────────────
d_abund_class <- d_plot %>%
  mutate(est_count=exp(log_est), lo_count=exp(lo89), hi_count=exp(hi89)) %>%
  group_by(Year, SiteType, Class) %>%
  summarise(total=sum(est_count, na.rm=TRUE), total_lo=sum(lo_count, na.rm=TRUE),
            total_hi=sum(hi_count, na.rm=TRUE), .groups="drop")
p7 <- ggplot(d_abund_class %>% dplyr::filter(Class != "Pup" | SiteType=="Breeding"),
             aes(Year, total, colour=SiteType, fill=SiteType)) +
  geom_ribbon(aes(ymin=total_lo, ymax=total_hi), alpha=0.15, colour=NA) +
  geom_line(linewidth=1.1) +
  scale_y_continuous(limits=c(0,NA), labels=scales::comma, expand=c(0,0)) +
  scale_colour_manual(values=c("Breeding"="#2166AC","Haul-out only"="#984EA3")) +
  scale_fill_manual(values=c("Breeding"="#AECDE8","Haul-out only"="#D9B3E8")) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="Estimated abundance",
       title="Total Estimated Abundance by Site Type and Class",
       subtitle=paste0("Breeding vs Haul-out (DR/PB). ", CI_LABEL),
       colour=NULL, fill=NULL,
       caption="Pup panel = breeding sites only") +
  theme_seal()
ggsave("Output/Plots/8site_total_by_type.jpeg", p7, width=22, height=30, units="cm")

# ── PLOT 8: % change by site ────────────────────────────────────────────────
d_pct <- d_diff %>%
  mutate(pct_change=(exp(log_diff)-1)*100,
         pct_lo=(exp(lo89_d)-1)*100, pct_hi=(exp(hi89_d)-1)*100)
p8 <- ggplot(d_pct, aes(Year, pct_change, colour=Site, group=Site, linetype=Site)) +
  geom_hline(yintercept=0, linetype=2, colour="grey40") +
  geom_ribbon(aes(ymin=pct_lo, ymax=pct_hi, fill=Site), alpha=0.10, colour=NA) +
  geom_line(linewidth=0.9) +
  scale_colour_manual(values=SITE_COLS) + scale_fill_manual(values=SITE_COLS) +
  scale_linetype_manual(values=SITE_LTY) +
  facet_wrap(~Class, ncol=1, scales="free_y") +
  labs(x="Year", y="% change from 1997",
       title="Percentage Change in Estimated Abundance from 1997 — All 8 Sites",
       subtitle=paste0("Dashed = DR & PB (haul-out only). Bands = ", CI_LABEL)) +
  theme_seal() + theme(legend.position="right")
ggsave("Output/Plots/8site_pct_change.jpeg", p8, width=22, height=28, units="cm")

# ── Summary statistics ──────────────────────────────────────────────────────
cat("\n── Haul-out site statistics ────────────────────────────────────────\n")
dr_pb_corr <- cor_mat_8site[c("DR","PB"), c("BL","DE","DP","PRH","TB","TP")]
cat("  DR vs breeding (mean r):", round(mean(dr_pb_corr["DR",]), 3), "\n")
cat("  PB vs breeding (mean r):", round(mean(dr_pb_corr["PB",]), 3), "\n")
cat("  Mean haul-out % of total:", round(mean(haulout_prop$prop)*100, 1), "%\n")
cat("  Range:", round(min(haulout_prop$prop)*100, 1), "-",
    round(max(haulout_prop$prop)*100, 1), "%\n")

cat("\nAll 8-site plots -> Output/Plots/8site_*\n")
cat("Next: source(\"Code/12_data_prep_legacy.R\")  (optional long-term analysis)\n")