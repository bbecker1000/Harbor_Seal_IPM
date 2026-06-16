# ============================================================================
# HARBOR SEAL IPM v3.2 — PORTFOLIO ANALYSIS
# Run after model has completed. Requires results.real (or reload from RDS).
# ============================================================================

# ── Load fit if starting fresh ───────────────────────────────────────────────
# If results.real exists in session, skip this block.
# Otherwise load from saved RDS:
#
#   source("harbor_seal_ipm_v3.2.R")
#   source("Code/harbor_seal_ipm_v3.2_plots.R")
#   out      <- load_seal_results("IPM_v3.2_real")
#   fit      <- out$fit
#   sim_data <- out$sim_data
#
# Or if results.real is in session:
fit      <- Results.real$fit
sim_data <- Results.real$data

prefix     <- "IPM_v3.2_real"
site_names <- c("BL","DE","DP","PRH","TB","TP")

dir.create("Output",       showWarnings=FALSE)
dir.create("Output/Plots", showWarnings=FALSE)

# ── Extract lambda draws ─────────────────────────────────────────────────────
years  <- sim_data$years
S      <- length(site_names)
T      <- length(years)

ldraws <- fit$draws(variables="lambda", format="matrix")
lmat   <- matrix(NA, S, T-1, dimnames=list(site_names, years[1:(T-1)]))
for (s in 1:S) for (t in 1:(T-1)) {
  cn <- paste0("lambda[",s,",",t,"]")
  if (cn %in% colnames(ldraws)) lmat[s,t] <- mean(ldraws[,cn])
}

# ── Core portfolio statistics ────────────────────────────────────────────────
cv_meta  <- sd(colMeans(lmat, na.rm=TRUE)) / mean(colMeans(lmat, na.rm=TRUE))
cv_sites <- mean(apply(lmat, 1, function(x) sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE)))
per      <- cv_meta / cv_sites
lcor     <- cor(t(lmat), use="pairwise.complete.obs")

cat(sprintf("\nPortfolio Effect Ratio: %.3f  (< 1 = buffering)\n", per))
cat(sprintf("Mean site correlation:  %.3f\n", mean(lcor[lower.tri(lcor)])))

# ── Extract phi_pup ──────────────────────────────────────────────────────────
phi_d <- fit$draws(variables="phi_pup", format="matrix")
phi_m <- matrix(NA, S, T, dimnames=list(site_names, years))
for (s in 1:S) for (t in 1:T) {
  cn <- paste0("phi_pup[",s,",",t,"]")
  if (cn %in% colnames(phi_d)) phi_m[s,t] <- mean(phi_d[,cn])
}

# ── Best/worst sites each year ───────────────────────────────────────────────
bw <- tibble(
  Year       = years[1:(T-1)],
  Best_Site  = site_names[apply(lmat, 2, which.max)],
  Worst_Site = site_names[apply(lmat, 2, which.min)]
)

# =============================================================================
# TABLES
# =============================================================================

# Table 1: Portfolio summary statistics
port_summary_tbl <- tibble(
  Metric = c(
    "Portfolio Effect Ratio (CV_meta / CV_sites)",
    "Mean between-site correlation (lambda)",
    "Min between-site correlation (lambda)",
    "Max between-site correlation (lambda)",
    "Overall mean lambda (across sites and years)",
    "SD of annual mean lambda"
  ),
  Value = c(
    round(per, 3),
    round(mean(lcor[lower.tri(lcor)]), 3),
    round(min(lcor[lower.tri(lcor)]), 3),
    round(max(lcor[lower.tri(lcor)]), 3),
    round(mean(lmat, na.rm=TRUE), 3),
    round(sd(colMeans(lmat, na.rm=TRUE)), 3)
  ),
  Interpretation = c(
    if (per < 1) "Portfolio buffering present (< 1)" else "No buffering (>= 1)",
    if (mean(lcor[lower.tri(lcor)]) < 0.3) "Low synchrony — strong portfolio effect"
    else if (mean(lcor[lower.tri(lcor)]) < 0.6) "Moderate synchrony"
    else "High synchrony — weak portfolio effect",
    "—", "—",
    if (mean(lmat,na.rm=TRUE) > 1) "Population growing on average"
    else "Population declining on average",
    "—"
  )
)
cat("\n── Portfolio Summary Table ──────────────────────────────────────\n")
print(port_summary_tbl, n=nrow(port_summary_tbl))
write_csv(port_summary_tbl, paste0("Output/",prefix,"_portfolio_summary.csv"))

# Table 2: Per-site lambda statistics
lambda_site_tbl <- as.data.frame(lmat) |>
  tibble::rownames_to_column("Site") |>
  mutate(
    Mean_lambda = round(rowMeans(across(where(is.numeric)), na.rm=TRUE), 3),
    SD_lambda   = round(apply(across(where(is.numeric)), 1,
                              function(x) sd(x, na.rm=TRUE)), 3),
    CV_lambda   = round(SD_lambda / Mean_lambda, 3),
    Times_best  = as.integer(table(bw$Best_Site)[Site])  |> tidyr::replace_na(0L),
    Times_worst = as.integer(table(bw$Worst_Site)[Site]) |> tidyr::replace_na(0L)
  ) |>
  select(Site, Mean_lambda, SD_lambda, CV_lambda, Times_best, Times_worst)
cat("\n── Per-Site Lambda Table ────────────────────────────────────────\n")
print(lambda_site_tbl)
write_csv(lambda_site_tbl, paste0("Output/",prefix,"_portfolio_lambda_by_site.csv"))

# Table 3: Pairwise between-site correlations (long format)
cor_long <- as.data.frame(lcor) |>
  tibble::rownames_to_column("Site1") |>
  pivot_longer(-Site1, names_to="Site2", values_to="r") |>
  filter(Site1 < Site2) |>
  mutate(
    r = round(r, 3),
    Interpretation = case_when(
      abs(r) < 0.3 ~ "Low",
      abs(r) < 0.6 ~ "Moderate",
      TRUE          ~ "High"
    )
  ) |>
  arrange(desc(abs(r)))
cat("\n── Between-Site Correlation Table ───────────────────────────────\n")
print(cor_long)
write_csv(cor_long, paste0("Output/",prefix,"_portfolio_correlations.csv"))

# =============================================================================
# PLOTS
# =============================================================================

ldf <- expand.grid(Site=site_names, Year=years[1:(T-1)]) |>
  mutate(Site=factor(Site,levels=site_names), lambda=as.vector(t(lmat)))

# 1. Lambda heatmap
p_heat <- ggplot(ldf, aes(x=Year, y=Site, fill=lambda)) +
  geom_tile(color="white", linewidth=0.5) +
  scale_fill_gradient2(low="red3", mid="white", high="darkgreen",
                       midpoint=1, limits=c(0.7,1.3), oob=scales::squish, name="\u03bb") +
  geom_text(aes(label=sprintf("%.2f",lambda)), size=2.5) +
  labs(x="Year", y="Site", title="Site-Specific \u03bb by Year") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8), panel.grid=element_blank())
ggsave(paste0("Output/Plots/",prefix,"_portfolio_lambda_heatmap.jpeg"),
       p_heat, width=35, height=15, units="cm")

# 2. Pup survival heatmap
phidf <- expand.grid(Site=site_names, Year=years) |>
  mutate(Site=factor(Site,levels=site_names), phi=as.vector(t(phi_m)))
p_phi <- ggplot(phidf, aes(x=Year, y=Site, fill=phi)) +
  geom_tile(color="white", linewidth=0.5) +
  scale_fill_viridis_c(name="\u03c6_pup", option="plasma",
                       limits=c(0.10,0.50), oob=scales::squish) +
  geom_text(aes(label=sprintf("%.2f",phi)), size=2.2) +
  labs(x="Year", y="Site", title="Site-Specific Pup Survival by Year") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8), panel.grid=element_blank())
ggsave(paste0("Output/Plots/",prefix,"_portfolio_pup_survival_heatmap.jpeg"),
       p_phi, width=35, height=15, units="cm")

# 3. Asynchrony line plot
p_async <- ggplot(ldf, aes(x=Year, y=lambda, color=Site, group=Site)) +
  geom_hline(yintercept=1, linetype=2, color="gray50") +
  geom_line(linewidth=1, alpha=0.7) + geom_point(size=2) +
  scale_color_brewer(palette="Dark2") +
  labs(x="Year", y=expression(lambda),
       title="Site-Level Population Growth \u2014 Asynchrony") +
  theme_seal() + theme(legend.position="bottom")
ggsave(paste0("Output/Plots/",prefix,"_portfolio_asynchrony.jpeg"),
       p_async, width=30, height=18, units="cm")

# 4. Best/worst sites tile
p_bw <- bw |>
  pivot_longer(c(Best_Site,Worst_Site), names_to="Performance", values_to="Site") |>
  mutate(Performance=recode(Performance, Best_Site="Best", Worst_Site="Worst"),
         Site=factor(Site, levels=site_names)) |>
  ggplot(aes(x=Year, y=Site, fill=Performance)) +
  geom_tile(color="white", linewidth=1, alpha=0.8) +
  scale_fill_manual(values=c(Best="darkgreen", Worst="red3")) +
  labs(x="Year", y="Site", title="Best and Worst Performing Sites Each Year") +
  theme_seal() +
  theme(axis.text.x=element_text(angle=45,hjust=1,size=8), panel.grid=element_blank())
ggsave(paste0("Output/Plots/",prefix,"_portfolio_best_worst.jpeg"),
       p_bw, width=35, height=12, units="cm")

# 5. Correlation matrix tile
cordf <- expand.grid(Site1=site_names, Site2=site_names) |>
  mutate(r=as.vector(lcor))
p_cor <- ggplot(cordf, aes(x=Site1, y=Site2, fill=r)) +
  geom_tile(color="white") +
  geom_text(aes(label=sprintf("%.2f",r)), size=4) +
  scale_fill_gradient2(low="blue", mid="white", high="red",
                       midpoint=0, limits=c(-1,1)) +
  labs(x="", y="", title="Between-Site Correlation in \u03bb") +
  theme_seal() + theme(panel.grid=element_blank()) + coord_fixed()
ggsave(paste0("Output/Plots/",prefix,"_portfolio_correlation.jpeg"),
       p_cor, width=20, height=18, units="cm")

# 6. CV by site bar chart
cv_site <- lambda_site_tbl |>
  ggplot(aes(x=reorder(Site, CV_lambda), y=CV_lambda, fill=Mean_lambda)) +
  geom_col(alpha=0.85) +
  geom_text(aes(label=sprintf("%.3f", CV_lambda)), vjust=-0.3, size=4) +
  scale_fill_gradient2(low="red3", mid="white", high="darkgreen",
                       midpoint=1, name=expression(bar(lambda))) +
  labs(x="Site", y="CV of \u03bb",
       title="Among-Year Variability in Population Growth by Site",
       subtitle="Lower CV = more stable site") +
  theme_seal() + theme(legend.position="right")
ggsave(paste0("Output/Plots/",prefix,"_portfolio_cv_by_site.jpeg"),
       cv_site, width=22, height=14, units="cm")

# 7. Pairwise lambda scatter grid
ldf_wide   <- as.data.frame(t(lmat)) |> setNames(site_names)
site_pairs <- combn(site_names, 2, simplify=FALSE)
pair_plots <- lapply(site_pairs, function(pair) {
  d <- data.frame(x_val=ldf_wide[[pair[1]]], y_val=ldf_wide[[pair[2]]])
  r <- cor(d$x_val, d$y_val, use="complete.obs")
  ggplot(d, aes(x=.data[["x_val"]], y=.data[["y_val"]])) +
    geom_point(alpha=0.6, colour=SEAL_COLS$pop, size=2) +
    geom_smooth(method="lm", se=FALSE, colour="gray40", linewidth=0.7) +
    geom_abline(slope=1, intercept=0, linetype=2, colour="gray70") +
    labs(x=paste0(pair[1]," \u03bb"), y=paste0(pair[2]," \u03bb"),
         title=sprintf("%s vs %s  (r=%.2f)", pair[1], pair[2], r)) +
    theme_seal(base_size=12)
})
p_pairs <- wrap_plots(pair_plots, ncol=3) +
  plot_annotation(title="Pairwise Site Synchrony in Population Growth (\u03bb)")
ggsave(paste0("Output/Plots/",prefix,"_portfolio_pairwise_lambda.jpeg"),
       p_pairs, width=36, height=30, units="cm")

# =============================================================================
# COLLECT RESULTS
# =============================================================================

portfolio.real <- list(
  # Plots
  lambda_heatmap     = p_heat,
  phi_heatmap        = p_phi,
  asynchrony         = p_async,
  best_worst         = p_bw,
  correlation        = p_cor,
  cv_by_site         = cv_site,
  pairwise           = p_pairs,
  # Tables
  summary_table      = port_summary_tbl,
  lambda_site_table  = lambda_site_tbl,
  correlation_table  = cor_long,
  # Raw data
  lambda_matrix      = lmat,
  phi_matrix         = phi_m,
  # Scalars
  portfolio_effect_ratio = per,
  mean_correlation       = mean(lcor[lower.tri(lcor)])
)

cat(sprintf("\n── Portfolio analysis complete ──────────────────────────────\n"))
cat(sprintf("   7 plots  -> Output/Plots/%s_portfolio_*.jpeg\n", prefix))
cat(sprintf("   3 tables -> Output/%s_portfolio_*.csv\n", prefix))
cat(sprintf("   Portfolio Effect Ratio: %.3f\n", per))
cat(sprintf("   Mean site correlation:  %.3f\n", mean(lcor[lower.tri(lcor)])))

