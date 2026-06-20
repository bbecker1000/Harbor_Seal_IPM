# 06a_amj_molt_plot



draws  <- out$fit$draws(format = "df")
beta_v  <- draws$beta_moci_amj_molt
logit_v <- draws$detect_molt_logit
xr     <- seq(-2, 2, length.out = 100)
idx    <- sample(seq_along(beta_v), 500)

# Build summary directly — no aggregate()
sm <- do.call(rbind, lapply(xr, function(x) {
  vals <- plogis(logit_v[idx] + beta_v[idx] * x)
  data.frame(
    cov_val = x,
    mn = mean(vals),
    lo = quantile(vals, 0.055),
    hi = quantile(vals, 0.945)
  )
}))

p4e <- ggplot(sm, aes(x = cov_val)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "blue") +
  geom_line(aes(y = mn), linewidth = 1.2, colour = "blue3") +
  geom_hline(yintercept = mean(plogis(logit_v)),
             linetype = 2, colour = "gray50") +
  geom_vline(xintercept = 0, linetype = 2, colour = "gray50") +
  labs(x = "MOCI Spring AMJ (SD)",
       y = "Molt detection probability",
       title = "AMJ \u2192 Molt Detection") +
  theme_seal()

ggsave("Output/Plots/IPM_v3.2_real_effects_amj_molt.jpeg",
       p4e, width = 20, height = 12, units = "cm")