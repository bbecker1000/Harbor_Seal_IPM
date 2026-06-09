# =============================================================================
# Harbor Seal IPM v3.2 — Stage-Structured Life-Cycle Diagram
# Circular nodes version
# =============================================================================

library(tidyverse)
library(ggforce)   # geom_circle()

# ── Colour palette ────────────────────────────────────────────────────────────
fill_pup  <- "white"
fill_juv  <- "white"
fill_adF  <- "white"
fill_adM  <- "white"

# ── Node centres ─────────────────────────────────────────────────────────────
px  <- 1.00;  py  <-  0.00   # Pup
jfx <- 3.40;  jfy <-  1.00   # Juv Female
jmx <- 3.40;  jmy <- -1.00   # Juv Male
afx <- 5.90;  afy <-  1.00   # Adult Female
amx <- 5.90;  amy <- -1.00   # Adult Male

r        <- 0.52   # circle radius
gap      <- 0.09   # clearance: arrow tip to circle edge
loop_gap <- 0.06   # clearance: stasis loop to circle top/bottom

# ── Arrow styles ──────────────────────────────────────────────────────────────
ar       <- arrow(length = unit(8, "pt"), type = "closed")
ar_loop  <- arrow(length = unit(6, "pt"), type = "closed")
ar_repro <- arrow(length = unit(6, "pt"), type = "closed")

lw      <- 1
lw_loop <- 0.85
lw_rep  <- 0.85

# ── Helper: point on circle edge in direction of target, offset by gap ────────
edge <- function(x0, y0, x1, y1, radius = r, g = gap) {
  theta <- atan2(y1 - y0, x1 - x0)
  c(x0 + (radius + g) * cos(theta),
    y0 + (radius + g) * sin(theta))
}

# Pre-compute arrow endpoints
pup_to_jf_d <- edge(px,  py,  jfx, jfy)
pup_to_jf_a <- edge(jfx, jfy, px,  py)
pup_to_jm_d <- edge(px,  py,  jmx, jmy)
pup_to_jm_a <- edge(jmx, jmy, px,  py)
jf_to_af_d  <- edge(jfx, jfy, afx, afy)
jf_to_af_a  <- edge(afx, afy, jfx, jfy)
jm_to_am_d  <- edge(jmx, jmy, amx, amy)
jm_to_am_a  <- edge(amx, amy, jmx, jmy)

# Fecundity arc: depart bottom of AdultF, arrive right of Pup
fec_d <- c(afx - 0.20, afy - r - gap)
fec_a <- c(px  + r + gap, py + 0.10)

life_cycle <- ggplot() +
  
  # ── Circles ───────────────────────────────────────────────────────────────────
  geom_circle(aes(x0=px,  y0=py,  r=r), fill=fill_pup, color=NA,      linewidth=0) +
  geom_circle(aes(x0=jfx, y0=jfy, r=r), fill=fill_juv, color=NA,      linewidth=0) +
  geom_circle(aes(x0=jmx, y0=jmy, r=r), fill=fill_juv, color=NA,      linewidth=0) +
  geom_circle(aes(x0=afx, y0=afy, r=r), fill=fill_adF, color=NA,      linewidth=0) +
  geom_circle(aes(x0=amx, y0=amy, r=r), fill=fill_adM, color=NA,      linewidth=0) +
  geom_circle(aes(x0=px,  y0=py,  r=r), fill=NA, color="black",       linewidth=1.4) +
  geom_circle(aes(x0=jfx, y0=jfy, r=r), fill=NA, color="black",       linewidth=1.4) +
  geom_circle(aes(x0=jmx, y0=jmy, r=r), fill=NA, color="black",       linewidth=1.4) +
  geom_circle(aes(x0=afx, y0=afy, r=r), fill=NA, color="black",       linewidth=1.4) +
  geom_circle(aes(x0=amx, y0=amy, r=r), fill=NA, color="black",       linewidth=1.4) +
  
  # ── Node labels ──────────────────────────────────────────────────────────────
  annotate("text", x=px,  y=py+0.13,  label="bold('N'['pup'])",       parse=TRUE, size=6) +
  annotate("text", x=px,  y=py-0.15,  label="italic('sex-neutral')",   parse=TRUE, size=4) +
  annotate("text", x=jfx, y=jfy+0.13, label="bold('N'['juv,F'])",     parse=TRUE, size=6) +
  annotate("text", x=jfx, y=jfy-0.15, label="italic('3-yr  \u2640')",  parse=TRUE, size=4) +
  annotate("text", x=jmx, y=jmy+0.13, label="bold('N'['juv,M'])",     parse=TRUE, size=6) +
  annotate("text", x=jmx, y=jmy-0.15, label="italic('3-yr  \u2642')",  parse=TRUE, size=4) +
  annotate("text", x=afx, y=afy+0.13, label="bold('N'['adult,F'])",   parse=TRUE, size=6) +
  annotate("text", x=afx, y=afy-0.15, label="'\u2640'",                parse=TRUE, size=5) +
  annotate("text", x=amx, y=amy+0.13, label="bold('N'['adult,M'])",   parse=TRUE, size=6) +
  annotate("text", x=amx, y=amy-0.15, label="'\u2642'",                parse=TRUE, size=5) +
  
  # ── Stage-transition arrows ───────────────────────────────────────────────────
  annotate("segment",
           x=pup_to_jf_d[1], y=pup_to_jf_d[2],
           xend=pup_to_jf_a[1], yend=pup_to_jf_a[2],
           arrow=ar, color="black", linewidth=lw) +
  annotate("segment",
           x=pup_to_jm_d[1], y=pup_to_jm_d[2],
           xend=pup_to_jm_a[1], yend=pup_to_jm_a[2],
           arrow=ar, color="black", linewidth=lw) +
  annotate("segment",
           x=jf_to_af_d[1], y=jf_to_af_d[2],
           xend=jf_to_af_a[1], yend=jf_to_af_a[2],
           arrow=ar, color="black", linewidth=lw) +
  annotate("segment",
           x=jm_to_am_d[1], y=jm_to_am_d[2],
           xend=jm_to_am_a[1], yend=jm_to_am_a[2],
           arrow=ar, color="black", linewidth=lw) +
  
  # ── Stasis loops ──────────────────────────────────────────────────────────────
  annotate("curve",
           x=jfx+r*0.55, y=jfy+r+loop_gap,
           xend=jfx-r*0.55, yend=jfy+r+loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop) +
  annotate("curve",
           x=jmx-r*0.55, y=jmy-r-loop_gap,
           xend=jmx+r*0.55, yend=jmy-r-loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop) +
  annotate("curve",
           x=afx+r*0.55, y=afy+r+loop_gap,
           xend=afx-r*0.55, yend=afy+r+loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop) +
  annotate("curve",
           x=amx-r*0.55, y=amy-r-loop_gap,
           xend=amx+r*0.55, yend=amy-r-loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop) +
  
  # ── Fecundity arrow ───────────────────────────────────────────────────────────
  annotate("curve",
           x=fec_d[1]-0.2, y=fec_d[2]+0.2, xend=fec_a[1], yend=fec_a[2]-0.1,
           curvature=-0.1, color="grey40", arrow=ar_repro,
           linewidth=lw_rep, linetype="solid") +
  
  # ── Arrow labels ─────────────────────────────────────────────────────────────
  annotate("text",
           x=(pup_to_jf_d[1]+pup_to_jf_a[1])/2 - 0.25,
           y=(pup_to_jf_d[2]+pup_to_jf_a[2])/2 + 0.25,
           label="phi[pup]~rho[F]", parse=TRUE, size=5, fontface="bold") +
  annotate("text",
           x=(pup_to_jm_d[1]+pup_to_jm_a[1])/2 - 0.25,
           y=(pup_to_jm_d[2]+pup_to_jm_a[2])/2 - 0.25,
           label="phi[pup]~(1-rho[F])", parse=TRUE, size=5, fontface="bold") +
  annotate("text",
           x=(jf_to_af_d[1]+jf_to_af_a[1])/2, y=jfy+0.25,
           label="frac(1,3)~phi[juv]", parse=TRUE, size=5, fontface="bold") +
  annotate("text",
           x=(jm_to_am_d[1]+jm_to_am_a[1])/2, y=jmy-0.25,
           label="frac(1,3)~phi[juv]", parse=TRUE, size=5, fontface="bold") +
  
  # Stasis labels
  annotate("text", x=jfx, y=jfy+r+0.23,
           label="frac(2,3)~phi[juv]", parse=TRUE, size=4) +
  annotate("text", x=jmx, y=jmy-r-0.23,
           label="frac(2,3)~phi[juv]", parse=TRUE, size=4.0) +
  annotate("text", x=afx, y=afy+r+0.23,
           label="phi[adult~F]", parse=TRUE, size=4.0) +
  annotate("text", x=amx, y=amy-r-0.23,
           label="phi[adult~M]", parse=TRUE, size=4.0) +
  
  # Fecundity label
  annotate("text",
           x=(px+afx)/2 + 0.3, y=-0.1,
           label="italic(f[t])~'(fecundity)'", parse=TRUE,
           size=3.2, fontface="italic") +
  
  # ── Theme ─────────────────────────────────────────────────────────────────────
  theme_void(base_size=13) +
  theme(
    plot.background = element_rect(fill="white", color=NA),
    plot.margin     = margin(20, 20, 25, 20),
    plot.title      = element_text(face="bold", size=14, hjust=0.5, margin=margin(b=6)),
    plot.subtitle   = element_text(size=10, hjust=0.5, color="grey40", margin=margin(b=10))
  ) +
  labs(
    title    = "Harbor Seal IPM v3.2 — Stage-Structured Life Cycle",
    subtitle = expression(
      phi[pup]*': pup survival  |  '*phi[juv]*': juvenile survival  |  '*
        phi[adult~F]*': adult female  |  '*phi[adult~M]*' = '*phi[adult~F]*' - '*delta[adult])
  ) +
  coord_equal(xlim=c(0.0, 7.1), ylim=c(-2.2, 2.6)) +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0))

life_cycle

ggsave("Output/Plots/life_cycle_v3.2.jpeg",
       life_cycle, width=32, height=22, units="cm", dpi=300)

