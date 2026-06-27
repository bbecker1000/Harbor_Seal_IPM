# =============================================================================
# Harbor Seal IPM — Stage-Structured Life-Cycle Diagrams
# Two versions:
#   (A) Marin IPM v3.3   — MOCI seasons as small colored text labels
#   (B) Regional IPM     — adds county-type modifier notation on fecundity
#
# MOCI label placement:
#   φ_pup arrows        → OND · AMJ · JFM (stacked, both arrows)
#   2/3 φ_juv stasis    → JFM (both sexes)
#   1/3 φ_juv → adult  → JFM (both sexes) — same φ_juv, same MOCI effect
#   φ_adult stasis      → JFM (both sexes)
#   fecundity arc       → OND
#   Regional only: fecundity arc also notes county modifier
# =============================================================================

library(tidyverse)
library(ggforce)

# ── Node centres ──────────────────────────────────────────────────────────────
px  <- 1.00;  py  <-  0.00
jfx <- 3.40;  jfy <-  1.00
jmx <- 3.40;  jmy <- -1.00
afx <- 5.90;  afy <-  1.00
amx <- 5.90;  amy <- -1.00

r        <- 0.52
gap      <- 0.09
loop_gap <- 0.06

# ── Arrow styles ──────────────────────────────────────────────────────────────
ar       <- arrow(length = unit(8, "pt"), type = "closed")
ar_loop  <- arrow(length = unit(6, "pt"), type = "closed")
ar_repro <- arrow(length = unit(6, "pt"), type = "closed")

lw      <- 1.00
lw_loop <- 0.85
lw_rep  <- 0.85

# ── MOCI label colours ────────────────────────────────────────────────────────
col_ond <- "#B03A3E"   # muted red    — OND
col_amj <- "#2E7D32"   # muted green  — AMJ
col_jfm <- "#1A5EA8"   # muted blue   — JFM
moci_sz <- 3.4         # font size for MOCI labels

# ── Helper: point on circle edge toward target ────────────────────────────────
edge <- function(x0, y0, x1, y1, radius = r, g = gap) {
  theta <- atan2(y1 - y0, x1 - x0)
  c(x0 + (radius + g) * cos(theta),
    y0 + (radius + g) * sin(theta))
}

# ── Pre-compute arrow endpoints ───────────────────────────────────────────────
pup_to_jf_d <- edge(px,  py,  jfx, jfy)
pup_to_jf_a <- edge(jfx, jfy, px,  py )
pup_to_jm_d <- edge(px,  py,  jmx, jmy)
pup_to_jm_a <- edge(jmx, jmy, px,  py )
jf_to_af_d  <- edge(jfx, jfy, afx, afy)
jf_to_af_a  <- edge(afx, afy, jfx, jfy)
jm_to_am_d  <- edge(jmx, jmy, amx, amy)
jm_to_am_a  <- edge(amx, amy, jmx, jmy)

fec_d <- c(afx - 0.20, afy - r - gap)
fec_a <- c(px  + r + gap, py + 0.10)

# Midpoints of pup survival arrows
mid_pup_jf <- c((pup_to_jf_d[1] + pup_to_jf_a[1]) / 2,
                (pup_to_jf_d[2] + pup_to_jf_a[2]) / 2)
mid_pup_jm <- c((pup_to_jm_d[1] + pup_to_jm_a[1]) / 2,
                (pup_to_jm_d[2] + pup_to_jm_a[2]) / 2)

# Midpoints of juv → adult transition arrows
mid_jf_af <- c((jf_to_af_d[1] + jf_to_af_a[1]) / 2,
               (jf_to_af_d[2] + jf_to_af_a[2]) / 2)
mid_jm_am <- c((jm_to_am_d[1] + jm_to_am_a[1]) / 2,
               (jm_to_am_d[2] + jm_to_am_a[2]) / 2)

# ── Base life-cycle layers ────────────────────────────────────────────────────
base_layers <- list(
  
  # Circles — fill first, then border
  geom_circle(aes(x0=px,  y0=py,  r=r), fill="white", color=NA, linewidth=0),
  geom_circle(aes(x0=jfx, y0=jfy, r=r), fill="white", color=NA, linewidth=0),
  geom_circle(aes(x0=jmx, y0=jmy, r=r), fill="white", color=NA, linewidth=0),
  geom_circle(aes(x0=afx, y0=afy, r=r), fill="white", color=NA, linewidth=0),
  geom_circle(aes(x0=amx, y0=amy, r=r), fill="white", color=NA, linewidth=0),
  geom_circle(aes(x0=px,  y0=py,  r=r), fill=NA, color="black", linewidth=1.4),
  geom_circle(aes(x0=jfx, y0=jfy, r=r), fill=NA, color="black", linewidth=1.4),
  geom_circle(aes(x0=jmx, y0=jmy, r=r), fill=NA, color="black", linewidth=1.4),
  geom_circle(aes(x0=afx, y0=afy, r=r), fill=NA, color="black", linewidth=1.4),
  geom_circle(aes(x0=amx, y0=amy, r=r), fill=NA, color="black", linewidth=1.4),
  
  # Node labels
  annotate("text", x=px,  y=py +0.13, label="bold('N'['pup'])",      parse=TRUE, size=6),
  annotate("text", x=px,  y=py -0.15, label="italic('sex-neutral')",  parse=TRUE, size=4),
  annotate("text", x=jfx, y=jfy+0.13, label="bold('N'['juv,F'])",    parse=TRUE, size=6),
  annotate("text", x=jfx, y=jfy-0.15, label="italic('3-yr  \u2640')", parse=TRUE, size=4),
  annotate("text", x=jmx, y=jmy+0.13, label="bold('N'['juv,M'])",    parse=TRUE, size=6),
  annotate("text", x=jmx, y=jmy-0.15, label="italic('3-yr  \u2642')", parse=TRUE, size=4),
  annotate("text", x=afx, y=afy+0.13, label="bold('N'['adult,F'])",  parse=TRUE, size=6),
  annotate("text", x=afx, y=afy-0.15, label="'\u2640'",               parse=TRUE, size=5),
  annotate("text", x=amx, y=amy+0.13, label="bold('N'['adult,M'])",  parse=TRUE, size=6),
  annotate("text", x=amx, y=amy-0.15, label="'\u2642'",               parse=TRUE, size=5),
  
  # Stage-transition arrows
  annotate("segment",
           x=pup_to_jf_d[1], y=pup_to_jf_d[2],
           xend=pup_to_jf_a[1], yend=pup_to_jf_a[2],
           arrow=ar, color="black", linewidth=lw),
  annotate("segment",
           x=pup_to_jm_d[1], y=pup_to_jm_d[2],
           xend=pup_to_jm_a[1], yend=pup_to_jm_a[2],
           arrow=ar, color="black", linewidth=lw),
  annotate("segment",
           x=jf_to_af_d[1], y=jf_to_af_d[2],
           xend=jf_to_af_a[1], yend=jf_to_af_a[2],
           arrow=ar, color="black", linewidth=lw),
  annotate("segment",
           x=jm_to_am_d[1], y=jm_to_am_d[2],
           xend=jm_to_am_a[1], yend=jm_to_am_a[2],
           arrow=ar, color="black", linewidth=lw),
  
  # Stasis loops
  annotate("curve",
           x=jfx+r*0.55, y=jfy+r+loop_gap,
           xend=jfx-r*0.55, yend=jfy+r+loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop),
  annotate("curve",
           x=jmx-r*0.55, y=jmy-r-loop_gap,
           xend=jmx+r*0.55, yend=jmy-r-loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop),
  annotate("curve",
           x=afx+r*0.55, y=afy+r+loop_gap,
           xend=afx-r*0.55, yend=afy+r+loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop),
  annotate("curve",
           x=amx-r*0.55, y=amy-r-loop_gap,
           xend=amx+r*0.55, yend=amy-r-loop_gap,
           curvature=1.6, ncp=25, color="black",
           arrow=ar_loop, linewidth=lw_loop),
  
  # Fecundity arrow (Adult F → Pup)
  annotate("curve",
           x=fec_d[1]-0.2, y=fec_d[2]+0.2,
           xend=fec_a[1], yend=fec_a[2]-0.1,
           curvature=-0.1, color="grey40",
           arrow=ar_repro, linewidth=lw_rep),
  
  # ── Vital-rate labels on transitions ─────────────────────────────────────────
  
  # Pup survival arrows
  annotate("text",
           x = mid_pup_jf[1] - 0.25,
           y = mid_pup_jf[2] + 0.26,
           label = "phi[pup]~rho[F]", parse = TRUE,
           size = 5, fontface = "bold"),
  annotate("text",
           x = mid_pup_jm[1] - 0.25,
           y = mid_pup_jm[2] - 0.26,
           label = "phi[pup]~(1-rho[F])", parse = TRUE,
           size = 5, fontface = "bold"),
  
  # Juv → Adult transition arrows
  annotate("text",
           x = mid_jf_af[1], y = jfy + 0.25,
           label = "frac(1,3)~phi[juv]", parse = TRUE,
           size = 5, fontface = "bold"),
  annotate("text",
           x = mid_jm_am[1], y = jmy - 0.25,
           label = "frac(1,3)~phi[juv]", parse = TRUE,
           size = 5, fontface = "bold"),
  
  # Stasis labels
  annotate("text", x = jfx, y = jfy + r + 0.23,
           label = "frac(2,3)~phi[juv]", parse = TRUE, size = 4),
  annotate("text", x = jmx, y = jmy - r - 0.23,
           label = "frac(2,3)~phi[juv]", parse = TRUE, size = 4),
  annotate("text", x = afx, y = afy + r + 0.23,
           label = "phi[adult~F]", parse = TRUE, size = 4),
  annotate("text", x = amx, y = amy - r - 0.23,
           label = "phi[adult~M]", parse = TRUE, size = 4),
  
  # Fecundity label
  annotate("text",
           x = (px + afx)/2 + 0.30, y = -0.10,
           label = "italic(f[t])~'(fecundity)'", parse = TRUE,
           size = 3.2, fontface = "italic")
)

# ── MOCI text label layers ────────────────────────────────────────────────────
# Small italic coloured season abbreviations placed near each affected label.
# The 1/3 φ_juv advancement arrows carry the same JFM label as the stasis
# loops because both depend on the common annual juvenile survival φ_juv.

moci_labels <- list(
  
  # ── Pup survival: OND · AMJ · JFM (both arrows) ──────────────────────────
  
  # Upper arrow (φ_pup ρ_F) — stack below the vital-rate label
  annotate("text",
           x = mid_pup_jf[1] - 0.25, y = mid_pup_jf[2] + 0.04,
           label = "italic('OND')", parse = TRUE,
           size = moci_sz, color = col_ond),
  annotate("text",
           x = mid_pup_jf[1] - 0.25, y = mid_pup_jf[2] - 0.12,
           label = "italic('AMJ')", parse = TRUE,
           size = moci_sz, color = col_amj),
  annotate("text",
           x = mid_pup_jf[1] - 0.25, y = mid_pup_jf[2] - 0.28,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  
  # Lower arrow (φ_pup (1−ρ_F)) — stack above the vital-rate label
  annotate("text",
           x = mid_pup_jm[1] - 0.25, y = mid_pup_jm[2] - 0.04,
           label = "italic('OND')", parse = TRUE,
           size = moci_sz, color = col_ond),
  annotate("text",
           x = mid_pup_jm[1] - 0.25, y = mid_pup_jm[2] + 0.12,
           label = "italic('AMJ')", parse = TRUE,
           size = moci_sz, color = col_amj),
  annotate("text",
           x = mid_pup_jm[1] - 0.25, y = mid_pup_jm[2] + 0.28,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  
  # ── Juvenile stasis loops: JFM ────────────────────────────────────────────
  # Female — below "2/3 φ_juv" (at y = jfy + r + 0.23 ≈ 1.75)
  annotate("text",
           x = jfx, y = jfy + r + 0.48,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  
  # Male — above "2/3 φ_juv" (at y = jmy - r - 0.23 ≈ −1.75)
  annotate("text",
           x = jmx, y = jmy - r - 0.48,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  
  # ── Juv → Adult advancement arrows: JFM ──────────────────────────────────
  # Same φ_juv, same MOCI effect as stasis.
  # Label placed just below "1/3 φ_juv" (female) and above "1/3 φ_juv" (male).
  annotate("text",
           x = mid_jf_af[1], y = jfy + 0.02,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  annotate("text",
           x = mid_jm_am[1], y = jmy - 0.02,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  
  # ── Adult stasis loops: JFM ───────────────────────────────────────────────
  # Female — below "φ_adult F" (at y = afy + r + 0.23 ≈ 1.75)
  annotate("text",
           x = afx, y = afy + r + 0.48,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  
  # Male — above "φ_adult M" (at y = amy - r - 0.23 ≈ −1.75)
  annotate("text",
           x = amx, y = amy - r - 0.48,
           label = "italic('JFM')", parse = TRUE,
           size = moci_sz, color = col_jfm),
  
  # ── Fecundity arc: OND ────────────────────────────────────────────────────
  # Below "f_t (fecundity)" which is at y = −0.10
  annotate("text",
           x = (px + afx)/2 + 0.30, y = -0.32,
           label = "italic('OND')", parse = TRUE,
           size = moci_sz, color = col_ond)
)

# ── Shared theme ──────────────────────────────────────────────────────────────
shared_theme <- list(
  theme_void(base_size = 13),
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin     = margin(20, 20, 30, 20),
    plot.title      = element_text(face = "bold", size = 14, hjust = 0.5,
                                   margin = margin(b = 6)),
    plot.subtitle   = element_text(size = 9.5, hjust = 0.5, color = "grey30",
                                   margin = margin(b = 10))
  ),
  coord_equal(xlim = c(0.0, 7.10), ylim = c(-2.40, 2.65)),
  scale_x_continuous(expand = c(0, 0)),
  scale_y_continuous(expand = c(0, 0))
)

# =============================================================================
# VERSION A — Marin IPM v3.3
# =============================================================================

sub_marin <- paste0(
  "\u03c6pup: pup survival  |  \u03c6juv: juvenile survival (stasis and advancement)  |  ",
  "\u03c6adult F: adult female  |  \u03c6adult M = \u03c6adult F \u2212 \u03b4adult\n",
  "Italic season labels: ",
  "OND = Oct\u2013Dec  \u00b7  ",
  "AMJ = Apr\u2013Jun  \u00b7  ",
  "JFM = Jan\u2013Mar  ",
  "(MOCI effect applies equally to stasis and advancement fractions of \u03c6juv)"
)

life_cycle_marin <- ggplot() +
  base_layers +
  moci_labels +
  shared_theme +
  labs(
    title    = "Harbor Seal Marin IPM v3.3 \u2014 Stage-Structured Life Cycle",
    subtitle = sub_marin
  )

ggsave("Output/Plots/life_cycle_marin_v3.3.jpeg",
       life_cycle_marin,
       width = 32, height = 22, units = "cm", dpi = 300)

cat("Saved: Output/Plots/life_cycle_marin_v3.3.jpeg\n")

# =============================================================================
# VERSION B — Regional IPM
# Identical structure; adds county-type modifier note on fecundity arc.
# =============================================================================

regional_extras <- list(
  # County modifier text below the OND fecundity label
  annotate("text",
           x = (px + afx)/2 + 0.30, y = -0.54,
           label = paste0("italic('+ county modifier'~",
                          "(delta[Marin]*', '*delta[Bay]*', '*delta[SB])~",
                          "'[fecundity only]')"),
           parse = TRUE, size = 2.9, color = "grey45")
)

sub_regional <- paste0(
  "\u03c6pup: pup survival  |  \u03c6juv: juvenile survival (stasis and advancement)  |  ",
  "\u03c6adult F: adult female  |  \u03c6adult M = \u03c6adult F \u2212 \u03b4adult\n",
  "Italic season labels: OND = Oct\u2013Dec  \u00b7  AMJ = Apr\u2013Jun  \u00b7  JFM = Jan\u2013Mar  ",
  "  |  Fecundity additionally modulated by county-type MOCI modifiers"
)

life_cycle_regional <- ggplot() +
  base_layers +
  moci_labels +
  regional_extras +
  shared_theme +
  labs(
    title    = "Regional Harbor Seal IPM \u2014 Stage-Structured Life Cycle",
    subtitle = sub_regional
  )

ggsave("Output/Plots/life_cycle_regional.jpeg",
       life_cycle_regional,
       width = 32, height = 22, units = "cm", dpi = 300)

cat("Saved: Output/Plots/life_cycle_regional.jpeg\n")
cat("Done.\n")
