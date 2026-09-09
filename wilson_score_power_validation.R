# wilson_score_power_validation.R
#
# Produces two figures:
#
#   Figure 1  manuscript_power_curves.png
#     Publication-ready analytic power curves (Wilson score CLT formula).
#     Intended to replace Picture1 in the editorial manuscript.
#
#   Figure 2  clt_approximation_check.png
#     Two-panel diagnostic:
#       Upper — CLT (solid), exact binomial (dashed), simulation (dots)
#       Lower — CLT error = CLT power − exact power, in percentage points
#
# Hypothesis: H0: p <= p0  vs  H1: p > p0 (one-sided), alpha = 0.05
# p0 = 0.90,  p1 in {0.94, 0.95, 0.96, 0.97, 0.98, 0.99}
#
# Key distinction between methods:
#
#   CLT analytic — approximates P(X >= c | Bin(N, p1)) using the normal
#     distribution for p-hat. Fast and closed-form; can overstate power
#     at small N, especially when p0 is far from 0.5.
#
#   Exact power of Wilson test — P(X >= c | Bin(N, p1)) via pbinom.
#     Same critical count c as the CLT test; no approximation. This is
#     the true power of the test the CLT formula is trying to compute.
#
#   Simulation — empirical rejection rate from B = 10 000 replicates,
#     each testing whether the 90% Wilson CI lower bound exceeds p0.
#     Should converge to the exact power as B -> infinity.
#
# The CLT approximation is reliable when N*(1-p0) >= 10 (the rarer event —
# failures — satisfies the CLT condition). For p0 = 0.90 this means N >= 100.
# Below this threshold the CLT can substantially overstate power.
#
# References: Fleiss, Levin & Paik (2003); Lachin (1981)
# Dependencies: ggplot2, dplyr, tidyr, patchwork

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

set.seed(2026)


# ── 1. Parameters ─────────────────────────────────────────────────────────────

p0     <- 0.90
p1_vec <- seq(0.94, 0.99, by = 0.01)
alpha  <- 0.05
B      <- 10000   # simulation replicates (Figure 2 only)


# ── 2. Core functions ─────────────────────────────────────────────────────────

# Critical count: Wilson score test rejects H0 iff X >= c.
# Applying this rule to your data involves no approximation.
wilson_crit <- function(p0, N, alpha) {
  ceiling(N * p0 + qnorm(1 - alpha) * sqrt(N * p0 * (1 - p0)))
}

# CLT analytic power — approximates P(X >= c | Bin(N, p1)) via the normal
# distribution. Same formula as proportion-calculator/index.html.
wilson_power_clt <- function(p0, p1, N, alpha) {
  za <- qnorm(1 - alpha)
  s0 <- sqrt(p0 * (1 - p0))
  s1 <- sqrt(p1 * (1 - p1))
  pnorm(((p1 - p0) * sqrt(N) - za * s0) / s1)
}

# Exact power of the Wilson test — evaluates P(X >= c | Bin(N, p1)) directly.
# c is the same critical count as above. Returns NA when c > N (test cannot
# reject any outcome, so power is 0; CLT does not detect this).
wilson_power_exact <- function(p0, p1, N, alpha) {
  c <- wilson_crit(p0, N, alpha)
  if (c > N) return(0)
  1 - pbinom(c - 1, N, p1)
}

# Analytic sample size: smallest N achieving target power under CLT formula.
wilson_n <- function(p0, p1, alpha, target_power) {
  za <- qnorm(1 - alpha)
  zb <- qnorm(target_power)
  s0 <- sqrt(p0 * (1 - p0))
  s1 <- sqrt(p1 * (1 - p1))
  ceiling(((za * s0 + zb * s1) / (p1 - p0))^2)
}

# Simulation: reject H0 iff lower bound of 90% Wilson CI > p0.
# The lower bound of a (1-2*alpha) two-sided Wilson CI equals the lower bound
# of a one-sided (1-alpha) Wilson CI, implementing a one-sided alpha test.
wilson_power_sim <- function(p0, p1, N, alpha, B) {
  z  <- qnorm(1 - alpha)
  x  <- rbinom(B, N, p1)
  pt <- (x + z^2 / 2) / (N + z^2)
  lb <- pt - z * sqrt(pt * (1 - pt) / (N + z^2))
  mean(lb > p0)
}


# ── 3. Built-in reference check ───────────────────────────────────────────────

cat("── Reference check (p0=0.90, p1=0.95, alpha=0.05, target=0.80) ─────────\n")
ref_N <- wilson_n(0.90, 0.95, 0.05, 0.80)
cat(sprintf("  N (CLT)         = %d   (expected 184)\n", ref_N))
cat(sprintf("  Power CLT       = %.4f (expected 0.8017)\n",
            wilson_power_clt(0.90, 0.95, ref_N, 0.05)))
cat(sprintf("  Power exact     = %.4f (exact binomial power of Wilson test)\n",
            wilson_power_exact(0.90, 0.95, ref_N, 0.05)))
cat(sprintf("  Critical count  = %d   (expected 173)\n\n",
            wilson_crit(0.90, ref_N, 0.05)))


# ── 4. Sample size table ──────────────────────────────────────────────────────

n_table <- expand.grid(p1 = p1_vec, target_power = c(0.80, 0.85, 0.90)) |>
  as_tibble() |>
  mutate(
    N          = mapply(wilson_n, p0, p1, alpha, target_power),
    power_clt  = round(mapply(wilson_power_clt,   p0, p1, N, alpha), 4),
    power_exact= round(mapply(Vectorize(wilson_power_exact), p0, p1, N, alpha), 4),
    crit       = mapply(wilson_crit, p0, N, alpha)
  )

cat("── Sample sizes: Wilson score CLT (p0 = 0.90, alpha = 0.05) ────────────\n")
print(n_table, n = Inf)
cat("\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Manuscript power curves (CLT analytic only)
# ══════════════════════════════════════════════════════════════════════════════

N_fig1 <- seq(10, 500, by = 2)

fig1_df <- expand.grid(p1 = p1_vec, N = N_fig1) |>
  as_tibble() |>
  mutate(power = mapply(wilson_power_clt, p0, p1, N, alpha)) |>
  filter(power >= 0.70, power <= 1)

fig1 <- ggplot(fig1_df, aes(x = N, y = power, color = factor(p1))) +
  geom_line(linewidth = 1) +
  geom_hline(
    yintercept = c(0.80, 0.85, 0.90),
    linetype = "dashed", color = "grey50", linewidth = 0.35
  ) +
  annotate("text",
           x     = 498,
           y     = c(0.80, 0.85, 0.90) + 0.012,
           label = c("80%", "85%", "90%"),
           hjust = 1, size = 3, color = "grey40") +
  scale_y_continuous(
    limits = c(0.70, 1.00),
    breaks = seq(0.70, 1.00, by = 0.05),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_continuous(
    limits = c(10, 500),
    breaks = seq(0, 500, by = 100),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  labs(
    title    = "H₀: p ≤ 0.90  •  H₁: p > 0.90",
    subtitle = paste0(
      "Wilson score test, one-sided α = ", alpha,
      "  —  power from analytic CLT formula"
    ),
    x     = "Sample size (N)",
    y     = "Power",
    color = "Expected\nproportion (p₁)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(color = "grey40", size = 10, hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position  = "right"
  )

ggsave("manuscript_power_curves.png", fig1, width = 9, height = 6, dpi = 150)
cat("Figure 1 saved: manuscript_power_curves.png\n\n")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — CLT approximation check
# ══════════════════════════════════════════════════════════════════════════════

# The CLT approximation is reliable when N*(1-p0) >= 10 (the standard rule
# for the failure count to satisfy CLT conditions). For p0 = 0.90: N >= 100.
# Below this, and especially where the critical count c exceeds N (test cannot
# reject any outcome), the CLT can substantially overstate the true power.

n_reliable <- ceiling(10 / (1 - p0))   # 100 for p0 = 0.90

# Minimum N where the test can actually reject (c <= N)
N_all <- 10:500
n_min_useful <- min(N_all[mapply(wilson_crit, p0, N_all, alpha) <= N_all])

cat(sprintf(
  "── CLT reliability thresholds (p0 = %.2f, alpha = %.2f) ──────────────────\n",
  p0, alpha
))
cat(sprintf("  n_min_useful (c <= N):          N >= %d\n", n_min_useful))
cat(sprintf("  n_reliable (N*(1-p0) >= 10):    N >= %d\n\n", n_reliable))

# Dense N grid for smooth CLT and exact lines
N_dense <- seq(10, 500, by = 4)

# Coarser grid for simulation (slow)
N_sim <- c(seq(30, 100, by = 20), seq(125, 500, by = 25))

# ── CLT and exact over the dense grid ─────────────────────────────────────

fig2_curves <- expand.grid(p1 = p1_vec, N = N_dense) |>
  as_tibble() |>
  mutate(
    clt   = mapply(wilson_power_clt,              p0, p1, N, alpha),
    exact = mapply(Vectorize(wilson_power_exact),  p0, p1, N, alpha)
  )

# ── Simulation over the coarse grid ───────────────────────────────────────

n_cells <- length(p1_vec) * length(N_sim)
cat(sprintf(
  "── Running simulations: %d cells × B = %s replicates ────────────────────\n",
  n_cells, format(B, big.mark = ",")
))

fig2_sim <- expand.grid(p1 = p1_vec, N = N_sim) |>
  as_tibble() |>
  mutate(
    sim = mapply(wilson_power_sim, p0, p1, N, alpha, B),
    se  = sqrt(sim * (1 - sim) / B)
  )

cat("Done.\n\n")

# ── Analytic vs simulation at N_80 ────────────────────────────────────────

n80 <- n_table |>
  filter(target_power == 0.80) |>
  select(p1, N80 = N, power_clt)

comparison <- n80 |>
  mutate(
    power_exact = round(mapply(Vectorize(wilson_power_exact), p0, p1, N80, alpha), 4),
    power_sim   = round(mapply(wilson_power_sim, p0, p1, N80, alpha, B), 4),
    clt_vs_exact_pp = round(100 * (power_clt - power_exact), 1),
    sim_vs_exact_pp = round(100 * (power_sim - power_exact),  1)
  )

cat("── Analytic vs exact vs simulation at N for 80% power ──────────────────\n")
cat("clt_vs_exact_pp and sim_vs_exact_pp in percentage points (+ = overstates)\n\n")
print(comparison, n = Inf)
cat("\n")

# ── Upper panel: power vs N, three methods ────────────────────────────────

panel_a <- ggplot() +
  # CLT analytic (solid lines)
  geom_line(
    data    = fig2_curves,
    mapping = aes(x = N, y = clt, color = factor(p1)),
    linewidth = 0.85, linetype = "solid"
  ) +
  # Exact binomial power (dashed lines)
  geom_line(
    data    = fig2_curves,
    mapping = aes(x = N, y = exact, color = factor(p1)),
    linewidth = 0.75, linetype = "22"
  ) +
  # Simulation (points + 95% CI error bars)
  geom_linerange(
    data    = fig2_sim,
    mapping = aes(x = N, ymin = sim - 1.96 * se, ymax = sim + 1.96 * se,
                  color = factor(p1)),
    linewidth = 0.5, alpha = 0.45
  ) +
  geom_point(
    data    = fig2_sim,
    mapping = aes(x = N, y = sim, color = factor(p1)),
    size = 1.5, shape = 16, alpha = 0.75
  ) +
  # Reliability boundary
  geom_vline(
    xintercept = n_reliable, linetype = "dotted",
    color = "grey30", linewidth = 0.5
  ) +
  annotate("text",
           x = n_reliable + 4, y = 0.35,
           label = paste0("N = ", n_reliable, "\n(CLT reliable\nthreshold)"),
           hjust = 0, size = 2.7, color = "grey30", lineheight = 0.9) +
  # Method legend via annotate
  annotate("segment",
           x = 360, xend = 390, y = 0.13, yend = 0.13,
           color = "grey20", linewidth = 0.85, linetype = "solid") +
  annotate("segment",
           x = 360, xend = 390, y = 0.09, yend = 0.09,
           color = "grey20", linewidth = 0.75, linetype = "22") +
  annotate("point",  x = 375, y = 0.05, color = "grey20", size = 1.5) +
  annotate("text",
           x = c(395, 395, 395), y = c(0.13, 0.09, 0.05),
           label = c("CLT analytic", "Exact (pbinom)", "Simulation ± 95% CI"),
           hjust = 0, size = 2.8, color = "grey20") +
  scale_y_continuous(
    limits = c(0, 1),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_x_continuous(breaks = seq(0, 500, by = 100)) +
  labs(
    title    = "CLT analytic vs exact binomial power vs simulation",
    subtitle = paste0(
      "Solid = CLT formula | Dashed = exact P(X ≥ c | Bin(N, p₁)) | ",
      "Dots = Wilson CI simulation (B = ", format(B, big.mark = ","), ") ± 95% CI"
    ),
    x     = NULL,
    y     = "Power",
    color = "p₁"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "grey40", size = 9),
    panel.grid.minor = element_blank(),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )

# ── Lower panel: CLT approximation error ──────────────────────────────────

error_df <- fig2_curves |>
  mutate(error_pp = 100 * (clt - exact))   # positive = CLT overstates

panel_b <- ggplot(error_df, aes(x = N, y = error_pp, color = factor(p1))) +
  annotate("rect",
           xmin = n_reliable, xmax = 500, ymin = -Inf, ymax = Inf,
           fill = "#e4f2ef", alpha = 0.55) +
  annotate("text",
           x = n_reliable + 4, y = max(error_df$error_pp) * 0.88,
           label = paste0("N ≥ ", n_reliable, ": CLT\nreliable"),
           hjust = 0, size = 2.7, color = "#0F766E", lineheight = 0.9) +
  geom_vline(
    xintercept = n_reliable, linetype = "dotted",
    color = "grey30", linewidth = 0.5
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  scale_y_continuous(
    labels = function(x) paste0(ifelse(x > 0, "+", ""), round(x, 0), " pp")
  ) +
  scale_x_continuous(breaks = seq(0, 500, by = 100)) +
  labs(
    subtitle = paste0(
      "CLT error = CLT power − exact power (pp). ",
      "Positive = CLT overstates. Large errors at small N arise when ",
      "the critical count c exceeds N and the test cannot reject any outcome."
    ),
    x     = "Sample size (N)",
    y     = "Error (pp)",
    color = "p₁"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.subtitle    = element_text(color = "grey40", size = 9),
    panel.grid.minor = element_blank()
  )

# ── Combine and save ───────────────────────────────────────────────────────

fig2 <- panel_a / panel_b +
  plot_layout(heights = c(2, 1), guides = "collect") &
  theme(legend.position = "right")

ggsave("clt_approximation_check.png", fig2, width = 10, height = 8, dpi = 150)
cat("Figure 2 saved: clt_approximation_check.png\n")
