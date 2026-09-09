# wilson_score_power_validation.R
#
# Purpose
# -------
# Generate power curves for a one-sided one-sample proportion test and
# validate the Wilson score CLT analytic formula against Monte Carlo
# simulation using the Wilson score interval test.
#
# Hypothesis
# ----------
#   H0: p <= p0   vs   H1: p > p0
#   p0 = 0.90  |  alpha = 0.05 (one-sided)  |  p1 in {0.94, ..., 0.99}
#
# Two approaches are compared
# ---------------------------
#   Analytic   — Wilson score CLT closed-form power/sample-size formula.
#                These are the same formulas used in the proportion-calculator
#                website (proportion-calculator/index.html).
#                Reference: Fleiss, Levin & Paik (2003); Lachin (1981).
#
#   Simulation — empirical rejection rate using the Wilson score interval.
#                For each simulated dataset, reject H0 iff the lower bound
#                of the 90% Wilson CI exceeds p0.  The lower bound of a
#                (1-2*alpha) = 90% two-sided CI equals the lower bound of
#                a one-sided (1-alpha) = 95% CI, implementing a one-sided
#                alpha = 0.05 test without approximation.
#
# The arcsine approximation (pwr::pwr.p.test) is also included to show
# why it underestimates N when p is far from 0.5.
#
# Dependencies: pwr, ggplot2, dplyr, tidyr

library(pwr)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(2026)


# ── 1.  Parameters ────────────────────────────────────────────────────────────

p0     <- 0.90
p1_vec <- seq(0.94, 0.99, by = 0.01)
alpha  <- 0.05
B      <- 10000   # simulation replicates per (p1, N) cell


# ── 2.  Analytic functions — Wilson score CLT ─────────────────────────────────
#
# These reproduce the JavaScript in proportion-calculator/index.html exactly.
#
# Power formula:
#   power = Phi( [(p1 - p0) * sqrt(N)  -  z_alpha * sqrt(p0*q0)] / sqrt(p1*q1) )
#
# Sample size formula:
#   N = ceiling( ((z_alpha*s0 + z_beta*s1) / (p1 - p0))^2 )
#   where s0 = sqrt(p0*q0), s1 = sqrt(p1*q1)
#
# Critical count (integer threshold for the test decision):
#   c = ceiling( N*p0 + z_alpha * sqrt(N) * s0 )
#   Reject H0 iff observed successes X >= c.
#   No approximation is involved in applying this rule to data.

wilson_power <- function(p0, p1, N, alpha) {
  za <- qnorm(1 - alpha)
  s0 <- sqrt(p0 * (1 - p0))
  s1 <- sqrt(p1 * (1 - p1))
  pnorm(((p1 - p0) * sqrt(N) - za * s0) / s1)
}

wilson_n <- function(p0, p1, alpha, target_power) {
  za <- qnorm(1 - alpha)
  zb <- qnorm(target_power)
  s0 <- sqrt(p0 * (1 - p0))
  s1 <- sqrt(p1 * (1 - p1))
  ceiling(((za * s0 + zb * s1) / (p1 - p0))^2)
}

wilson_crit <- function(p0, N, alpha) {
  za <- qnorm(1 - alpha)
  ceiling(N * p0 + za * sqrt(N) * sqrt(p0 * (1 - p0)))
}


# ── 3.  Simulation function — Wilson score interval test ──────────────────────
#
# The Wilson lower bound is computed directly from the closed-form expression
# rather than calling binom.confint(), which is equivalent but slower.
#
# For x successes in N trials, the lower bound of the (1-2*alpha) Wilson CI is:
#   p_tilde = (x + z^2/2) / (N + z^2)
#   N_tilde = N + z^2
#   lower   = p_tilde - z * sqrt(p_tilde*(1-p_tilde) / N_tilde)
# where z = qnorm(1 - alpha).
#
# This is identical to binom.confint(x, N, conf.level = 1-2*alpha, "wilson")$lower.

wilson_lower_ci <- function(x, N, alpha) {
  z  <- qnorm(1 - alpha)
  pt <- (x + z^2 / 2) / (N + z^2)
  nt <- N + z^2
  pt - z * sqrt(pt * (1 - pt) / nt)
}

wilson_power_sim <- function(p0, p1, N, alpha, B) {
  x  <- rbinom(B, N, p1)
  lb <- wilson_lower_ci(x, N, alpha)
  mean(lb > p0)
}


# ── 4.  Built-in validation at reference values ───────────────────────────────
#
# Cross-check against the values in proportion-calculator/index.html and
# METHODS.md (p0=0.90, p1=0.95, alpha=0.05, power target=0.80).

cat("── Internal validation (p0=0.90, p1=0.95, alpha=0.05) ──────────────────\n")
ref_N    <- wilson_n(0.90, 0.95, 0.05, 0.80)       # expected: 184
ref_pow  <- wilson_power(0.90, 0.95, ref_N, 0.05)  # expected: 0.8017
ref_crit <- wilson_crit(0.90, ref_N, 0.05)          # expected: 173
cat(sprintf("  N = %d  (expected 184)\n", ref_N))
cat(sprintf("  power at N = %.4f  (expected 0.8017)\n", ref_pow))
cat(sprintf("  critical count c = %d  (expected 173)\n\n", ref_crit))


# ── 5.  Sample size table ─────────────────────────────────────────────────────

targets <- c(0.80, 0.85, 0.90)

n_table <- expand.grid(p1 = p1_vec, target_power = targets) |>
  as_tibble() |>
  mutate(
    N_wilson    = mapply(wilson_n, p0, p1, alpha, target_power),
    power_check = round(mapply(wilson_power, p0, p1, N_wilson, alpha), 4),
    crit        = mapply(wilson_crit, p0, N_wilson, alpha),
    N_arcsine   = ceiling(mapply(function(p1, tp) {
      pwr.p.test(
        h = ES.h(p1, p0), sig.level = alpha,
        power = tp, alternative = "greater"
      )$n
    }, p1, target_power)),
    arcsine_actual_power = round(
      mapply(wilson_power, p0, p1, N_arcsine, alpha), 4
    )
  )

cat("── Sample sizes: Wilson CLT vs arcsine (pwr.p.test) ─────────────────────\n")
cat("p0 =", p0, "  alpha =", alpha, "(one-sided)\n")
cat("Note: 'arcsine_actual_power' is the Wilson CLT power at the arcsine N.\n",
    "      It is below the target, showing the arcsine method underestimates N.\n\n")
print(n_table, n = Inf)
cat("\n")


# ── 6.  Power curves — analytic ───────────────────────────────────────────────

N_range   <- seq(30, 600, by = 5)

analytic_df <- expand.grid(p1 = p1_vec, N = N_range) |>
  as_tibble() |>
  mutate(
    power  = mapply(wilson_power, p0, p1, N, alpha),
    source = "Analytic (Wilson CLT)"
  )


# ── 7.  Power curves — simulation ─────────────────────────────────────────────
# Simulate at a coarser N grid to keep runtime manageable (~minutes).

N_sim_grid <- seq(30, 600, by = 30)
n_cells    <- length(p1_vec) * length(N_sim_grid)

cat(sprintf(
  "── Running simulations: %d cells, B = %s replicates each ────────────────\n",
  n_cells, format(B, big.mark = ",")
))

sim_df <- expand.grid(p1 = p1_vec, N = N_sim_grid) |>
  as_tibble() |>
  mutate(
    power  = mapply(wilson_power_sim, p0, p1, N, alpha, B),
    source = "Simulation"
  )

cat("Done.\n\n")


# ── 8.  Analytic vs simulation at N for 80% power ─────────────────────────────

n80 <- n_table |>
  filter(target_power == 0.80) |>
  select(p1, N80 = N_wilson)

comparison <- n80 |>
  mutate(
    power_analytic = round(mapply(wilson_power, p0, p1, N80, alpha), 4),
    power_sim      = round(mapply(wilson_power_sim, p0, p1, N80, alpha, B), 4),
    diff_pp        = round(100 * (power_sim - power_analytic), 1)
  )

cat("── Analytic vs simulation at N for 80% power ────────────────────────────\n")
cat("B =", format(B, big.mark = ","), "replicates per cell\n")
cat("diff_pp = 100*(power_sim - power_analytic); positive = simulation above analytic\n\n")
print(comparison, n = Inf)
cat("\n")


# ── 9.  Plot ──────────────────────────────────────────────────────────────────

# Reference lines: one horizontal per target power
ref_lines <- data.frame(yintercept = c(0.80, 0.85, 0.90),
                        label = c("80%", "85%", "90%"))

plot_df <- bind_rows(analytic_df, sim_df) |>
  mutate(p1_label = factor(paste0("p₁ = ", p1)))

p <- ggplot(plot_df, aes(x = N, y = power, color = p1_label)) +
  # analytic curves
  geom_line(
    data      = filter(plot_df, source == "Analytic (Wilson CLT)"),
    linewidth = 0.9
  ) +
  # simulation points
  geom_point(
    data  = filter(plot_df, source == "Simulation"),
    size  = 1.8,
    shape = 16,
    alpha = 0.75
  ) +
  # target power reference lines
  geom_hline(
    data     = ref_lines,
    aes(yintercept = yintercept),
    linetype = "dashed",
    color    = "grey40",
    linewidth = 0.4,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ref_lines,
    aes(x = 600, y = yintercept + 0.025, label = label),
    hjust = 1, size = 3, color = "grey40",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = function(x) paste0(round(100 * x), "%"),
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  scale_x_continuous(
    limits = c(30, 600),
    expand = expansion(mult = 0.02)
  ) +
  labs(
    title    = "Wilson score test power: analytic formula vs. simulation",
    subtitle = paste0(
      "H₀: p ≤ ", p0,
      "  │  one-sided α = ", alpha,
      "  │  lines = Wilson CLT analytic, dots = simulation (B = ",
      format(B, big.mark = ","), ")"
    ),
    x     = "Sample size (N)",
    y     = "Power",
    color = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position  = "right",
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(color = "grey45", size = 10),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave(
  "wilson_power_validation.png",
  plot   = p,
  width  = 9,
  height = 6,
  dpi    = 150
)
cat("Plot saved to wilson_power_validation.png\n")
