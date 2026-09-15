# bb_power_simulation.R
#
# Budget-constrained power simulation for a one-sided Beta-Binomial test.
# H0: mu <= 0.90  vs  H1: mu > 0.90  (one-sided, alpha = 0.05).
#
# N = entities (questions); k = draws per entity (responses per question).
# Entity i has true acceptable rate p_i ~ Beta(a, b); X_i | p_i ~ Bin(k, p_i).
#
# ── Cost model (one model, two cases) ─────────────────────────────────────────
#
# Every response costs one unit. Each query additionally incurs a fixed
# overhead c (in units of one response) for obtaining the query itself.
# Total cost of a design with N queries and k responses per query:
#     C = N * (c + k).
# Case 1 (Parts 1-2): c = 0, so C = N*k — the budget buys responses only.
# Case 2 (Part 3):    c > 0 — each new query carries an overhead that
#                     additional responses to an existing query avoid.
#
# PART 1 — Case 1 power (c = 0; C = N*k):
#   For a fixed budget, how does power vary with k?
#   Result: k = 1 always maximises power. VIF = 1 + (k-1)*rho >= 1 means an
#   additional response to an existing query is less informative than a
#   response to a new query.
#
# PART 2 — Case 1 type I error (same grid, null distribution mu = 0.90):
#   Which test controls alpha when data are clustered?
#   Result: glmmTMB holds type I error near alpha. The naive Wilson score
#   test inflates it when rho > 0 and k > 1 (treats N*k correlated outcomes
#   as independent); at k = 1 it is exactly valid for any rho.
#
# PART 3 — Case 2 power (c > 0; C = N*(c + k), C = 1000):
#   Does k > 1 ever maximise power once queries carry an overhead?
#   Result (simulation): for low-rho scenarios (Unimodal, J-shape) k = 10-20
#   is best at every c > 0 examined, and the advantage grows with c; for the
#   U-shaped scenarios k = 1 stays best or nearly best even at c = 20.
#
# ── Tests compared ─────────────────────────────────────────────────────────────
#
# glmmTMB   Beta-Binomial MLE (k >= 2); one-sided Wald test on logit(mu).
#           k = 1: the Beta-Binomial reduces exactly to the Binomial (rho is
#           unidentifiable from one binary draw per entity), so the Wilson
#           score test below is used. (A Wald logistic-GLM test at k = 1 has
#           zero power for N < ~60 because p_hat = 1 gives a degenerate logit;
#           this produced power = 0 in the c = 20 cells of Part 3.)
#
# naive     Wilson score test on all M = N*k pooled outcomes — the test
#           recommended in the manuscript and the calculator:
#             prop.test(X, M, p = p0, alternative = "greater", correct = FALSE)
#           (equivalently X >= ceiling(M*p0 + z_alpha*sqrt(M*p0*(1-p0)))).
#           Valid only when rho = 0 or k = 1. Shown in Parts 1-2 only.
#
# ── glmmTMB convergence failures ──────────────────────────────────────────────
#
# A fit is dropped (NA) when any of:
#   (a) Hessian not positive-definite: m$sdr$pdHess != TRUE
#   (b) Optimiser did not converge:    m$fit$convergence != 0
#   (c) SE is non-finite or <= 0
#
# Common causes:
#   - Small N: overdispersion poorly identified; flat likelihood ridge in phi
#   - rho near 0 (J-shape, Unimodal): phi -> infinity at boundary; flat region
#   - rho near 1 (Strong-U): phi -> 0 at boundary; Hessian ill-conditioned
#
# Effect: dropped fits are excluded from the power mean. Near rho = 0 the
# dropped datasets are the ones showing NO overdispersion (phi -> infinity),
# i.e. those with the smallest SE and highest power, so retained power is
# biased DOWNWARD. In addition, the Wald test on the logit scale is markedly
# conservative when the effective sample size N_eff = N*k/VIF is small
# (< ~60) and p_hat is near 1 (Hauck-Donner effect). Together these produce
# the dip at k = 3 for c = 20 in Part 3 (that design has more effectively
# independent observations than k = 1, so the drop is not a loss of
# information). Cells with fail_glmmTMB > 0.10 or N_eff < 60 are flagged in the
# convergence summary (`unreliable`) and should not be over-interpreted.
#
# ── Outputs (written to out_dir) ──────────────────────────────────────────────
#
#   results_part1.csv / results_part2.csv / results_part3.csv   full grids
#   power_equal_costs.png, type1_equal_costs.png, power_unequal_costs.png
#   convergence_failures.csv   fail_glmmTMB for every cell of all three parts
#   convergence_failures.png   heatmaps of fail_glmmTMB (scenario x k, by design)
#
# Dependencies: glmmTMB, ggplot2
# Runtime: ~1-3 hours for n_sims = 1000. Set n_sims = 200 for a quick check.

library(glmmTMB)
library(ggplot2)

set.seed(88)


# ── Global parameters ──────────────────────────────────────────────────────────

threshold   <- 0.90
alpha_level <- 0.05
n_sims      <- 1000   # replications per simulation cell

# Output directory: works whether the working directory is the repo root or
# the beta-binomial-simulations folder itself.
out_dir <- if (dir.exists("beta-binomial-simulations")) "beta-binomial-simulations" else "."
save_fig <- function(name, plot, width, height)
  ggsave(file.path(out_dir, name), plot, width = width, height = height, dpi = 200)

# ── Plot styling (manuscript figures) ─────────────────────────────────────────
# Scenario strips ordered by intra-query correlation rho, labelled with the
# Appendix terminology (N = queries, k = responses per query).
scenario_levels <- c("Binomial (rho=0)", "Unimodal (rho=0.024)", "J-shape (rho=0.048)",
                     "Mild-U (rho=0.500)", "Strong-U (rho=0.833)")
scenario_labels <- c("Binomial\n(ρ = 0)", "Unimodal\n(ρ = 0.02)", "J-shape\n(ρ = 0.05)",
                     "Mild U-shape\n(ρ = 0.50)", "Strong U-shape\n(ρ = 0.83)")
ms_scenario <- function(x, oneline = FALSE) {
  lab <- if (oneline) gsub("\n", " ", scenario_labels) else scenario_labels
  factor(x, levels = scenario_levels, labels = lab)
}
theme_ms <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(plot.title       = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          strip.text       = element_text(face = "bold"),
          legend.position  = "right")
}


# ── Beta scenarios ─────────────────────────────────────────────────────────────

# H1: mean = 0.95, rho determined by concentration (a+b)
betas_h1 <- data.frame(
  scenario = c("J-shape (rho=0.048)",
               "Unimodal (rho=0.024)",
               "Mild-U (rho=0.500)",
               "Strong-U (rho=0.833)"),
  a = c(19,   38,   0.95, 0.19),
  b = c(1,    2,    0.05, 0.01),
  stringsAsFactors = FALSE
)

# H0 null boundary: mean = 0.90, same rho (identical a+b concentration)
betas_h0 <- data.frame(
  scenario = betas_h1$scenario,
  a = c(18,   36,   0.90, 0.18),
  b = c(2,    4,    0.10, 0.02),
  stringsAsFactors = FALSE
)

print_betas <- function(df, label) {
  cat(label, "\n")
  for (i in seq_len(nrow(df))) {
    a <- df$a[i]; b <- df$b[i]
    cat(sprintf("  %-26s  Beta(%-.2f, %-.2f)  mean=%.4f  rho=%.4f\n",
                df$scenario[i], a, b, a / (a + b), 1 / (a + b + 1)))
  }
  cat("\n")
}

print_betas(betas_h1, "── H1 scenarios (mean = 0.95) ───────────────────────────────────────────")
print_betas(betas_h0, "── H0 null boundary (mean = 0.90, same rho) ─────────────────────────────")


# ── Density plots ──────────────────────────────────────────────────────────────

xx <- seq(0.001, 0.999, length.out = 500)

make_density_df <- function(df) {
  do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    data.frame(scenario = df$scenario[i], x = xx,
               density  = dbeta(xx, df$a[i], df$b[i]))
  }))
}

plot_densities <- function(df, title) {
  d <- make_density_df(df); d$scenario_ms <- ms_scenario(d$scenario, oneline = TRUE)
  ggplot(d, aes(x, density, colour = scenario_ms)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey40") +
    coord_cartesian(ylim = c(0, 20)) +
    labs(title = title, x = "Query-wise acceptability probability (pᵢ)",
         y = "Density", colour = NULL) +
    guides(colour = guide_legend(nrow = 2)) +
    theme_ms() + theme(legend.position = "top")
}

print(plot_densities(betas_h1, "H1 distributions: mean = 0.95"))
print(plot_densities(betas_h0, "H0 null boundary: mean = 0.90"))
save_fig("beta_distributions.png", plot_densities(betas_h1, "Query-level acceptability distributions (p = 0.95)"),
         width = 8, height = 4.5)


# ── Core simulation functions ──────────────────────────────────────────────────

# Fit Beta-Binomial (k >= 2 only; at k = 1 sim_once uses the exact test).
# Returns c(logit_mu_hat, SE) or c(NA, NA) on convergence failure.
fit_bb <- function(successes, k) {
  d <- data.frame(succ = successes, fail = k - successes)

  m <- try(suppressWarnings(
    glmmTMB(cbind(succ, fail) ~ 1, family = glmmTMB::betabinomial(), data = d)
  ), silent = TRUE)
  if (inherits(m, "try-error")) return(c(NA_real_, NA_real_))

  converged <- isTRUE(m$sdr$pdHess) && isTRUE(m$fit$convergence == 0L)
  se        <- suppressWarnings(sqrt(vcov(m)$cond[1L, 1L]))
  if (!converged || !is.finite(se) || se <= 0) return(c(NA_real_, NA_real_))

  c(fixef(m)$cond[["(Intercept)"]], se)
}

# One-sided Wald test on logit(mu). Returns logical or NA (propagates failure).
reject_bb <- function(fit) {
  if (anyNA(fit) || fit[2L] <= 0) return(NA)
  pnorm((fit[1L] - qlogis(threshold)) / fit[2L], lower.tail = FALSE) < alpha_level
}

# Naive Wilson score test on all M = N*k pooled outcomes (the manuscript's
# recommended test): base-R prop.test() without continuity correction, i.e.
# the one-sided Wilson score test — reject iff the lower one-sided Wilson
# bound exceeds p0. Equivalent to the critical-count form used in
# proportion-calculator/index.html, X >= ceiling(M*p0 + z_alpha*sqrt(M*p0*q0)).
reject_naive <- function(successes, k) {
  prop.test(sum(successes), length(successes) * k, p = threshold,
            alternative = "greater", correct = FALSE)$p.value < alpha_level
}

# One dataset -> named rejection vector c(glmmTMB = T/F/NA, naive = T/F).
# At k = 1 the Beta-Binomial is exactly Binomial, so both columns use the
# Wilson score test (see header: the Wald GLM is degenerate at small N).
sim_once <- function(N, k, a, b, dgp) {
  p         <- if (dgp == "betabinom") rbeta(N, a, b) else rep(a / (a + b), N)
  successes <- rbinom(N, size = k, prob = p)
  naive     <- reject_naive(successes, k)
  c(glmmTMB = if (k == 1L) naive else reject_bb(fit_bb(successes, k)),
    naive   = naive)
}

# n_sims replications -> list(power, fail_glmmTMB).
run_sim <- function(N, k, a, b, dgp) {
  out <- replicate(n_sims, sim_once(N, k, a, b, dgp))
  list(power        = rowMeans(out, na.rm = TRUE),
       fail_glmmTMB = mean(is.na(out["glmmTMB", ])))
}

# Run a prepared grid (must have columns: scenario, a, b, N, k, dgp).
# Appends power_glmmTMB, power_naive, fail_glmmTMB. Extra grid columns pass through.
# Convergence failures > 10% trigger a console warning.
run_grid <- function(grid) {
  n   <- nrow(grid)
  out <- do.call(rbind, lapply(seq_len(n), function(i) {
    row <- grid[i, ]
    cat(sprintf("  [%3d/%d]  N=%4d  k=%2d  %s\n", i, n, row$N, row$k, row$scenario))
    res <- run_sim(row$N, row$k, row$a, row$b, row$dgp)
    if (res$fail_glmmTMB > 0.10)
      message(sprintf("  !! fail_glmmTMB = %.0f%% — power estimate may be biased upward",
                      100 * res$fail_glmmTMB))
    data.frame(power_glmmTMB = unname(res$power["glmmTMB"]),
               power_naive   = unname(res$power["naive"]),
               fail_glmmTMB  = res$fail_glmmTMB,
               row.names     = NULL)
  }))
  cbind(grid, out)
}

pct_fmt <- function(x) paste0(round(100 * x), "%")

# Pivot power_glmmTMB / power_naive to long format for faceted plots.
to_long <- function(res) {
  out <- rbind(
    data.frame(res, method = "GLMM (Beta Binomial)", value = res$power_glmmTMB),
    data.frame(res, method = "Wilson score test",    value = res$power_naive)
  )
  out$method      <- factor(out$method, levels = c("GLMM (Beta Binomial)", "Wilson score test"))
  out$scenario_ms <- ms_scenario(out$scenario)
  out
}


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1 — Case 1 (c = 0): power  (C = N × k)
# ═══════════════════════════════════════════════════════════════════════════════

C_vals <- c(200, 300, 400)
k_vals <- c(1, 3, 5, 10, 20)
ck     <- expand.grid(C = C_vals, k = k_vals, stringsAsFactors = FALSE)

# Beta-Binomial DGP for all four scenarios
g1_bb      <- merge(betas_h1, ck, by = NULL)
g1_bb$N    <- as.integer(g1_bb$C / g1_bb$k)
g1_bb$dgp  <- "betabinom"

# Binomial control (rho=0, p fixed at 0.95); reuse first scenario's a/b for DGP
g1_bin          <- merge(betas_h1[1L, ], ck, by = NULL)
g1_bin$N        <- as.integer(g1_bin$C / g1_bin$k)
g1_bin$dgp      <- "binomial"
g1_bin$scenario <- "Binomial (rho=0)"

grid1 <- rbind(g1_bb, g1_bin)

cat("══ PART 1: Case 1 — no per-query overhead (c = 0; C = N × k): power ═════\n")
res1 <- run_grid(grid1)
cat("\n"); print(res1[, c("scenario","C","k","N","power_glmmTMB","power_naive","fail_glmmTMB")],
                digits = 3, row.names = FALSE)

write.csv(res1, file.path(out_dir, "results_part1.csv"), row.names = FALSE)

long1   <- to_long(res1)
long1$C <- factor(long1$C)

p1 <- ggplot(long1, aes(k, value, colour = C, group = C)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  facet_grid(scenario_ms ~ method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = pct_fmt) +
  scale_x_continuous(breaks = k_vals) +
  labs(title = "Power  —  no per-query overhead  (c = 0;  C = N × k)",
       x = "Responses per query (k)", y = "Power", colour = "Budget (C)") +
  theme_ms()
print(p1)
save_fig("power_equal_costs.png", p1, width = 8, height = 10)

cat("\n── Convergence failures > 5% (Part 1) ──────────────────────────────────\n")
f1 <- res1[res1$fail_glmmTMB > 0.05,
           c("scenario","C","k","N","fail_glmmTMB","power_glmmTMB")]
if (nrow(f1)){
  print(f1[order(-f1$fail_glmmTMB), ], digits = 3, row.names = FALSE)
} else {cat("  None\n")}
cat("  Cells with fail_glmmTMB > 0.10 are unreliable (upward-biased power).\n\n")


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2 — Case 1 (c = 0): type I error  (mu = 0.90)
# ═══════════════════════════════════════════════════════════════════════════════

g2_bb      <- merge(betas_h0, ck, by = NULL)
g2_bb$N    <- as.integer(g2_bb$C / g2_bb$k)
g2_bb$dgp  <- "betabinom"

g2_bin          <- merge(betas_h0[1L, ], ck, by = NULL)
g2_bin$N        <- as.integer(g2_bin$C / g2_bin$k)
g2_bin$dgp      <- "binomial"
g2_bin$scenario <- "Binomial (rho=0)"

grid2 <- rbind(g2_bb, g2_bin)

cat("══ PART 2: Type I error under H0  (mu = 0.90) ════════════════════════════\n")
res2 <- run_grid(grid2)
cat("\n"); print(res2[, c("scenario","C","k","N","power_glmmTMB","power_naive","fail_glmmTMB")],
                digits = 3, row.names = FALSE)

write.csv(res2, file.path(out_dir, "results_part2.csv"), row.names = FALSE)

long2   <- to_long(res2)
long2$C <- factor(long2$C)

p2 <- ggplot(long2, aes(k, value, colour = C, group = C)) +
  geom_hline(yintercept = alpha_level, linetype = "dashed", colour = "red") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  facet_grid(scenario_ms ~ method) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1), labels = pct_fmt) +
  scale_x_continuous(breaks = k_vals) +
  labs(title = "Type I error rate  —  no per-query overhead  (c = 0;  C = N × k)",
       x = "Responses per query (k)", y = "Type I error rate", colour = "Budget (C)") +
  theme_ms()
print(p2)
save_fig("type1_equal_costs.png", p2, width = 8, height = 10)

cat("\n── Convergence failures > 5% (Part 2) ──────────────────────────────────\n")
f2 <- res2[res2$fail_glmmTMB > 0.05,
           c("scenario","C","k","N","fail_glmmTMB","power_glmmTMB")]
if (nrow(f2)) {print(f2[order(-f2$fail_glmmTMB), ], digits = 3, row.names = FALSE)
} else {cat("  None\n")}
cat("\n")


# ═══════════════════════════════════════════════════════════════════════════════
# PART 3 — Case 2 (c > 0): power  (C = N × (c + k))
# ═══════════════════════════════════════════════════════════════════════════════
#
# Cost model (see header): C = N*(c + k), with c the per-query overhead in
# units of one response. c = 5 means obtaining one new query costs as much
# as generating and rating five responses. N = floor(C / (c + k)).
#
# c = 0   reproduces Case 1 at a larger budget (C = N*k).
# c = 3, 5, 10, 20   increasing overhead per query.
#
# Only glmmTMB power is shown (naive is excluded: it is invalid under
# clustering for k > 1 and adds no information about optimal design).
# The empirically best k per scenario x c is tabulated in `opt3` below.
#
# Cells with N < 5 are excluded; at very high c all N are small so
# convergence failures will be common and power will be low throughout.

C_total <- 1000
k_vals3 <- c(1, 3, 5, 10, 20)
cN_vals <- c(0, 3, 5, 10, 20)

ck3         <- expand.grid(cN = cN_vals, k = k_vals3, stringsAsFactors = FALSE)
grid3       <- merge(betas_h1, ck3, by = NULL)
grid3$N     <- as.integer(floor(C_total / (grid3$cN + grid3$k)))
grid3$dgp   <- "betabinom"
grid3       <- grid3[grid3$N >= 5L, ]

cat("══ PART 3: Case 2 — per-query overhead  (C = N × (c + k),  C =", C_total, ") ═══\n")
cat("   c = per-query overhead (units of one response)\n")
cat("   N = floor(C / (c + k));  cells with N < 5 excluded\n\n")
res3 <- run_grid(grid3)
cat("\n"); print(res3[, c("scenario","cN","k","N","power_glmmTMB","fail_glmmTMB")],
                digits = 3, row.names = FALSE)

write.csv(res3, file.path(out_dir, "results_part3.csv"), row.names = FALSE)

res3$cost_label  <- factor(paste0("c = ", res3$cN),
                            levels = paste0("c = ", sort(unique(res3$cN))))
res3$scenario_ms <- ms_scenario(res3$scenario, oneline = TRUE)

p3 <- ggplot(res3, aes(k, power_glmmTMB, colour = cost_label, group = cost_label)) +
  geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 0.8) + geom_point(size = 1.6) +
  facet_wrap(~ scenario_ms, ncol = 2) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = pct_fmt) +
  scale_x_continuous(breaks = k_vals3) +
  labs(title = paste0("Power  —  per-query overhead  (c > 0;  C = N × (c + k) = ", C_total, ")"),
       x = "Responses per query (k)", y = "Power (GLMM)", colour = "Per-query\noverhead (c)") +
  theme_ms()
print(p3)
save_fig("power_unequal_costs.png", p3, width = 8, height = 7)

cat("\n── Part 3: empirically best k by scenario and per-query overhead ────────\n")
cat("   k_best = argmax power_glmmTMB within each scenario × c cell\n\n")

opt3 <- do.call(rbind, lapply(
  split(res3, list(res3$scenario, res3$cN), drop = TRUE),
  function(d) {
    best_i <- which.max(d$power_glmmTMB)
    if (!length(best_i)) return(NULL)
    best  <- d[best_i, ]
    pw_k1 <- d$power_glmmTMB[d$k == 1L]
    data.frame(
      scenario    = best$scenario,
      c           = best$cN,
      k_best      = best$k,
      N_at_k_best = best$N,
      power_best  = round(best$power_glmmTMB, 3),
      power_k1    = round(if (length(pw_k1)) pw_k1 else NA_real_, 3),
      fail_k_best = round(best$fail_glmmTMB, 3)
    )
  }
))
opt3 <- opt3[order(opt3$scenario, opt3$c), ]
print(opt3, row.names = FALSE)
write.csv(opt3, file.path(out_dir, "results_part3_optimal_k.csv"), row.names = FALSE)
cat("\n  Cells with fail_k_best > 0.10 are flagged above during run_grid.\n")


# ═══════════════════════════════════════════════════════════════════════════════
# CONVERGENCE FAILURE SUMMARY  (all three parts)
# ═══════════════════════════════════════════════════════════════════════════════
#
# fail_glmmTMB = share of the n_sims replications whose Beta-Binomial fit was
# dropped (see header). Written as a CSV table and plotted as heatmaps
# (scenario x k, one panel per part x design). A cell is `unreliable` when
# fail_glmmTMB > 0.10 or N_eff = N*k/VIF < 60 (small-sample Wald-test
# conservatism); such cells are starred in the heatmap.

conv_cols <- c("scenario", "a", "b", "k", "N", "fail_glmmTMB", "power_glmmTMB")
conv <- rbind(
  data.frame(part = "Case 1: power",        design = paste0("C = ", res1$C),  res1[, conv_cols]),
  data.frame(part = "Case 1: type I error", design = paste0("C = ", res2$C),  res2[, conv_cols]),
  data.frame(part = "Case 2: power",        design = paste0("c = ", res3$cN), res3[, conv_cols])
)
conv$rho   <- ifelse(grepl("^Binomial", conv$scenario), 0, 1 / (conv$a + conv$b + 1))
conv$N_eff <- round(conv$N * conv$k / (1 + (conv$k - 1) * conv$rho))
conv$a <- NULL; conv$b <- NULL
conv$flag_fail  <- conv$fail_glmmTMB > 0.10
conv$flag_small <- conv$k > 1L & conv$N_eff < 60   # k = 1 uses the Wilson test: no Wald issue
conv$unreliable <- conv$flag_fail | conv$flag_small
write.csv(conv, file.path(out_dir, "convergence_failures.csv"), row.names = FALSE)

cat("\n══ Convergence failure summary ══════════════════════════════════════════\n")
cat(sprintf("  %d of %d cells have fail_glmmTMB > 10%%; %d have N_eff < 60; %d unreliable overall\n\n",
            sum(conv$flag_fail), nrow(conv), sum(conv$flag_small), sum(conv$unreliable)))
conv_by_scn <- aggregate(fail_glmmTMB ~ part + scenario, conv[conv$k > 1L, ],
                         function(x) c(mean = mean(x), max = max(x)))
conv_by_scn <- do.call(data.frame, conv_by_scn)
names(conv_by_scn)[3:4] <- c("mean_fail", "max_fail")
print(conv_by_scn[order(conv_by_scn$part, -conv_by_scn$mean_fail), ],
      digits = 2, row.names = FALSE)
cat("  (k = 1 cells excluded: exact binomial test, no fitting)\n")

panel_lab    <- paste0(conv$part, "\n", conv$design)   # two lines: avoids strip clipping
conv$panel   <- factor(panel_lab, levels = unique(panel_lab))
conv$k_f     <- factor(conv$k, levels = sort(unique(conv$k)))
conv$label   <- paste0(sprintf("%.0f%%", 100 * conv$fail_glmmTMB), ifelse(conv$flag_small, "*", ""))

conv$scenario_ms <- ms_scenario(conv$scenario, oneline = TRUE)
conv$scenario_ms <- factor(conv$scenario_ms, levels = rev(levels(conv$scenario_ms)))  # top-to-bottom by rho

p_conv <- ggplot(conv, aes(k_f, scenario_ms, fill = fail_glmmTMB)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = label, colour = fail_glmmTMB > 0.35), size = 2.8) +
  scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "white"), guide = "none") +
  scale_fill_gradient(low = "#f7fbff", high = "#b2182b", limits = c(0, 1),
                      labels = pct_fmt, name = "GLMM fits\ndiscarded") +
  facet_wrap(~ panel, ncol = 3, scales = "free_y") +
  labs(title = "GLMM convergence failures by simulation condition",
       x = "Responses per query (k)", y = NULL,
       caption = "* fewer than 60 effective observations (N × k / VIF)") +
  theme_ms(base_size = 11) +
  theme(panel.grid = element_blank(), strip.text = element_text(size = 8.5, face = "bold"))
print(p_conv)
save_fig("convergence_failures.png", p_conv, width = 12, height = 10)

cat(sprintf("\nOutputs written to '%s'.\n", normalizePath(out_dir)))
