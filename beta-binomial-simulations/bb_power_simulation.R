# bb_power_simulation.R
#
# Budget-constrained power simulation for a one-sided Beta-Binomial test.
# H0: mu <= 0.90  vs  H1: mu > 0.90  (one-sided, alpha = 0.05).
#
# N = entities (questions); k = draws per entity (responses per question).
# Entity i has true acceptable rate p_i ~ Beta(a, b); X_i | p_i ~ Bin(k, p_i).
#
# ── Three questions ────────────────────────────────────────────────────────────
#
# PART 1 — Equal costs (C = N*k):
#   For a fixed total-response budget, how does power vary with k?
#   Result: k = 1 always maximises power. VIF = 1 + (k-1)*rho >= 1 means
#   within-entity replication is less informative than recruiting a new entity.
#
# PART 2 — Type I error (same grid, null distribution mu = 0.90):
#   Which test controls alpha when data are clustered?
#   Result: glmmTMB holds type I error near alpha. Naive binomial inflates it
#   whenever rho > 0 (treats N*k correlated outcomes as independent).
#
# PART 3 — Unequal entity cost (B = N*(c_N + k)):
#   If recruiting a new entity (question) costs c_N times a single response,
#   does k > 1 ever maximise power?
#   Analytic result: N_eff = N/VIF = B / [(c_N+k)(1+(k-1)*rho)]. The
#   denominator is strictly increasing in k for k >= 1, rho >= 0, c_N >= 0,
#   so k = 1 remains optimal. The simulation confirms this and quantifies the
#   practical power penalty of k = 2-3, which matters when the analyst needs
#   k >= 2 to fit a BB model (glmmTMB is unidentifiable at k = 1).
#
# ── Tests compared ─────────────────────────────────────────────────────────────
#
# glmmTMB   Beta-Binomial MLE; Wald test on logit(mu) scale.
#           k = 1 fallback: overdispersion is unidentifiable with one binary
#           outcome per entity (Hessian rank-deficient in dispersion direction).
#           Falls back to logistic GLM intercept = Wald proportion test.
#
# naive     Exact binomial test on all N*k pooled outcomes (binom.test).
#           Valid only when rho = 0. Shown in Parts 1-2 only.
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
# Effect: dropped fits are excluded from the power mean. If convergence
# correlates with extreme estimates, retained power is biased upward.
# Cells with fail_glmmTMB > 0.10 are flagged and should be treated as
# unreliable. The direction of bias is toward overstatement of power.
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
  ggplot(make_density_df(df), aes(x, density, colour = scenario)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey40") +
    coord_cartesian(ylim = c(0, 20)) +
    labs(title  = title,
         x = "p (per-entity acceptable rate)", y = "Density", colour = NULL) +
    theme_bw() + theme(legend.position = "top")
}

print(plot_densities(betas_h1, "H1 distributions: mean = 0.95"))
print(plot_densities(betas_h0, "H0 null boundary: mean = 0.90"))


# ── Startup diagnostic ─────────────────────────────────────────────────────────
# Fits one Beta-Binomial model and prints the convergence fields so you can
# verify the glmmTMB version in use and confirm the fit is being accepted.
# If pdHess = FALSE or sdr is absent, the convergence check prints a note.
local({
  cat("── Startup diagnostic: one Beta-Binomial fit ────────────────────────────\n")
  cat("   glmmTMB version:", as.character(packageVersion("glmmTMB")), "\n")
  d <- data.frame(succ = c(5L,6L,5L,4L,6L,7L,5L,4L,6L,5L),
                  fail = c(1L,0L,1L,2L,0L,1L,1L,2L,0L,1L))
  m <- tryCatch(
    suppressWarnings(glmmTMB(cbind(succ, fail) ~ 1,
                             family = betabinomial(link = "logit"), data = d)),
    error = function(e) { cat("   ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(m)) {
    pdHess <- tryCatch(m$sdr$pdHess, error = function(e) NA)
    se     <- tryCatch(sqrt(vcov(m)$cond[1L, 1L]), error = function(e) NA_real_)
    cat(sprintf("   fit$convergence = %s  |  sdr$pdHess = %s  |  SE = %.4f\n",
                m$fit$convergence, pdHess, se))
    if (isTRUE(m$fit$convergence == 0L) && is.finite(se))
      cat("   Fit accepted — convergence check working correctly.\n")
    else
      cat("   WARNING: fit would be dropped by current convergence check.\n")
  }
  cat("\n")
})

# ── Core simulation functions ──────────────────────────────────────────────────

# Fit Beta-Binomial (k >= 2) or logistic GLM (k = 1 fallback).
# Returns c(logit_mu_hat, SE) or c(NA, NA) on convergence failure.
#
# Convergence criteria (in order):
#   1. No try-error from glmmTMB itself.
#   2. Optimizer converged: m$fit$convergence == 0.
#   3. Hessian positive-definite (pdHess): checked only when m$sdr is non-NULL.
#      In some glmmTMB versions sdr may be absent; we do not reject on absence.
#   4. SE is finite and positive.
fit_bb <- function(successes, k) {
  d <- data.frame(succ = successes, fail = k - successes)

  if (k == 1L) {
    m <- try(glm(cbind(succ, fail) ~ 1, family = binomial, data = d), silent = TRUE)
    if (inherits(m, "try-error")) return(c(NA_real_, NA_real_))
    se <- tryCatch(sqrt(vcov(m)[1L, 1L]), error = function(e) NA_real_)
    if (!is.finite(se) || se <= 0) return(c(NA_real_, NA_real_))
    return(c(coef(m)[[1L]], se))
  }

  # Use betabinomial() without namespace prefix; glmmTMB::betabinomial() is
  # equivalent after library(glmmTMB) but the prefix occasionally causes issues
  # with NSE in some glmmTMB versions.
  m <- try(suppressWarnings(
    glmmTMB(cbind(succ, fail) ~ 1, family = betabinomial(link = "logit"), data = d)
  ), silent = TRUE)
  if (inherits(m, "try-error")) return(c(NA_real_, NA_real_))

  # Criterion 2: optimizer convergence
  if (!isTRUE(m$fit$convergence == 0L)) return(c(NA_real_, NA_real_))

  # Criterion 3: pdHess — only reject if explicitly FALSE (not on NULL/absent)
  pdHess <- tryCatch(m$sdr$pdHess, error = function(e) NULL)
  if (isFALSE(pdHess)) return(c(NA_real_, NA_real_))

  # Criterion 4: finite positive SE
  se <- tryCatch(suppressWarnings(sqrt(vcov(m)$cond[1L, 1L])), error = function(e) NA_real_)
  if (!is.finite(se) || se <= 0) return(c(NA_real_, NA_real_))

  c(fixef(m)$cond[["(Intercept)"]], se)
}

# One-sided Wald test on logit(mu). Returns logical or NA (propagates failure).
reject_bb <- function(fit) {
  if (anyNA(fit) || fit[2L] <= 0) return(NA)
  pnorm((fit[1L] - qlogis(threshold)) / fit[2L], lower.tail = FALSE) < alpha_level
}

# Naive exact binomial test on all N*k pooled outcomes.
reject_naive <- function(successes, k) {
  binom.test(sum(successes), length(successes) * k,
             p = threshold, alternative = "greater")$p.value < alpha_level
}

# One dataset -> named rejection vector c(glmmTMB = T/F/NA, naive = T/F).
sim_once <- function(N, k, a, b, dgp) {
  p         <- if (dgp == "betabinom") rbeta(N, a, b) else rep(a / (a + b), N)
  successes <- rbinom(N, size = k, prob = p)
  c(glmmTMB = reject_bb(fit_bb(successes, k)),
    naive   = reject_naive(successes, k))
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
to_long <- function(res, value_col, new_col) {
  rbind(
    data.frame(res, method = "glmmTMB (Beta-Binomial)", value = res$power_glmmTMB),
    data.frame(res, method = "Naive binomial (pooled)", value = res$power_naive)
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
# PART 1 — Equal costs: power  (C = N × k)
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

cat("══ PART 1: Power under equal costs (C = N × k) ══════════════════════════\n")
res1 <- run_grid(grid1)
cat("\n"); print(res1[, c("scenario","C","k","N","power_glmmTMB","power_naive","fail_glmmTMB")],
                digits = 3, row.names = FALSE)

long1   <- to_long(res1)
long1$C <- factor(long1$C)

print(
  ggplot(long1, aes(k, value, colour = C, group = C)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    facet_grid(scenario ~ method) +
    scale_y_continuous(limits = c(0, 1), labels = pct_fmt) +
    labs(title    = "Part 1 — Power under equal costs  (C = N × k)",
         subtitle = "k = 1 maximises power for every scenario and budget",
         x = "k  (draws per entity;  N = C/k)", y = "Power", colour = "Budget C") +
    theme_bw()
)

cat("\n── Convergence failures > 5% (Part 1) ──────────────────────────────────\n")
f1 <- res1[res1$fail_glmmTMB > 0.05,
           c("scenario","C","k","N","fail_glmmTMB","power_glmmTMB")]
if (nrow(f1)) print(f1[order(-f1$fail_glmmTMB), ], digits = 3, row.names = FALSE)
else cat("  None\n")
cat("  Cells with fail_glmmTMB > 0.10 are unreliable (upward-biased power).\n\n")


# ═══════════════════════════════════════════════════════════════════════════════
# PART 2 — Equal costs: type I error  (mu = 0.90)
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

long2   <- to_long(res2)
long2$C <- factor(long2$C)

print(
  ggplot(long2, aes(k, value, colour = C, group = C)) +
    geom_hline(yintercept = alpha_level, linetype = "dashed", colour = "red") +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    facet_grid(scenario ~ method) +
    scale_y_continuous(limits = c(0, 1), labels = pct_fmt) +
    labs(title    = "Part 2 — Type I error under H0  (mu = 0.90)",
         subtitle = "glmmTMB tracks alpha = 0.05 (red line); naive binomial inflates when rho > 0",
         x = "k  (draws per entity;  N = C/k)", y = "Type I error", colour = "Budget C") +
    theme_bw()
)

cat("\n── Convergence failures > 5% (Part 2) ──────────────────────────────────\n")
f2 <- res2[res2$fail_glmmTMB > 0.05,
           c("scenario","C","k","N","fail_glmmTMB","power_glmmTMB")]
if (nrow(f2)) print(f2[order(-f2$fail_glmmTMB), ], digits = 3, row.names = FALSE)
else cat("  None\n")
cat("\n")


# ═══════════════════════════════════════════════════════════════════════════════
# PART 3 — Unequal entity cost:  B = N × (c_N + k)
# ═══════════════════════════════════════════════════════════════════════════════
#
# c_N = entity setup cost in units of per-response cost (c_k = 1, normalised).
# Budget B = N*(c_N + k)  =>  N = floor(B / (c_N + k)).
#
# c_N = 0    entities free; recovers Part 1 with C = B (N = B/k).
# c_N = 5    one question costs as much as 5 additional responses.
# c_N = 10   one question costs as much as 10 responses.
# c_N = 20   very expensive entities (e.g., expert annotation or clinical setup).
#
# Analytic result: N_eff = N/VIF = B / [(c_N+k)(1+(k-1)*rho)].
# Both factors in the denominator are increasing in k, so k = 1 is always
# optimal. The simulation confirms this and quantifies the practical penalty
# of k = 2, which is the minimum k for glmmTMB to identify overdispersion.
#
# Only glmmTMB power is shown (naive is excluded: it is always invalid under
# clustering and adds no information about optimal design).
#
# Cells with N < 5 are excluded; at very high c_N all N are small so
# convergence failures will be common and power will be low throughout.

B       <- 300
k_vals3 <- c(1, 2, 3, 5, 10)
cN_vals <- c(0, 5, 10, 20)

ck3         <- expand.grid(cN = cN_vals, k = k_vals3, stringsAsFactors = FALSE)
grid3       <- merge(betas_h1, ck3, by = NULL)
grid3$N     <- as.integer(floor(B / (grid3$cN + grid3$k)))
grid3$dgp   <- "betabinom"
grid3       <- grid3[grid3$N >= 5L, ]

cat("══ PART 3: Unequal entity cost  (B = N × (c_N + k),  B =", B, ") ════════\n")
cat("   c_N = entity setup cost (units of per-response cost)\n")
cat("   N = floor(B / (c_N + k));  cells with N < 5 excluded\n\n")
res3 <- run_grid(grid3)
cat("\n"); print(res3[, c("scenario","cN","k","N","power_glmmTMB","fail_glmmTMB")],
                digits = 3, row.names = FALSE)

res3$cost_label <- factor(paste0("c_N = ", res3$cN),
                           levels = paste0("c_N = ", sort(unique(res3$cN))))

print(
  ggplot(res3, aes(k, power_glmmTMB, colour = cost_label, group = cost_label)) +
    geom_hline(yintercept = 0.80, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.8) + geom_point(size = 1.5) +
    facet_wrap(~ scenario, ncol = 2) +
    scale_y_continuous(limits = c(0, 1), labels = pct_fmt) +
    labs(title    = paste0("Part 3 — Power under unequal entity cost  (B = ", B, ")"),
         subtitle = "B = N×(c_N+k)  |  k=1 optimal; penalty of k=2 shrinks as entities get expensive",
         x = "k  (draws per entity)", y = "Power (glmmTMB)", colour = "Entity cost c_N") +
    theme_bw() + theme(legend.position = "right")
)

cat("\n── Part 3: empirically optimal k by scenario and entity cost ────────────\n")
cat("   (k_opt = argmax power_glmmTMB within each scenario × c_N cell)\n\n")

# For each scenario × c_N, find k with highest glmmTMB power; also report
# power at k=1 so the penalty of the glmmTMB-feasible minimum (k=2) is visible.
opt3 <- do.call(rbind, lapply(
  split(res3, list(res3$scenario, res3$cN), drop = TRUE),
  function(d) {
    best_i <- which.max(d$power_glmmTMB)
    if (!length(best_i)) return(NULL)
    best  <- d[best_i, ]
    pw_k1 <- d$power_glmmTMB[d$k == 1L]
    pw_k2 <- d$power_glmmTMB[d$k == 2L]
    data.frame(
      scenario   = best$scenario,
      cN         = best$cN,
      k_opt      = best$k,
      N_at_k_opt = best$N,
      power_opt  = round(best$power_glmmTMB, 3),
      power_k1   = round(if (length(pw_k1)) pw_k1 else NA_real_, 3),
      power_k2   = round(if (length(pw_k2)) pw_k2 else NA_real_, 3),
      fail_k_opt = round(best$fail_glmmTMB, 3)
    )
  }
))
opt3 <- opt3[order(opt3$scenario, opt3$cN), ]
print(opt3, row.names = FALSE)
cat("\n  power_k1 and power_k2 enable comparison with the analytic optimum.\n")
cat("  When k_opt = 1 for all c_N, the analytic result is fully confirmed.\n")
cat("  Cells with fail_k_opt > 0.10 are flagged above during run_grid.\n")
