# Analytic Calculator vs. Simulation: Sources of Discrepancy

*Companion document to METHODS.md. Explains why the browser calculator's power
and type I error predictions can differ by several percentage points from Monte
Carlo results. For k ≥ 2 the simulation fits a Beta-Binomial model with glmmTMB
and tests on the logit scale. For k = 1, glmmTMB cannot identify the
overdispersion parameter (one observation per cluster), so the simulation falls
back to a simple one-sample proportion test (OLS intercept) instead (as in
`v5_beta_binomial.R`).*

---

## 1. Summary of the Four Structural Differences

| Source | Calculator | Simulation (k ≥ 2, glmmTMB) | Simulation (k = 1) |
|---|---|---|---|
| **Test statistic** | Score test on proportion scale | Wald test on logit(μ) scale | t-test on proportion (OLS intercept) |
| **ρ (ICC)** | Specified by user, treated as known | Estimated jointly with μ from data | Not applicable (VIF = 1 always) |
| **Power formula** | Analytic CLT approximation | Empirical rejection rate over 1 000 simulated datasets | Empirical rejection rate |
| **What is being tested** | H₀: μ ≤ p₀ with user-fixed VIF | H₀: logit(μ) ≤ logit(p₀) with MLE-estimated overdispersion | H₀: μ ≤ p₀ via simple proportion test |

**Note on the k = 1 column.** When k = 1 there is exactly one binary response per question, so glmmTMB cannot identify the overdispersion parameter (no within-cluster replication means the Hessian is rank-deficient or the optimisation converges to a boundary). The simulation therefore falls back to a simple one-sample proportion test — equivalent to OLS linear regression of the binary outcomes on an intercept — rather than glmmTMB. This means the k = 1 column is a clean comparison of the score test (calculator) against the Wald t-test on the proportion scale (simulation), with no confounding from ρ estimation or convergence failures.

None of these differences makes one approach "wrong." They answer slightly
different questions. The discrepancies that arise are predictable, directional,
and explainable from first principles.

---

## 2. Test Statistic Scale: Proportion vs. Logit

**The single largest contributor to the discrepancy.**

The calculator's Beta-Binomial score test rejects when

    T_score = (p̂ − p₀) / √(p₀q₀·VIF/N)  >  z_α,

evaluated on the proportion scale. The glmmTMB analysis fits the full
Beta-Binomial likelihood, then performs a Wald test on the logit scale:

    T_Wald = (logit(μ̂) − logit(p₀)) / SE(logit(μ̂))  >  z_α.

These two statistics are asymptotically equivalent but differ in finite samples.
The difference is amplified when p₀ is far from 0.5 — exactly the regime of LLM
evaluation (p₀ = 0.90, μ₁ = 0.95), where the logit transformation is
increasingly non-linear.

**Quantitative illustration (k=1, purely binomial, N=200).**

With k=1 there is one binary response per question, so glmmTMB cannot estimate
overdispersion and the simulation falls back to a simple one-sample proportion
test (OLS regression of binary outcomes on an intercept, equivalent to a Wald
t-test on p̂). There is therefore no ρ estimation to worry about, and the
comparison is clean:

| Quantity | Value |
|---|---|
| logit(0.95) − logit(0.90) | 2.944 − 2.197 = **0.747** |
| SE of logit(p̂) under H₁ at N=200 | 1/√(200×0.95×0.05) = **0.324** |
| Wald power (logit / t-test scale) | Φ(0.747/0.324 − 1.645) = Φ(0.66) ≈ **0.745** |
| Score test power (proportion scale) | Φ((0.05√200 − 1.645×0.3)/0.218) ≈ **0.837** |
| Simulated power (simple proportion test) | ≈ **0.77** |

The calculator predicts 0.837; the simulation lands near 0.77; the logit Wald
formula gives 0.745. The gap at k=1 is entirely the test-scale effect — the
simulation does not use glmmTMB here, so convergence and ρ estimation play no role.

**Direction:** Because the score test on the proportion scale uses the null SE
(√(p₀q₀/N)) as denominator while the Wald test on the logit scale uses the
estimated SE (which converges to the SE under the alternative), the score test
is somewhat more powerful for p₀ near 1. The calculator will therefore
consistently *overstate* achievable power compared to a glmmTMB analysis, and
the gap grows as p₀ → 1.

---

## 3. ICC Assumed Known vs. Simultaneously Estimated

**The calculator treats ρ as a user-specified planning constant.** Power and
type I error are computed at the exact ρ the researcher enters. The VIF is a
deterministic function: VIF = 1 + (k−1)ρ.

**glmmTMB estimates ρ from the same data used to estimate μ.** This adds three
complications:

1. **Estimation variance.** The MLE of ρ has its own sampling distribution,
   especially at small n. This inflates the variance of μ̂ beyond what the
   assumed-known formula predicts, reducing realized power below the analytic
   prediction.

2. **Convergence failures.** For small n (which happens at large k under a fixed
   budget C = N·k) or when true ρ is near 0 or near 1, the Beta-Binomial
   likelihood can have a flat ridge or be poorly identified. The R script drops
   non-converged fits (`isTRUE(m$sdr$pdHess)`), which creates selection bias: if
   convergence is more likely when the data look extreme, the retained fits will
   have a biased power estimate. This is the main explanation for the non-monotone
   (V-shaped) power curves visible in the glmmTMB column of the power plot at
   small N.

3. **Bias of the MLE for ρ near zero.** The Beta-Binomial overdispersion
   parameter is bounded below at 0. The MLE is biased away from zero in small
   samples, which means ρ is systematically overestimated, inflating the
   estimated VIF and deflating power. This makes glmmTMB conservative for the
   J-shape and Unimodal scenarios (ρ ≈ 0.024–0.048).

**Practical implication.** When specifying ρ in the calculator, use a value
you believe is realistic or conservative. The calculator gives the power of the
ideal test with known ρ; real glmmTMB analyses will be slightly less powerful
due to estimation uncertainty.

---

## 4. Analytic vs. Empirical Type I Error

The calculator's naive-binomial inflation formula (§4.3 of METHODS.md) is:

    true α_naive = 1 − Φ(z_α / √VIF).

This is a normal (score-test) approximation. The simulation uses the exact
binomial test (`binom.test`) on the pooled N·k observations. The two agree
closely when N·k is large enough for the CLT to apply, and diverge modestly at
small N·k.

**Verification against simulation Type I error plots.** For the strong-U case
(Beta(0.18, 0.02), ρ = 0.833):

| k | N (at C=300) | VIF | Calculator true α (naive) | Simulation (approx.) |
|---|---|---|---|---|
| 3 | 100 | 2.67 | 0.157 | ≈ 0.15–0.18 |
| 5 | 60 | 4.33 | 0.215 | ≈ 0.22–0.26 |
| 10 | 30 | 8.50 | 0.286 | ≈ 0.28–0.32 |
| 20 | 15 | 16.8 | 0.344 | ≈ 0.35–0.45 |

Agreement is good at moderate k; the simulation runs somewhat higher at k=20
where N is very small (15 questions, 300 total observations) and the CLT
approximation weakens. **The glmmTMB column, by contrast, stays near 0.05 across
all k and C values**, confirming that the BB model with estimated ρ properly
controls the type I error — it is the correctly-calibrated test.

---

## 5. Analytic Power Predictions at Fixed Budget C = N·k

The table below gives the calculator's BB score test power prediction for each
simulation scenario at budget C = 200, with the key parameters. Compare these
to the glmmTMB power curves in the simulation plots (green line, C = 200).

**p₀ = 0.90, μ₁ = 0.95, α = 0.05 (one-sided).**

### Beta(19, 1) — J-shape, ρ = 0.048

| k | N=C/k | VIF | BB score power | Naive stated power | Naive true α |
|---|---|---|---|---|---|
| 1 | 200 | 1.000 | 0.837 | 0.837 | 0.050 |
| 3 | 67 | 1.095 | 0.318 | 0.354 | 0.058 |
| 5 | 40 | 1.190 | 0.175 | 0.228 | 0.066 |
| 10 | 20 | 1.429 | 0.080 | 0.150 | 0.084 |
| 20 | 10 | 1.905 | 0.041 | 0.132 | 0.117 |

The simulation's glmmTMB power is generally 3–12 pp lower than the BB score
test column, primarily due to the test-scale difference (§2) and ρ estimation
uncertainty for this low-ρ scenario. The non-monotone dip at k=3–5 in glmmTMB
is a convergence artifact (§3, point 2).

### Beta(38, 2) — Unimodal, ρ = 0.024

| k | N=C/k | VIF | BB score power | Naive stated power | Naive true α |
|---|---|---|---|---|---|
| 1 | 200 | 1.000 | 0.837 | 0.837 | 0.050 |
| 3 | 67 | 1.049 | 0.332 | 0.351 | 0.054 |
| 5 | 40 | 1.098 | 0.190 | 0.219 | 0.058 |
| 10 | 20 | 1.220 | 0.091 | 0.131 | 0.068 |
| 20 | 10 | 1.463 | 0.048 | 0.102 | 0.087 |

Very small ρ means VIF ≈ 1 for small k; the single-proportion (naive) test and
BB test give nearly the same power, and the type I error inflation is minor.
glmmTMB simulation values track close to the BB score test column for k ≤ 5,
with larger deviations at k=10–20 where n is very small.

### Beta(0.95, 0.05) — Mild-U, ρ = 0.513

| k | N=C/k | VIF | BB score power | Naive stated power | Naive true α |
|---|---|---|---|---|---|
| 1 | 200 | 1.000 | 0.837 | 0.837 | 0.050 |
| 3 | 67 | 2.026 | 0.172 | 0.392 | 0.124 |
| 5 | 40 | 3.051 | 0.076 | 0.321 | 0.173 |
| 10 | 20 | 5.615 | 0.034 | 0.301 | 0.244 |
| 20 | 10 | 10.74 | 0.021 | 0.319 | 0.308 |

The BB score test power drops steeply with k because the large ρ makes
within-question replication nearly useless — each additional response per
question adds little new information. The naive test's stated power (0.39 at
k=3) is three times the true BB power (0.17) because it ignores that the
effective sample size is N/VIF, not N.

### Beta(0.19, 0.01) — Strong-U, ρ = 0.833

| k | N=C/k | VIF | BB score power | Naive stated power | Naive true α |
|---|---|---|---|---|---|
| 1 | 200 | 1.000 | 0.837 | 0.837 | 0.050 |
| 3 | 67 | 2.667 | 0.132 | 0.405 | 0.157 |
| 5 | 40 | 4.333 | 0.059 | 0.348 | 0.215 |
| 10 | 20 | 8.500 | 0.028 | 0.336 | 0.286 |
| 20 | 10 | 16.83 | 0.018 | 0.354 | 0.344 |

Notably, the glmmTMB simulation shows higher power here than the BB score test
predicts (e.g., ~0.30 vs 0.132 at k=3, C=200). This inversion occurs because
the bimodal Beta (most questions either very easy or very hard, with mean 0.95)
gives the MLE of μ a particularly clean signal: the logit-scale Wald test in
glmmTMB benefits from the fact that logit(0.95) is far from logit(0.90) and
the MLE concentrates sharply there. The score test on the proportion scale does
not exploit this asymmetry as effectively. The BB score test's CLT approximation
also breaks down more severely for bimodal ρ near 1 at small n.

---

## 6. What Each Tool Is Appropriate For

| Use case | Use the calculator | Use Monte Carlo / glmmTMB |
|---|---|---|
| Quick planning under assumed ρ | ✓ | |
| Understanding how VIF degrades effective N | ✓ | |
| Validating calibration of a specific test | | ✓ |
| Accounting for ρ estimation uncertainty | | ✓ |
| Budget-constrained optimization over k | ✓ (via manual scan) | ✓ |
| Non-equal k per question | | ✓ |
| Verifying calculator predictions | | ✓ |

The calculator is best understood as a planning tool for the **ideal scenario**:
the score test with known VIF gives an upper bound on what is achievable. Real
analyses with glmmTMB will be slightly less powerful. A conservative planning
rule is to target power 5–10 pp above the minimum acceptable level to absorb
the gap between the analytic score-test prediction and the empirical glmmTMB
Wald-test power.

---

## 7. The Binomial Control Scenario

The bottom rows of the simulation plots use a truly binomial data-generating
process (p fixed at 0.95 or 0.90, no question-level variation, ρ=0). With ρ=0,
VIF=1, and the calculator's BB formula collapses exactly to the standard
proportion test. Yet glmmTMB's power curve for the binomial scenario shows the
same V-shape dip at k=3–5 seen in the J-shape case. This cannot be a VIF
effect (VIF=1 always), confirming that the non-monotonicity is a
**convergence artifact** of fitting the Beta-Binomial model when the true ρ is
near zero: glmmTMB attempts to estimate a near-zero overdispersion from a small
number of questions, the Hessian is poorly conditioned, and many fits are
dropped. This is also visible in the `glmmTMB_fail` column of the printed
results table, which is highest at intermediate k values in these low-ρ
scenarios.

---

## 8. Practical Recommendation

When using the calculator for planning:

1. **At k=1**, glmmTMB cannot be used (overdispersion is unidentifiable with one
   observation per cluster). The simulation uses a simple one-sample proportion
   test instead. The only discrepancy vs. the calculator is therefore the
   test-scale difference (score on proportion scale vs. Wald t-test on proportion
   scale); expect realized power ~5–10 pp lower than the calculator's prediction.
   For your own k=1 analysis, use an exact binomial test or Wilson score test —
   not glmmTMB.

2. **At k≥3 with ρ > 0.1**, the VIF correction dominates, and the calculator's
   type I error formula for the naive test is accurate. The BB score test power
   is a reasonable upper bound for what glmmTMB will achieve.

3. **For ρ near 0 (J-shape, Unimodal)**, glmmTMB power at moderate k is
   dominated by convergence behavior and n; avoid k > C/30 (fewer than ~30
   questions) where the overdispersion is poorly identified.

4. **For ρ near 1 (strong-U, bimodal)**, the logit-scale Wald test in glmmTMB
   can substantially exceed the calculator's score-test power prediction at
   moderate k. Run a targeted simulation using the specific ρ and k values to
   verify.

5. **In all cases**, set the target power in the calculator at least 10 pp above
   the minimum acceptable value to provide a buffer for estimation uncertainty
   and the test-scale gap.

---

*Document version: September 2026.*
