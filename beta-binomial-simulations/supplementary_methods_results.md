# Supplementary Note: How many prompts, and how many responses per prompt?

An LLM evaluation judges responses as *acceptable* or *not* and asks whether the
acceptable rate exceeds a threshold (here p₀ = 0.90). If the same prompt can be answered
several times, a study has two knobs: the number of distinct prompts **N** and the number
of responses per prompt **k**. This note explains how to set them and when the simple
one-proportion test used in the main text remains valid.

## Key takeaways

1. **Spend the budget on new prompts, not repeated responses.** For a fixed number of
   judged responses, one response per prompt (k = 1) gives the highest power, whatever the
   spread of prompt difficulty. Responses to the same prompt are partly redundant; a
   response on a new prompt never is.
2. **With one response per prompt, the simple test is exactly valid.** Every judged
   outcome is then an independent draw with probability equal to the mean acceptable rate,
   so the Wilson score test (and its sample-size formula) applies without any correction.
3. **If prompts are repeated, the simple test over-rejects — badly when prompts are
   polarised, mildly otherwise.** When most prompts are effectively always-pass or
   always-fail, the false-positive rate of the pooled test rises from 5% to 10–40% as k
   grows. When prompt difficulty clusters around the mean, it stays near 5% for k ≤ 5 and
   reaches only about 10% at k = 20.
4. **Repeat prompts only when new prompts are expensive and difficulty is homogeneous.**
   If creating a new prompt costs as much as *c* judged responses, the cost-optimal number
   of responses per prompt is k\* = √(c(1−ρ)/ρ). This exceeds 1 only when prompts are costly
   and the between-prompt correlation ρ is small.

## Why: responses to the same prompt are partly redundant

Each prompt *i* has its own acceptable rate *p_i*; a hard prompt tends to fail every time
and an easy one to pass every time. Writing *p_i* ~ Beta(a, b) and X_i | p_i ~ Binomial(k, p_i)
gives the Beta-Binomial model, summarised by the mean rate μ = a/(a+b) — the quantity the
test is about — and the between-prompt correlation ρ = 1/(a+b+1), which measures how
much prompts differ in difficulty [6].

Clustering inflates the variance of the pooled proportion by the *design effect*, or
variance inflation factor [1–3]:

    VIF = 1 + (k − 1)·ρ

so N·k judged responses carry only N·k / VIF responses' worth of information. For a
fixed total of C = N·k responses, N·k/VIF = C/VIF is largest when k = 1 (takeaway 1). With
k = 1, VIF = 1 for any ρ (takeaway 2). Applying the pooled test to clustered data uses a
standard error that is too small by √VIF, so its true false-positive rate is
1 − Φ(z_α/√VIF) instead of α (takeaway 3).

When a new prompt costs *c* responses, a budget B buys N = B/(c + k) prompts and the
effective number of independent responses is

    N_eff = B·k / [(c + k)(1 + (k − 1)ρ)],

which is maximised at k\* = √(c(1−ρ)/ρ) — the classical optimal cluster-size result from
two-stage sampling and cluster-randomised design [4, 5] (takeaway 4).

## How we checked this

We simulated evaluations under the Beta-Binomial model and compared two analyses:

- **Pooled Wilson score test** (the main-text test): all N·k outcomes treated as
  independent; reject H₀: μ ≤ 0.90 if X ≥ ⌈N·k·p₀ + z_α√(N·k·p₀(1−p₀))⌉, α = 0.05.
- **Beta-Binomial model**: maximum-likelihood fit (`glmmTMB`, R) with a one-sided Wald
  test on logit(μ); the reference analysis that accounts for ρ. At k = 1 the model reduces
  to the binomial, so both analyses coincide.

Five prompt-difficulty distributions were used, all with μ = 0.95 under the alternative
(μ = 0.90 for false-positive checks):

| Scenario | Beta(a, b) | ρ | Reading |
|---|---|---|---|
| Binomial control | all p_i = 0.95 | 0 | every prompt equally hard |
| Unimodal | Beta(38, 2) | 0.024 | difficulty clusters near the mean |
| J-shape | Beta(19, 1) | 0.048 | most prompts near 1, a tail of harder ones |
| Mild-U | Beta(0.95, 0.05) | 0.500 | most prompts always-pass or always-fail |
| Strong-U | Beta(0.19, 0.01) | 0.833 | almost all prompts always-pass or always-fail |

Designs: (Parts 1–2) fixed total responses C = N·k ∈ {200, 300, 400} with
k ∈ {1, 3, 5, 10, 20}; (Part 3) budget B = 1000 with prompt cost
c ∈ {0, 3, 5, 10, 20} and N = ⌊B/(c + k)⌋. 1000 replications per cell (Monte Carlo SE
≤ 0.016). Script and outputs: `bb_power_simulation.R`, `results_part*.csv`, figures
`power_equal_costs.png`, `type1_equal_costs.png`, `power_unequal_costs.png`.

## What we found

**Equal cost (Figure `power_equal_costs.png`).** k = 1 gave the highest power in every
scenario and budget. At C = 300 responses, Beta-Binomial power fell from 0.95 at k = 1 to
0.76 (Unimodal), 0.69 (J-shape), 0.27 (Mild-U) and 0.10 (Strong-U) at k = 20. In the
Binomial control (ρ = 0) the pooled test's power was flat in k, as VIF = 1 predicts.

**False-positive rate (Figure `type1_equal_costs.png`).** At k = 1 the pooled test stayed
close to its 5% level in every scenario (0.04–0.07). For k ≥ 3 it reached 0.11–0.38
(Mild-U) and 0.15–0.40 (Strong-U), versus at most 0.14 (J-shape), 0.09 (Unimodal) and
0.07 (Binomial). The Beta-Binomial model stayed within 0.01–0.10 throughout (highest at
N ≤ 20 prompts).

**Unequal cost (Figure `power_unequal_costs.png`).** The best k in the simulation tracked
the analytic optimum k\* = √(c(1−ρ)/ρ):

| Scenario (ρ) | c = 3 | c = 5 | c = 10 | c = 20 |
|---|---|---|---|---|
| Unimodal (0.024) | k\* 11 · best 10 | 14 · 10 | 20 · 20 | 28 · 20 |
| J-shape (0.048) | 7.7 · 10 | 10 · 10 | 14 · 10 | 20 · 20 |
| Mild-U (0.500) | 1.7 · 1 | 2.2 · 1 | 3.2 · 5 | 4.5 · 10 |
| Strong-U (0.833) | 0.8 · 1 | 1.0 · 1 | 1.4 · 1 | 2.0 · 1 |

*(k\* · empirically best k over the grid {1, 3, 5, 10, 20}; c in responses per prompt.
Where power is near 1, neighbouring k are tied within Monte Carlo error.)*
With expensive prompts (c = 20) and homogeneous difficulty, repeating each of 25 prompts
20 times gave power 0.91–0.94, versus 0.31–0.32 for 47 single-response prompts. With
polarised prompts (Strong-U) repeating never helped.

## Caveats

- The Beta-Binomial fit fails when ρ is near 0 or 1 (the correlation parameter sits on the
  boundary); such fits were dropped and power computed over the remainder. Failure rates
  per cell are in `convergence_failures.csv` / `convergence_failures.png`. Near ρ ≈ 0 the
  dropped datasets are those showing no overdispersion — the ones with the smallest
  standard errors — so retained power is, if anything, understated.
- **Small effective samples.** With fewer than about 60 effective responses
  (N·k/VIF < 60) the Beta-Binomial Wald test is markedly conservative (the logit-scale
  standard error balloons as the observed rate approaches 1) and up to 64% of fits fail.
  This is why, in Figure `power_unequal_costs.png`, the k = 3 point at c = 20 dips *below*
  k = 1 before rising again: the analytic power for those cells is monotone in k
  (e.g. J-shape: 0.25 → 0.59 → 0.76 → 0.89 → 0.93), and the k = 1 point uses the
  Wilson score test, which does not suffer from this. Cells flagged `unreliable` in
  `convergence_failures.csv` (failures > 10% or N·k/VIF < 60) should not be
  over-interpreted; they do not affect takeaways 1–4.
- The simulation fixes μ = 0.95 vs 0.90 and α = 0.05; absolute power depends on these,
  the ordering of designs does not.
- Planning formulas treat ρ as known. In practice estimate it from pilot data (fit a
  Beta-Binomial) and, when unsure, plan with a larger ρ.

## References

1. Kish, L. (1965). *Survey Sampling*. Wiley. (Design effect.)
2. Donner, A., Birkett, N., & Buck, C. (1981). Randomization by cluster: sample size
   requirements and analysis. *American Journal of Epidemiology*, 114(6), 906–914.
   https://doi.org/10.1093/oxfordjournals.aje.a113261
3. Donner, A., & Klar, N. (2000). *Design and Analysis of Cluster Randomization Trials in
   Health Research*. Arnold.
4. Cochran, W. G. (1977). *Sampling Techniques* (3rd ed.), Ch. 10 (two-stage sampling;
   optimal subsample size). Wiley.
5. Raudenbush, S. W. (1997). Statistical analysis and optimal design for cluster
   randomized trials. *Psychological Methods*, 2(2), 173–185.
   https://doi.org/10.1037/1082-989X.2.2.173
6. Ridout, M. S., Demétrio, C. G. B., & Firth, D. (1999). Estimating intraclass correlation
   for binary data. *Biometrics*, 55(1), 137–148.
   https://doi.org/10.1111/j.0006-341X.1999.00137.x
7. Wilson, E. B. (1927). Probable inference, the law of succession, and statistical
   inference. *Journal of the American Statistical Association*, 22(158), 209–212.
