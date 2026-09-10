# Statistical Methods: LLM Evaluation Power Calculator

*Methods documentation for the browser-based power calculator. All formulas are stated with sufficient precision that any reader can reproduce every output the calculator produces.*

---

## 1. Scope and Estimand

The calculator addresses a canonical evaluation paradigm: a set of LLM-generated responses is judged by human raters (or an automated rubric), each response receiving a binary outcome — acceptable (1) or not (0). The primary estimand is the true acceptable-response probability, tested against a pre-specified null value p₀.

Two design families are supported.

**Tab 1 — Independent Bernoulli trials.** Each of N responses is an independent draw. Three methods are provided: the exact binomial test, the score (Wilson) z-test, and the score z-test with continuity correction.

**Tab 2 — Clustered (Beta-Binomial) design.** Responses are grouped by question. Each question i has its own latent acceptability probability pᵢ drawn from a Beta distribution. Within a question, k responses are sampled and independently judged. This introduces intra-question correlation that must be accounted for in power calculations.

---

## 2. Notation

| Symbol | Meaning |
|--------|---------|
| N | Total number of responses (observations) |
| n | Number of questions (clusters); N = n·k |
| k | Responses per question |
| p₀ | Null acceptable-response probability |
| p₁, μ₁ | Alternative acceptable-response probability |
| δ | Effect size: δ = μ₁ − p₀ |
| α | One-sided type I error rate |
| β | Type II error rate; target power = 1 − β |
| zα | Upper-α quantile of N(0,1): Φ⁻¹(1 − α) |
| zβ | Upper-β quantile of N(0,1): Φ⁻¹(1 − β) |
| q = 1 − p | Complement probability |
| Φ(·) | Standard normal CDF |
| αB, βB | Shape parameters of the Beta prior on pᵢ |
| φ | Beta concentration: φ = αB + βB |
| μ | Beta mean: μ = αB / (αB + βB) |
| ρ | Intra-class correlation: ρ = 1 / (1 + φ) |
| VIF | Variance inflation factor: VIF = 1 + (k − 1)ρ |

All tests are one-sided (upper tail: H₁: p > p₀). The calculator assumes p₁ > p₀.

---

## 3. Tab 1 — One-Proportion Tests

### 3.1 Exact Binomial Test

Under the null, X ~ Binomial(N, p₀). The test rejects when X ≥ c, where the critical count c is the smallest integer such that

    P(X ≥ c | Bin(N, p₀)) ≤ α.

**Computation of the critical value.** To avoid numerical underflow at large N or extreme p, all probability mass function values are computed in log-space using a downward recurrence from k = N. The starting value is

    log P(X = N) = N · log(p₀).

The recurrence step uses the ratio of successive PMF values:

    P(X = k−1) / P(X = k) = [k / (N − k + 1)] · [(1 − p₀) / p₀],

so

    log P(X = k−1) = log P(X = k) + log[k / (N − k + 1)] + log[(1 − p₀) / p₀].

The tail probability P(X ≥ c) is then the sum of exp(log P(X = k)) for k = c, c+1, ..., N. The critical count c is the smallest c such that this sum does not exceed α.

**Achieved type I error.** The achieved α is P(X ≥ c | Bin(N, p₀)), which may be strictly less than the nominal α because the binomial is discrete.

**Power.** With the critical count c fixed, power is

    Power = P(X ≥ c | Bin(N, p₁)),

computed by the same log-space tail sum but with p₁ in place of p₀.

**Sample size.** N is found by scanning N = 1, 2, 3, ... and, for each N, computing c and then power, stopping at the first N for which power ≥ target.

---

### 3.2 Score (Wilson) Z-Test

The score test evaluates the standard error at the null parameter value p₀ rather than at the sample proportion p̂. This is equivalent (in the one-sided case) to asking whether the lower bound of the one-sided Wilson confidence interval exceeds p₀. The interval itself was introduced by Wilson (1927) [1]; the power and sample size formulas below are the standard analytic results derived from that score statistic — see Fleiss, Levin & Paik (2003) [4] §4.2 and Lachin (1981) [5] for derivations in the clinical-trial context.

**Closed-form sample size.** The standard score-test formula is

    N = ⌈ (zα · √(p₀q₀) + zβ · √(p₁q₁))² / (p₁ − p₀)² ⌉,

where qᵢ = 1 − pᵢ and ⌈·⌉ denotes the ceiling function.

**Power at a given N.** Given fixed N,

    Power = Φ( [(p₁ − p₀) · √N − zα · √(p₀q₀)] / √(p₁q₁) ).

**Critical count (for display).** The test rejects when the observed count X satisfies

    (X/N − p₀) / √(p₀q₀/N) > zα,

which rearranges to X > N·p₀ + zα · √N · √(p₀q₀). Because X is an integer, the effective rejection region is X ≥ c where

    c = ⌈ N·p₀ + zα · √N · √(p₀q₀) ⌉.

---

### 3.3 Score Z-Test with Continuity Correction

The continuity correction shifts the critical boundary by ½ to improve the normal approximation to the discrete binomial.

**Critical count.**

    c = ⌈ N·p₀ + zα · √N · √(p₀q₀) + ½ ⌉.

**Achieved α.**

    achieved α = 1 − Φ( (c − 0.5 − N·p₀) / (√N · √(p₀q₀)) ).

**Power.**

    Power = Φ( (N·p₁ − c + 0.5) / (√N · √(p₁q₁)) ).

The ±½ correction treats the discrete count X as if it were drawn from a continuous distribution on (X − ½, X + ½), so the rejection event X ≥ c corresponds to the continuous event Y > c − ½ under H₁.

**Sample size.** Because the continuity correction disrupts the closed-form relationship between N and power, sample size is found by searching upward from the uncorrected score-test N (minus 2, to allow for the CC making the test more conservative), stopping at the first N where the corrected power formula reaches the target.

---

## 4. Tab 2 — Beta-Binomial Clustered Design

### 4.1 Model and Variance Structure

Each question i (i = 1, ..., n) has a latent acceptability probability

    pᵢ ~ Beta(αB, βB),

and conditional on pᵢ the k responses to question i are

    Xᵢ | pᵢ ~ Binomial(k, pᵢ),  i = 1, ..., n,

independently across questions.

The marginal distribution of Xᵢ is Beta-Binomial(k, αB, βB). Its mean and variance are

    E[Xᵢ] = k·μ,
    Var(Xᵢ) = k·μ·(1−μ) · [1 + (k−1)·ρ],

where μ = αB / (αB + βB) is the population mean acceptability probability and ρ = 1/(1 + αB + βB) = 1/(1 + φ) is the intra-class correlation (ICC).

The ICC arises because two responses to the same question share the common latent probability pᵢ; their correlation is

    Corr(Xᵢⱼ, Xᵢⱼ') = ρ  for j ≠ j',

which follows directly from the law of total variance applied to the Beta-Binomial model.

The variance inflation factor relative to a simple binomial with the same N = n·k observations is

    VIF = 1 + (k − 1)·ρ.

When k = 1 (one response per question), VIF = 1 regardless of ρ, recovering the standard binomial; this is the correct edge-case behaviour because with only one observation per cluster there is no within-cluster replication over which correlation can accumulate.

### 4.2 Score Test for the Beta-Binomial Mean

The total count X = Σᵢ Xᵢ has mean N·μ and variance N·μ·(1−μ)·VIF under the Beta-Binomial model. A score test of H₀: μ = p₀ versus H₁: μ > p₀ evaluates the test statistic

    T = (p̂ − p₀) / √(p₀·q₀·VIF / N),

where p̂ = X/N. Under H₀, T is approximately N(0,1).

**Power.** Under H₁ with true mean μ₁, the test statistic has approximate distribution N(λ, 1) where the non-centrality parameter is

    λ = (μ₁ − p₀) / √(μ₁·q₁·VIF / N).

The rejection condition T > zα then gives

    Power = Φ( [(μ₁ − p₀)·√(N/VIF) − zα·√(p₀·q₀)] / √(μ₁·q₁) ).

**Critical count.** Rejection occurs when X ≥ c where

    c = ⌈ N·p₀ + zα · √(p₀·q₀) · √N · √VIF ⌉.

**Sample size.** Using the power formula, solving for N (with N = n·k and VIF = 1 + (k−1)ρ) yields the number of questions:

    n = ⌈ VIF · (zα·√(p₀·q₀) + zβ·√(μ₁·q₁))² / (δ² · k) ⌉,

and total responses N = n·k.

### 4.3 Inflation When Naive Binomial Test Is Applied to Beta-Binomial Data

A common error is to apply the standard one-proportion z-test — which assumes independent Bernoulli observations — to clustered data. Under the Beta-Binomial model, the naive test statistic

    T_naive = (p̂ − p₀) / √(p₀·q₀ / N)

is not N(0,1) under H₀; instead it is approximately N(0, VIF). This follows because Var(p̂) = Var(X)/N² = μ·(1−μ)·VIF/N under the Beta-Binomial, whereas the naive denominator assumes variance p₀·q₀/N.

**True type I error rate.** The naive test rejects at threshold zα, but the actual rejection probability under H₀ is

    true α = 1 − Φ(zα / √VIF).

Because VIF > 1 for k > 1 and ρ > 0, the actual α exceeds the nominal α.

**True power.** Under H₁ with true mean μ₁, the naive test statistic has approximate distribution N((μ₁−p₀)/√(μ₁·q₁·VIF/N), VIF). The true rejection probability is therefore

    true power = Φ( [(μ₁ − p₀)·√N − zα·√(p₀·q₀)] / [√(μ₁·q₁)·√VIF] ).

This exceeds the correctly-calibrated power, meaning the naive approach simultaneously understates the actual type I error and overstates the actual power relative to the stated α — a doubly misleading result.

---

## 5. Numerical Implementation

**Normal CDF (Φ).** Implemented via the complementary error function erfc using a Horner-form Chebyshev polynomial approximation (10 terms), giving absolute errors below 10⁻¹⁵ across the range of inputs used.

**Inverse normal (Φ⁻¹).** Implemented using the rational approximations of Acklam (a variant of the Beasley-Springer-Moro algorithm), followed by one Newton-Raphson refinement step. This achieves relative accuracy suitable for all displayed decimal places.

**Binomial tail probability.** Uses the log-space downward recurrence described in §3.1. This avoids underflow for large N or probabilities near 0 or 1. The implementation has been verified at N = 3886, p = 0.50 versus p = 0.52 (a regime where naive direct computation of PMF values underflows double precision).

**Log-gamma (for Beta PDF display).** Uses the Lanczos 6-term approximation from Numerical Recipes, providing accurate evaluation of lgamma for all positive real arguments arising in the Beta distribution parameter range of the calculator.

---

## 6. Reliability of the Normal Approximation

The score test and Beta-Binomial score test rely on the central limit theorem. The approximation is generally reliable when:

- N·p₀ ≥ 10 and N·(1 − p₀) ≥ 10 (standard rule of thumb for binomial CLT).
- For the Beta-Binomial, the effective sample size N/VIF replaces N; the condition becomes N·p₀/VIF ≥ 10 and N·q₀/VIF ≥ 10.

For the primary use case (p₀ near 0.90 and moderate effect sizes), the test tends to be one-sided near the upper tail, and p₀ close to 1 makes the approximation less accurate on the lower tail side. The exact binomial test (Tab 1) should be preferred when N < 100 or when p₀ < 0.10 or p₀ > 0.90 and high precision is required. The score test with continuity correction performs better than the uncorrected score test in these boundary regions.

---

## 7. Reference Validation Values

The following values are produced by the calculator at default settings: p₀ = 0.90, p₁ = 0.95, α = 0.05 (one-sided), target power = 0.80. They serve as a reproducibility check.

| Method | N (total) | n (questions) | k | VIF | Critical count c | Achieved α | Power |
|---|---|---|---|---|---|---|---|
| Exact binomial | 179 | — | — | — | 168 | 0.0473 | 0.8011 |
| Score z-test | 184 | — | — | — | 173 | (nominal) | 0.8017 |
| BB score test | 222 | 74 | 3 | 1.20 | — | — | 0.8041 |

Additional spot checks (ρ = 0.10, k = 3):
- naiveTrueAlpha = 0.0666 (versus nominal 0.05)
- Edge case k = 1: VIF = 1 for any ρ; BB result matches standard binomial

---

## 8. Reproducibility Checklist

A methods section reporting results from this calculator should specify all of the following parameters. Omitting any one makes exact reproduction impossible.

1. **Test method**: exact binomial / score z-test / score z-test with CC / Beta-Binomial score test
2. **Null hypothesis value** p₀
3. **Alternative hypothesis value** p₁ (or μ₁ for Beta-Binomial)
4. **One-sided type I error rate** α
5. **Target power** (1 − β)
6. **For Beta-Binomial only**: intra-class correlation ρ and responses per question k (from which VIF = 1 + (k−1)ρ is derived)
7. **Reported outputs**: N (and n, k for BB), critical count c, achieved α, and achieved power
8. **Software**: "browser-based LLM Evaluation Power Calculator (open-source, available at [URL])"

Example reporting sentence: "Sample size was determined using an exact binomial one-sided test (p₀ = 0.90, p₁ = 0.95, α = 0.05, power = 0.80), yielding N = 179 responses with critical count c = 168, achieved α = 0.047, and achieved power = 0.801."

---

## 9. Method Choice: Wilson Score vs. Clopper-Pearson

The two methods answer the same question but make different trade-offs.

**Clopper-Pearson (exact binomial).** The type I error is guaranteed to be at or below the nominal α for every N and p₀. This is the sense in which the method is "exact." The cost is conservatism: because the binomial is discrete, the achieved α is often strictly below the nominal (e.g., 4.86% instead of 5%), which reduces power slightly. For regulatory or clinical submissions where a strict upper bound on α is required, CP is the standard choice.

**Wilson score.** The score test evaluates the variance at the null (p₀q₀/N) rather than at the observed proportion (p̂q̂/N). This makes it better calibrated than the Wald test across the full range of p, including near 0 and 1 — a regime relevant to LLM evaluation. Three key evaluations support Wilson as the default for general planning:

- **Wilson (1927) [1]** introduced the score interval and showed it has better finite-sample behaviour than the normal approximation with observed-proportion variance (Wald).
- **Agresti & Coull (1998) [2]** demonstrated that the Wilson interval (and its "add-2" variant) outperforms the "exact" CP interval in terms of average coverage, arguing that CP's over-coverage is itself a defect rather than a virtue when the goal is to attain the nominal level. Their title — "Approximate is better than 'exact'" — summarises the finding.
- **Brown, Cai & DasGupta (2001) [3]** provided a comprehensive coverage-probability analysis showing that the Wald interval fails badly near p = 0 and p = 1 (the regime of this calculator), while Wilson and Jeffreys intervals maintain near-nominal coverage throughout. The CP interval works well but systematically over-covers.

**Practical summary.** For p₀ near 0.90 (the default for LLM evaluation) and N in the range 100–400, the power difference between Wilson and CP is small (typically 1–3 pp). The CP method guarantees α control; the Wilson method has better average coverage and a closed-form sample size formula. Either can be reported; state which you used and cite accordingly.

**What not to use: the arcsine (Cohen's h) approximation.** The `pwr.p.test()` function in R uses the arcsine-transformation effect size h = 2 arcsin(√p₁) − 2 arcsin(√p₀), which is a variance-stabilising transform designed for proportions near 0.5. For p₀ = 0.90, p₁ = 0.95, it gives N = 167, but at N = 167 the Wilson test achieves only 75.8% power (below the 80% target). The correct sample sizes are N = 184 (Wilson) and N = 179 (CP exact). Use the arcsine approximation only when p is near 0.5.

---

## 10. References

[1] Wilson, E. B. (1927). Probable inference, the law of succession, and statistical inference. *Journal of the American Statistical Association*, 22(158), 209–212. https://doi.org/10.1080/01621459.1927.10502953

[2] Agresti, A., & Coull, B. A. (1998). Approximate is better than "exact" for interval estimation of binomial proportions. *The American Statistician*, 52(2), 119–126. https://doi.org/10.1080/00031305.1998.10480550

[3] Brown, L. D., Cai, T. T., & DasGupta, A. (2001). Interval estimation for a binomial proportion. *Statistical Science*, 16(2), 101–133. https://doi.org/10.1214/ss/1009213286

[4] Fleiss, J. L., Levin, B., & Paik, M. C. (2003). *Statistical Methods for Rates and Proportions* (3rd ed.). John Wiley & Sons. https://doi.org/10.1002/0471445428

[5] Lachin, J. M. (1981). Introduction to sample size determination and power analysis for clinical trials. *Controlled Clinical Trials*, 2(2), 93–113. https://doi.org/10.1016/0197-2456(81)90001-5

---

*Document version: September 2026. All formulas have been verified against the calculator's JavaScript implementation.*
