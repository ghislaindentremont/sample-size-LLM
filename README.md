# How Many Queries Are Enough? — Sample Size Calculator

Companion repository to the *JMIR AI* editorial **How Many Queries Are Enough? Sample Size Calculations in Early-Stage LLM Validation**. It contains a single-page, dependency-free browser calculator ([`proportion-calculator/index.html`](proportion-calculator/index.html)), the R scripts behind the editorial's figures, and the Monte Carlo simulation reported in the Appendix.

## What the calculator does

Domain experts rate each of N LLM responses to unique user queries (query-response pairs) as acceptable or unacceptable, and the study tests whether the true proportion of acceptable LLM responses exceeds a pre-defined absolute performance threshold (default p₀ = 0.90; one-sided). Enter your parameters to get:

- the **sample size** (number of query-response pairs) needed to reach a target power, or
- the **power** achieved at a fixed N.

Two methods are computed side by side:

| Method | Approach | Role |
|---|---|---|
| **Wilson score** | Normal-approximation score test formula; closed-form | Recommended (better coverage than Wald and Clopper-Pearson) |
| **Clopper-Pearson** | Exact binomial; no approximation | Shown for comparison; conservative |

Copy-ready R and Python code is generated for every result.

A second tab, **Repetitions (Beta Binomial)**, covers studies that generate k LLM responses per query. It applies the variance inflation factor VIF = 1 + (k − 1)ρ (ρ = intra-query correlation) to the Wilson score formulas, shows what the uncorrected test would report versus what it actually delivers, and gives the cost-optimal number of repetitions per query.

## Running it

Open `proportion-calculator/index.html` in a browser. Everything runs client-side; no build step or server needed.

To host on GitHub Pages: Settings → Pages → deploy from the `main` branch, folder `/proportion-calculator`.

## Reference values

For the defaults (p₀ = 0.90, p₁ = 0.95, α = 0.05, 80% power):

| Method | N | Reject H₀ if |
|---|---|---|
| Wilson score | 184 | X ≥ 173 |
| Clopper-Pearson exact | 179 | X ≥ 168 |

## Repository layout

```
proportion-calculator/   active website
analysis/                R validation script and manuscript figures
.archive/                archived materials (not part of the active site)
```
