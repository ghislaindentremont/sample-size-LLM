# LLM Evaluation Sample Size Calculator

A single-page, dependency-free browser tool for planning LLM response evaluation studies. The calculator lives at [`proportion-calculator/index.html`](proportion-calculator/index.html).

## What it does

Researchers judge N LLM responses as acceptable or not, and want to show that the true acceptable rate exceeds a performance threshold (default p₀ = 0.90) at a given significance level. Enter your parameters to get:

- the **sample size** needed to reach a target power, or
- the **power** achieved at a fixed N.

Two methods are computed side by side:

| Method | Approach | When to use |
|---|---|---|
| **Wilson score** | Normal (CLT) approximation; closed-form | Default — reliable for N ≥ 100 at p₀ = 0.90 |
| **Clopper-Pearson exact** | Exact binomial; no approximation | Strict type I error guarantee, or N < 100 |

Copy-ready R and Python code is generated for every result.

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
