# LLM Evaluation Power Calculator

A single-page, dependency-free website for planning evaluation studies of LLM
responses. Researchers judge N responses as acceptable or not, and want to show
that the true acceptable rate exceeds a performance threshold (default 0.90)
using a one-sided, one-sample proportion test. Given an expected proportion
(default 0.95) and a significance level, the page computes either

- the **sample size** needed to reach a target power, or
- the **power** achieved with a fixed N,

using either an exact binomial test (recommended) or the normal approximation.
It also states the statistical assumptions under which the test is valid.

## Running it

Open `index.html` in a browser. Everything is computed client-side; there is no
build step and no server.

To host on GitHub Pages: Settings → Pages → deploy from the `main` branch root.

## Method

- **Hypotheses:** H₀: p ≤ p₀ vs H₁: p > p₀, where p is the true probability a
  response is judged acceptable.
- **Exact binomial:** the critical count c is the smallest integer with
  P(X ≥ c | N, p₀) ≤ α; power is P(X ≥ c | N, p₁). Sample size is the smallest
  N whose power reaches the target.
- **Normal approximation:**
  N = ⌈((z_α√(p₀q₀) + z_β√(p₁q₁)) / (p₁ − p₀))²⌉ and
  power = Φ(((p₁ − p₀)√N − z_α√(p₀q₀)) / √(p₁q₁)).

Reference values for the defaults (p₀ = 0.90, p₁ = 0.95, α = 0.05, power = 0.80):
exact N = 179 (reject if X ≥ 168), normal N = 184.
