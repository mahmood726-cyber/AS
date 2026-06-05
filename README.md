# AS-Logic — Asymptomatic Aortic Stenosis Evidence Synthesis

A single-file, **fully offline** dashboard that pools hazard ratios for early
intervention (TAVR / SAVR) versus surveillance in asymptomatic severe aortic
stenosis, using a DerSimonian–Laird random-effects model with a
Hartung–Knapp–Sidik–Jonkman small-sample adjustment, a forest plot, a
subgroup-interaction test, an NNT/waffle translation panel, and a GRADE summary.

**Live app:** open `index.html` (or the GitHub Pages link). No build step, no
network, no external CDN.

## Layout

```
index.html   single-file UI (loads engine.js)
engine.js    pure statistical core — runs in Node and the browser
tests.js     Node test harness, 42 assertions
LICENSE      Apache-2.0
```

## Statistical core (`engine.js`)

| Function | What it does |
|---|---|
| `calculateMetaAnalysis(trials)` | DerSimonian–Laird random-effects pooling of logHR: derives per-study SE from the reported 95% CI, computes Q, I² = max(0,(Q−df)/Q), τ² = (Q−df)/C, an HKSJ SE floor (`finalSE = max(RE-SE, HKSJ-SE)`), a `t_{k−1}` confidence interval and a `t_{k−2}` prediction interval |
| `normalCdf(x)` | standard-normal **upper tail** Q(x)=1−Φ(x), valid for x ≥ 0 |
| `chiSqPval(chi, df)` | two-sided χ²(1) tail probability (subgroup-interaction p-value) |
| `getTScore(df)` | tabulated two-sided 0.975 t critical values |

## Fixes applied during revival (2026-06-05)

- **Offline:** removed the Google Fonts `<link>` (the only external resource);
  the page now loads nothing over the network (system fonts fall back).
- **Single source of truth:** extracted the inline statistical core into a pure
  `engine.js` and deleted the inline duplicates; the page now loads `engine.js`.
- **Correctness bug fixed (`chiSqPval`, df=1):** the original returned
  `2·Φ(√chi)` instead of `2·(1−Φ(√chi))`. Because `normalCdf` here is the
  *upper* tail Q(x)=1−Φ(x), the original expression evaluated to `2·Φ(√chi)`,
  which produces impossible p-values **> 1** (e.g. χ²=3.841 → 1.95 instead of
  0.05) and corrupted every displayed subgroup-interaction p-value. It now
  returns `2·normalCdf(√chi)` = `2·(1−Φ(√chi))`, verified against the χ²(1)
  critical values (3.841→0.05, 6.635→0.01), and locked by tests.
- Added `tests.js` (42 assertions, all passing), `.nojekyll`, `.gitignore`,
  this README, and `E156-PROTOCOL.md`; renamed `AS.html` → `index.html`.

The `normalCdf` function was checked and is **correct** (Q(0)=0.5,
Q(1.96)=0.025) — only its *caller* `chiSqPval` had the sign error. The pooling
math (τ², I², HKSJ floor, weights) was verified against hand computation and
left unchanged.

## Tests

```
node tests.js
# 42 passed, 0 failed
```

Checks include normal-tail reference points (Q(0)=0.5, Q(1.96)=0.025),
χ²(1) p-values against critical values, the k=1 single-trial passthrough, a
two-identical-trial case (τ²=0, I²=0, 50/50 weights), and a fully hand-worked
heterogeneous two-study case (se₁≈0.1804, se₂≈0.3296, Q≈3.014, I²≈66.8%,
τ²≈0.1422, pooled HR≈0.654, weights ≈58.9/41.1).

## Caveats

The displayed model derives each study's variance from its **published CI**
(z=3.92 for 95%), not from raw event counts. DerSimonian–Laird under-estimates
τ² for small *k* (REML/Paule–Mandel are preferred for k<10); the HKSJ
adjustment here is implemented as a conservative SE floor with `t_{k−1}` CIs.
With only a handful of asymptomatic-AS trials the prediction intervals are very
wide — treat pooled HRs as hypothesis-generating, not as a clinical decision
rule. Apache-2.0 licensed.
