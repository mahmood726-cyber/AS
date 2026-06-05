# E156-PROTOCOL — AS-Logic (Asymptomatic Aortic Stenosis Synthesis)

- **Project:** AS (GitHub repo `AS`, user `mahmood726-cyber`)
- **Revived:** 2026-06-05 (from a single-file `AS.html` dump)
- **Type:** single-file offline browser tool + Node-testable engine
- **Dashboard:** GitHub Pages (`index.html`)

## What changed in the revival

- Made **fully offline**: removed the Google Fonts CDN `<link>` (the only
  external resource); the app now loads nothing over the network.
- Extracted the statistical core into a pure `engine.js` (single source of
  truth; the inline duplicates were removed and the page now loads `engine.js`).
- **Fixed a correctness bug** in `chiSqPval` (df=1): it returned `2·Φ(√chi)`
  instead of `2·(1−Φ(√chi))`, yielding impossible interaction p-values > 1
  (χ²=3.841 gave 1.95, not 0.05). Now correct and test-locked.
- Added `tests.js` (42 assertions, all passing), `.nojekyll`, `.gitignore`,
  README; renamed `AS.html` → `index.html`.

## Body (E156 draft — CURRENT BODY)

Does early valve intervention reduce death or hospitalisation in asymptomatic
severe aortic stenosis, and how stable is that pooled signal across surgical and
transcatheter approaches? This dashboard collates published hazard ratios from
the RECOVERY, AVATAR, EVOLVED and EARLY-TAVR trials. Each study's variance is
derived from its reported 95% confidence interval on the log scale. It then
pools log hazard ratios with a DerSimonian–Laird random-effects model, adding a
Hartung–Knapp–Sidik–Jonkman standard-error floor, t-based confidence and
prediction intervals, and a subgroup-interaction test. The bundled set yields a
pooled hazard ratio near 0.65 favouring intervention, but with high
heterogeneity and prediction intervals that approach unity once between-trial
variance is honoured. A revival
audit found and fixed a chi-square bug that had produced impossible
subgroup-interaction p-values above one, and locked the core behind a
42-assertion hand-checked test suite. The honest read is a promising but
heterogeneity-limited benefit, presented as a transparent exploratory synthesis
rather than a clinical decision rule.

SUBMITTED: [ ]
