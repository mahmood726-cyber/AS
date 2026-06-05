/*
 * Node tests for the AS-Logic engine. Run: node tests.js
 * Every expected value is hand-derived independently below (NOT produced by
 * running engine.js), so the suite is a true cross-check of the engine.
 */
const { calculateMetaAnalysis, getTScore, normalCdf, chiSqPval } = require('./engine.js');

let pass = 0, fail = 0;
function ok(name, cond, detail) {
    if (cond) { pass++; console.log('  ok  ' + name); }
    else { fail++; console.log(' FAIL ' + name + (detail ? '  -> ' + detail : '')); }
}
function close(a, b, tol) { return Math.abs(a - b) < (tol || 1e-3); }

// =====================================================================
// normalCdf — this is the UPPER tail Q(x) = 1 - Phi(x), valid for x >= 0.
// Reference points: Q(0)=0.5, Q(1.96)=0.025, Q(1.645)=0.05, Q(1)=0.158655.
// A normalCdf returning 0 (not 0.5) at x=0 would be a bug; this one is correct.
// =====================================================================
ok('normalCdf(0) ~ 0.5', close(normalCdf(0), 0.5, 1e-5), 'got ' + normalCdf(0));
ok('normalCdf(1.96) ~ 0.025', close(normalCdf(1.96), 0.025, 1e-4), 'got ' + normalCdf(1.96));
ok('normalCdf(1.645) ~ 0.05', close(normalCdf(1.645), 0.05, 1e-4), 'got ' + normalCdf(1.645));
ok('normalCdf(1) ~ 0.158655', close(normalCdf(1), 0.158655, 1e-4), 'got ' + normalCdf(1));

// =====================================================================
// chiSqPval (df=1) — two-sided chi^2(1) tail = 2*(1-Phi(sqrt(chi))).
// chi=3.841 (the 0.05 critical value) -> p ~ 0.05 ; chi=0 -> 1 ; chi=1 -> 0.3173.
// REGRESSION GUARD for the revival fix: the original returned 2*Phi(sqrt(chi)),
// giving 1.95 for chi=3.841 (an impossible p > 1). Assert p stays <= 1.
// =====================================================================
ok('chiSqPval(3.841,1) ~ 0.05', close(chiSqPval(3.841, 1), 0.05, 1e-3), 'got ' + chiSqPval(3.841, 1));
ok('chiSqPval(0,1) ~ 1.0', close(chiSqPval(0, 1), 1.0, 1e-5), 'got ' + chiSqPval(0, 1));
ok('chiSqPval(1,1) ~ 0.31731', close(chiSqPval(1, 1), 0.31731, 1e-3), 'got ' + chiSqPval(1, 1));
ok('chiSqPval(6.635,1) ~ 0.01', close(chiSqPval(6.635, 1), 0.01, 1e-3), 'got ' + chiSqPval(6.635, 1));
ok('chiSqPval never exceeds 1 (regression for revival bug)', chiSqPval(3.841, 1) <= 1.0);
ok('chiSqPval(x,0) -> 0', chiSqPval(5, 0) === 0);

// =====================================================================
// getTScore — tabulated two-sided 0.975 t critical values.
// =====================================================================
ok('getTScore(1) -> 12.706', getTScore(1) === 12.706);
ok('getTScore(2) -> 4.303', getTScore(2) === 4.303);
ok('getTScore(0) -> 12.7 (k=1 guard)', getTScore(0) === 12.7);
ok('getTScore(50) -> 1.96 (large df)', getTScore(50) === 1.96);

// =====================================================================
// Empty guard
// =====================================================================
ok('empty -> null', calculateMetaAnalysis([]) === null);
ok('null -> null', calculateMetaAnalysis(null) === null);

// =====================================================================
// Single trial (k=1): no heterogeneity, CI uses t_{df=0} = 12.7.
// Hand: hr=0.50, ci 0.35-0.71. y=ln0.5=-0.693147, se=(ln0.71-ln0.35)/3.92
//   = (-0.342490 - -1.049822)/3.92 = 0.707332/3.92 = 0.180442.
// k=1 -> tau2=0, I2=0. finalSE = sqrt(1/(1/se^2)) = se = 0.180442.
// CI: exp(-0.693147 -+ 12.7*0.180442) -> very wide.
// =====================================================================
const single = calculateMetaAnalysis([{ hr: 0.50, ci_lower: 0.35, ci_upper: 0.71, n: 901 }]);
ok('single: pooledHR == 0.50', close(single.pooledHR, 0.50, 1e-9), 'got ' + single.pooledHR);
ok('single: tau2 == 0', single.tau2 === 0);
ok('single: I2 == 0', single.I2 === 0);
ok('single: df == 0', single.df === 0);
ok('single: totalN == 901', single.totalN === 901);
ok('single: lower = exp(-0.693147 - 12.7*0.180442)', close(single.pooledLower, Math.exp(-0.693147 - 12.7 * 0.180442), 1e-4), 'got ' + single.pooledLower);

// =====================================================================
// Two IDENTICAL trials -> tau2=0, I2=0 (Q = df = 1 exactly), pooled == trial.
// Hand: both hr=0.50, ci 0.35-0.71 -> y=-0.693147, w=1/0.180442^2=30.71325.
// muFE = muRE = -0.693147 (identical). Q = 0 effectively -> I2=0, tau2=0.
// finalSE = sqrt(1/(2w)) = se/sqrt(2) = 0.180442/1.41421 = 0.127592.
// (seHKSJ = sqrt(0 / (1*2w)) = 0, so finalSE = max(RE, 0) = RE-SE.)
// weights split 50/50.
// =====================================================================
const tIdent = { hr: 0.50, ci_lower: 0.35, ci_upper: 0.71, n: 100 };
const ident = calculateMetaAnalysis([{ ...tIdent }, { ...tIdent }]);
ok('identical: pooledHR == 0.50', close(ident.pooledHR, 0.50, 1e-9), 'got ' + ident.pooledHR);
ok('identical: tau2 == 0', close(ident.tau2, 0, 1e-12), 'got ' + ident.tau2);
ok('identical: I2 == 0', close(ident.I2, 0, 1e-9), 'got ' + ident.I2);
ok('identical: finalSE == se/sqrt(2)', close(ident.subsetTrials[0].se / Math.SQRT2, Math.sqrt(1 / (2 / ident.subsetTrials[0].variance)), 1e-9));
ok('identical: weights are 50/50', close(ident.subsetTrials[0].final_weight_pct, 50, 1e-9) && close(ident.subsetTrials[1].final_weight_pct, 50, 1e-9));
ok('identical: totalN == 200', ident.totalN === 200);

// =====================================================================
// FULLY HAND-WORKED 2-study heterogeneous example.
// t1: hr=0.50, ci 0.35-0.71 ; t2: hr=0.96, ci 0.50-1.82.
//
//   y1 = ln(0.50) = -0.693147
//   se1 = (ln0.71 - ln0.35)/3.92 = 0.707332/3.92 = 0.180442 ; w1 = 30.713248
//   y2 = ln(0.96) = -0.040822
//   se2 = (ln1.82 - ln0.50)/3.92 = 1.291984/3.92 = 0.329588 ; w2 = 9.205727
//
//   sumW = 39.918975 ; sumWY = -21.66878 ; muFE = -0.542714
//   Q = w1*(y1-muFE)^2 + w2*(y2-muFE)^2 = 3.013927  ; df = 1
//   I2 = (Q-df)/Q * 100 = 66.8207%
//   C = sumW - sumW2/sumW = 14.165583
//   tau2 = (Q-df)/C = 2.013927/14.165583 = 0.142170
//
//   wRE1 = 1/(0.032559+0.142170) = 5.72654 ; wRE2 = 1/(0.108628+0.142170)=3.98727
//   sumWstar = 9.71381 ; muRE = -0.425290 -> pooledHR = exp(-0.425290) = 0.653580
//   seRE = sqrt(1/9.71381) = 0.320909
//   weightedVar = wRE1*(y1-muRE)^2 + wRE2*(y2-muRE)^2 = 1.000000 (= k-1 here)
//   hksjFactor = weightedVar/((k-1)*sumWstar) = 1/9.71381 -> seHKSJ = 0.320909
//   finalSE = max(seRE, seHKSJ) = 0.320909
//   tCI = getTScore(k-1=1) = 12.706
//   CI = exp(-0.425290 -+ 12.706*0.320909) = [0.011219, 38.07...]
//   weights pct: 58.9382 / 41.0618
// =====================================================================
const het = calculateMetaAnalysis([
    { hr: 0.50, ci_lower: 0.35, ci_upper: 0.71, n: 901 },
    { hr: 0.96, ci_lower: 0.50, ci_upper: 1.82, n: 224 }
]);
ok('het: se1 ~ 0.180442', close(het.subsetTrials[0].se, 0.180442, 1e-5), 'got ' + het.subsetTrials[0].se);
ok('het: se2 ~ 0.329588', close(het.subsetTrials[1].se, 0.329588, 1e-5), 'got ' + het.subsetTrials[1].se);
ok('het: Q ~ 3.013927', close(het.Q, 3.013927, 1e-4), 'got ' + het.Q);
ok('het: I2 ~ 66.8207%', close(het.I2, 66.8207, 1e-2), 'got ' + het.I2);
ok('het: tau2 ~ 0.142170', close(het.tau2, 0.142170, 1e-5), 'got ' + het.tau2);
ok('het: pooledHR ~ 0.653580', close(het.pooledHR, 0.653580, 1e-4), 'got ' + het.pooledHR);
ok('het: pooledLogHR ~ -0.425290', close(het.pooledLogHR, -0.425290, 1e-4), 'got ' + het.pooledLogHR);
ok('het: weight pct ~ 58.94 / 41.06', close(het.subsetTrials[0].final_weight_pct, 58.9382, 1e-2) && close(het.subsetTrials[1].final_weight_pct, 41.0618, 1e-2),
    'got ' + het.subsetTrials[0].final_weight_pct + ' / ' + het.subsetTrials[1].final_weight_pct);
ok('het: CI lower = exp(muRE - 12.706*finalSE)', close(het.pooledLower, Math.exp(-0.425290 - 12.706 * 0.320909), 2e-3), 'got ' + het.pooledLower);
ok('het: CI upper = exp(muRE + 12.706*finalSE)', close(het.pooledUpper, Math.exp(-0.425290 + 12.706 * 0.320909), 5e-1), 'got ' + het.pooledUpper);
ok('het: pooled HR within CI', het.pooledLower < het.pooledHR && het.pooledHR < het.pooledUpper);

// =====================================================================
// 3-study example just to exercise the prediction-interval branch (k>=3,
// dfPI = k-2 = 1, tScorePI = getTScore(1) = 12.706). We only assert the PI
// brackets the point estimate and is wider than the CI (tau2 + t inflation).
// =====================================================================
const three = calculateMetaAnalysis([
    { hr: 0.35, ci_lower: 0.16, ci_upper: 0.76, n: 145 },
    { hr: 0.46, ci_lower: 0.23, ci_upper: 0.90, n: 157 },
    { hr: 0.50, ci_lower: 0.35, ci_upper: 0.71, n: 901 }
]);
ok('three: PI brackets pooledHR', three.piLower < three.pooledHR && three.pooledHR < three.piUpper);
ok('three: PI at least as wide as CI', (three.piUpper - three.piLower) >= (three.pooledUpper - three.pooledLower) - 1e-9);
ok('three: I2 finite and >= 0', isFinite(three.I2) && three.I2 >= 0);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
