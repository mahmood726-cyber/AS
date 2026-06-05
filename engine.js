/*
 * AS-Logic engine — pure meta-analysis core for the asymptomatic aortic
 * stenosis (early-intervention vs surveillance) evidence-synthesis dashboard.
 *
 * Extracted verbatim from the dashboard's inline <script> so the statistical
 * core is a single source of truth, importable under Node for testing.
 * Browser: functions are globals (plain declarations). Node: module.exports.
 *
 * Method: DerSimonian-Laird tau^2, I^2/Q on the log-HR scale, with a
 * Hartung-Knapp-Sidik-Jonkman small-sample SE adjustment applied as a variance
 * floor (finalSE = max(RE-SE, HKSJ-SE)) and t_{k-1} critical values for the CI.
 * Faithful to the shipped app — no methodology changes.
 *
 * One correctness fix applied during the 2026-06 revival (see chiSqPval): the
 * original returned 2*Phi(sqrt(chi)) for df=1, which yields impossible p-values
 * > 1 (e.g. chi=3.841 -> 1.95 instead of 0.05) and corrupted every displayed
 * subgroup-interaction p-value. It now returns 2*(1-Phi(sqrt(chi))) correctly.
 */

// T-Distribution lookup (Two-tailed, alpha=0.05). Keys are degrees of freedom.
const tTable = {
    1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
    6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
    100: 1.984 // fallback for large N
};

function getTScore(df) {
    if (df <= 0) return 12.7; // Ultra conservative for 1 study
    if (df > 10) return 1.96;
    return tTable[df] || 2.0;
}

// Standard normal UPPER-TAIL probability Q(x) = 1 - Phi(x), valid for x >= 0
// (Zelen & Severo 26.2.17 polynomial). Q(0)=0.5, Q(1.96)=0.025, Q(1.645)=0.05.
function normalCdf(x) {
    var t = 1 / (1 + .2316419 * x);
    var d = .3989423 * Math.exp(-x * x / 2);
    var p = d * t * (.3193815 + t * (-.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))));
    return p;
}

// Two-sided chi-square(df) tail probability (only df=1 is used by the app).
// FIX (2026-06 revival): normalCdf is the UPPER tail Q(z)=1-Phi(z), so the
// two-sided chi^2(1) p-value is 2*Q(sqrt(chi)) = 2*normalCdf(sqrt(chi)). The
// original wrote 2*(1 - normalCdf(...)) = 2*Phi(...), producing p-values > 1.
function chiSqPval(chi, df) {
    if (df <= 0) return 0;
    // Approximation for df=1 (most common here)
    if (df === 1) {
        // ChiSq(1) is the square of a standard Normal.
        return 2 * normalCdf(Math.sqrt(chi));
    }
    return 1; // Placeholder for high DF, not needed for this app
}

function calculateMetaAnalysis(subsetTrials) {
    if (!subsetTrials || subsetTrials.length === 0) return null;
    const k = subsetTrials.length;

    // 1. Calculate Log HR and Variance per study
    let sumWi = 0, sumWiYi = 0, sumWiYi2 = 0, sumWi2 = 0;

    subsetTrials.forEach(t => {
        // Derive variance from the reported CI (log scale, z=3.92 for 95%).
        t.logHR = Math.log(t.hr);
        t.se = (Math.log(t.ci_upper) - Math.log(t.ci_lower)) / 3.92;
        t.variance = t.se * t.se;

        // Fixed Effect Weight (Inverse Variance)
        t.fe_weight = 1 / t.variance;

        sumWi += t.fe_weight;
        sumWiYi += t.fe_weight * t.logHR;
        sumWiYi2 += t.fe_weight * Math.pow(t.logHR, 2);
        sumWi2 += (t.fe_weight * t.fe_weight);
    });

    // 2. Heterogeneity (Q and I^2)
    const Q = subsetTrials.reduce((acc, t) => acc + t.fe_weight * Math.pow(t.logHR - (sumWiYi / sumWi), 2), 0);
    const df = k - 1;
    const I2 = k > 1 ? Math.max(0, (Q - df) / Q) * 100 : 0;

    // 3. Tau^2 (DerSimonian-Laird)
    let tau2 = 0;
    if (k > 1 && Q > df) {
        const c = sumWi - (sumWi2 / sumWi);
        tau2 = (Q - df) / c;
    }

    // 4. Random Effects Weights
    let sumWiStar = 0, sumWiStarYi = 0;
    subsetTrials.forEach(t => {
        t.re_weight = 1 / (t.variance + tau2);
        sumWiStar += t.re_weight;
        sumWiStarYi += t.re_weight * t.logHR;
    });

    // Normalize weights
    subsetTrials.forEach(t => t.final_weight_pct = (t.re_weight / sumWiStar) * 100);

    // 5. Pooled Point Estimate
    const pooledLogHR = sumWiStarYi / sumWiStar;

    // 6. HKSJ Standard Error Adjustment (Small Sample Correction)
    let finalSE = Math.sqrt(1 / sumWiStar); // Default standard RE SE

    if (k >= 2) {
        // Weighted variance of estimates
        const weightedVar = subsetTrials.reduce((acc, t) => {
            return acc + t.re_weight * Math.pow(t.logHR - pooledLogHR, 2);
        }, 0);

        const hksjFactor = weightedVar / ((k - 1) * sumWiStar);
        const seHKSJ = Math.sqrt(hksjFactor);

        // Use the larger of the two to be conservative (Knapp-Hartung logic)
        finalSE = Math.max(finalSE, seHKSJ);
    }

    // 7. Confidence Intervals using T-distribution (Dynamic DF)
    // Use k-1 for the CI t-score
    const tScoreCI = getTScore(k - 1);
    const pooledLower = Math.exp(pooledLogHR - tScoreCI * finalSE);
    const pooledUpper = Math.exp(pooledLogHR + tScoreCI * finalSE);
    const pooledHR = Math.exp(pooledLogHR);

    // 8. Prediction Interval (Dynamic)
    // PI incorporates Tau2 and uses t-distribution
    const sePred = Math.sqrt((finalSE * finalSE) + tau2);
    // Use k-2 for PI t-score (standard convention for PI), but if k=2 floor at 1.
    const dfPI = Math.max(1, k - 2);
    const tScorePI = getTScore(dfPI);

    const piLower = Math.exp(pooledLogHR - tScorePI * sePred);
    const piUpper = Math.exp(pooledLogHR + tScorePI * sePred);

    // 9. Significance (checked against the CI, as in the shipped app)
    const isSignificant = (pooledLower < 1.0 && pooledUpper < 1.0) || (pooledLower > 1.0 && pooledUpper > 1.0);

    return {
        pooledHR, pooledLower, pooledUpper, pooledLogHR,
        piLower, piUpper,
        I2, tau2, Q, df,
        subsetTrials,
        totalN: subsetTrials.reduce((sum, t) => sum + t.n, 0),
        isSignificant
    };
}

if (typeof module !== 'undefined' && module.exports) {
    module.exports = { calculateMetaAnalysis, getTScore, normalCdf, chiSqPval, tTable };
}
