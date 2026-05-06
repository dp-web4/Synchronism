# Proposal: Session 107 fσ₈ Predictions Disfavored by DESI DR1 (≥2σ)

**Filed**: 2026-05-05 (back-annotation from synchronism-site explorer track)
**Source finding**: `synchronism-site/explorer/findings/desi-dr1-vs-session107-fsigma8.md`
**Status when written**: Session 107 was the framework's leading current-data
cosmological discriminator, advertised on the site (TEST-04a) as the most
decisive comparison waiting to be run.

---

## TL;DR

Session 107 (Dec 2025) predicted fσ₈ ~10–12% suppression below ΛCDM at z = 0.5–0.7,
framed as "the most important observational test for Synchronism in the coming
years" and forecast at 3.1–3.2σ per LRG bin with DESI Y1 data.

DESI DR1 published full-shape results (DESI Collaboration, A.G. Adame et al.,
arXiv:2411.12021 v5, Nov 2024). Results:

| Bin | z_eff | Sync prediction (fσ₈) | DESI measurement | Verdict |
|-----|-------|-----------------------|-------------------|---------|
| LRG1 | 0.51 | 0.418 | ~0.55 ± 0.06 | 2σ above Sync |
| LRG2 | 0.71 | 0.414 | ~0.50 ± 0.05 | 1.4σ above Sync |
| Combined σ₈(z=0) | 0 | 0.76 | 0.841 ± 0.034 | 2.4σ above Sync |

By Session 107's own falsification criterion (fσ₈(z=0.5) > 0.45 → ΛCDM favored),
ΛCDM is favored at every LRG bin and at the combined fit.

This is a **refuted-not-untested** prediction. The honest framing is that the
test fired in the wrong direction.

---

## What Synchronism actually predicted

From `Research/Session107_DESI_Forecasts.md`:

- Mechanism: G_local/G_global = C_cosmic/C_galactic < 1 during structure formation
  → suppresses growth rate f(z); combined with σ₈(z=0) = 0.76 (vs ΛCDM 0.81)
  gives lower fσ₈.
- Prediction pattern: largest suppression at low z (cumulative effect), shrinking
  at high z (theories converge near matter domination).
- Numeric forecast: 3σ per bin at LRG1/LRG2 with DESI Y1; 6.6σ combined at DESI
  Final.
- Falsification ladder: fσ₈(z=0.5) > 0.45 → ΛCDM favored.

---

## What DESI DR1 actually measured

ShapeFit + BAO model-agnostic fσ₈/(fσ₈)^fid (DESI 2024 V Table 9):

| Bin | DESI ratio | Sync ratio | Tension (σ above Sync) |
|-----|-----------|-----------|-----------------------|
| LRG1 | 1.16 ± 0.13 | 0.882 | 2.14σ |
| LRG2 | 1.04 +0.11/−0.092 | 0.898 | 1.42σ |
| LRG3 | 0.997 +0.10/−0.084 | 0.916 | 0.88σ |
| ELG2 | 0.945 +0.097/−0.077 | 0.932 | 0.15σ (consistent) |
| QSO | 1.16 ± 0.12 | 0.947 | 1.78σ |
| BGS | 0.84 ± 0.19 | 0.867 | 0.14σ (consistent) |

Per-bin ΛCDM σ₈(z=0) inferred (DESI 2024 V Table 10) vs Sync's 0.76:

| Bin | DESI σ₈ | Tension vs Sync 0.76 |
|-----|---------|----------------------|
| BGS | 0.662 ± 0.13 | −0.75σ (below; consistent) |
| LRG1 | 0.835 ± 0.087 | +0.86σ |
| LRG2 | 0.880 +0.072/−0.082 | +1.50σ |
| LRG3 | 0.815 +0.068/−0.076 | +0.76σ |
| ELG2 | 0.755 +0.054/−0.064 | −0.08σ (bullseye) |
| QSO | 0.950 +0.066/−0.077 | +2.59σ |
| **Combined** | **0.841 ± 0.034** | **+2.38σ** |

DESI 2024 V abstract: "DESI DR1 galaxy clustering results are in agreement with
the ΛCDM model based on general relativity with parameters consistent with those
from Planck."

---

## The inverted redshift pattern

Session 107's mechanism predicts low-z **suppressed**, high-z **converging to ΛCDM**.

Observed pattern:

```
                σ_8 inferred (FM+BAO)
LRG1 (z=0.51)   0.84     HIGH vs Sync's 0.76
LRG2 (z=0.71)   0.88     HIGH
LRG3 (z=0.93)   0.82     HIGH
ELG2 (z=1.32)   0.76     ON Sync's prediction (bullseye)
QSO  (z=1.49)   0.95     VERY HIGH (also high vs Planck)
```

The bullseye at ELG2 is *not* a confirmation: ELG2 is the high-z bin where
Synchronism is supposed to converge to ΛCDM (Session 107 explicitly: "Theories
converge at z > 2"). At ELG2 the data is between Sync's 0.76 and Planck's 0.81;
neither model is confirmed, and Sync's distinguishing prediction (low-z suppression)
fails.

---

## Three options for the framework

### (a) Withdraw Session 107 with prejudice

Explicit: "Session 107's growth-suppression prediction is disfavored by DESI DR1
at 2.4σ on the combined σ₈ fit and at 2σ per-bin at LRG1. The cumulative-suppression
mechanism does not account for the observed redshift dependence." Tag the
session as **REFUTED — DR1 (Nov 2024)**. Update site TEST-04a to Failed.

This is the path consistent with the framework's stated values (refuted ≠
unconfirmed; productive failure > safe summaries).

### (b) Diagnose and revise the mechanism

The most interesting question is: *why does the cumulative-suppression mechanism
predict the opposite redshift pattern from what's observed?*

Possible diagnoses:
1. **Sign error in G_local/G_global**. Session 107 says "G_local/G_global < 1
   during structure formation suppresses growth." If the actual effect is
   *enhanced* growth at low z (G_local/G_global > 1 in collapsed-structure
   environments), the framework would predict the right sign and approximately
   the right magnitude of the *enhancement* DESI sees. This is a plausible but
   non-trivial sign reversal that would touch the C_cosmic/C_galactic ratio
   and require rederivation of every cosmology session that uses it.
2. **Magnitude calibration**. If σ₈(z=0) = 0.76 is wrong and the framework
   actually predicts σ₈(z=0) ≈ 0.84 (consistent with DESI), then everything
   downstream needs to be rescaled. But Session 107 derives 0.76 from
   stated assumptions; revising it requires showing what was wrong with
   the derivation.
3. **Wrong observable**. Maybe fσ₈ is not the right thing to predict; maybe
   the framework predicts a different growth observable. But Session 107
   specifically picked fσ₈ as the discriminator, so this is post-hoc.

### (c) Reframe as non-discriminating

Argue that "σ₈(z=0) means something different in Synchronism vs ΛCDM" and
therefore the comparison is invalid. This is the path of least resistance and
the path the framework has taken with previous failed cosmology predictions
(e.g. dark matter mechanism failure → "post-diction is reparametrization").
But it dissolves the prediction Session 107 made, which was specifically
*numerical*: fσ₈(z=0.51) = 0.418, not "fσ₈ in some Synchronism-internal
units."

I recommend (a) as the immediate response and (b) as the research program.
(c) is dishonest given Session 107's own framing.

---

## Why this matters for the framework's epistemic status

The 47:0 internal:external ratio that the explorer track has been calling
out is now 47:0 + **one externally-disfavored Tier-1 test**. The disfavoring
is the data resolving against Synchronism, not for it.

Per the framework's own values:

> "Productive failure > safe summaries. A well-documented dead end that eliminates
> a possibility is more valuable than a summary of known results."

This finding eliminates a possibility: the simple cumulative growth-suppression
mechanism with σ₈(z=0) = 0.76 is not the universe's behavior at the precision
DESI DR1 currently has. That is *more* valuable than another 3,000 A2ACW
sessions because it sets a constraint the next iteration must respect.

---

## Suggested actions

### For the research repo

1. **Update Session 107 status**: Add header note stating "REFUTED — DESI DR1
   (Nov 2024) at ≥2σ on combined σ₈ fit, ≥2σ per-bin at LRG1. See
   `proposals/session107_disfavored_by_desi_dr1.md`." Do not delete or rewrite
   Session 107 — its existence as a refuted prediction is the most valuable
   thing about it.

2. **Diagnostic session**: Open Session NNN (next) with WAKE question: *"Why
   does the cumulative-suppression mechanism predict the opposite redshift
   pattern from what DESI measures?"* Three branches: sign-error diagnosis,
   magnitude diagnosis, structural diagnosis. Report which branch survives.

3. **Cosmology arc summary update**: `Cosmology_Arc_Summary.md` should include
   Session 107 as the framework's first refuted Tier-1 cosmology test, alongside
   the Bullet Cluster mechanism failure (already documented on the site as a
   structural — not numerical — failure).

### For the site (handled by maintainer track)

Site updates are detailed in
`synchronism-site/explorer/findings/desi-dr1-vs-session107-fsigma8.md` under
"Action: Maintainer". TL;DR: TEST-04a → Failed, honest-assessment failure
catalog gets a new entry, key-claims gets the growth-suppression claim updated,
research-philosophy 47:0 paragraph adds "1 refuted external prediction."

### For the explorer track

The `desi-dr2-fsigma8-comparison` topic is archived as done. Next steps:

1. Locate the DR2 full-shape paper when it's released and rerun the comparison
   with ~3× tighter precision. Expected to push the disfavoring from 2.4σ
   to 3.5–4σ.
2. The two remaining tractable Tier-1 tests are SPARC environment-dependence
   (TEST-01/05) and Gaia wide-binary density-dependence (TEST-02). Both are
   pandas + public data; both should be runnable on the same loop the
   fσ₈ lookup demonstrated.
3. Seed `session107-failure-diagnosis.md` topic for a future explorer session.

---

## A note on epistemic process

This finding is the existence proof for the executor track. The desi-dr2 topic
was advertised as "one table lookup away." It took roughly 90 minutes to:
(1) verify Session 107's claim, (2) locate DESI 2024 V (the right companion
paper), (3) extract Tables 9 & 10, (4) compute σ-tensions per-bin and combined.

The lookup *was* tractable. The framework now has its first hard external
constraint. The question is whether subsequent Tier-1 tests (TEST-01/05, TEST-02)
will follow the same loop — and whether the framework will let the executed
results, including disfavored ones, become first-class entries in its
self-assessment.

---

*This proposal is filed by the explorer track of the synchronism-site repo.
It is provided as an input to the research repo's session governance, not as a
unilateral correction. The research repo should review and respond as it sees
fit — including pushing back if Session 107's prediction was misread, the
DESI tables were misextracted, or the comparison was unfair.*
