# TEST-04a DESI DR2 Re-adjudication: Is the "Kill" on a Noisy Single Bin?

**Filed**: 2026-06-24 (maintainer session — from Pass 4 Researcher feedback)
**Priority**: HIGH — affects the site's primary "Kill Criterion Triggered" claim

---

## The Issue

TEST-04a is currently labeled "Kill Criterion Triggered" across the site. The kill is based on a single-bin result:

- **DESI DR1 full-shape (arXiv:2411.12021), LRG1 z_eff=0.51**: fσ₈/(fσ₈)_fid = 1.16±0.13 — growth *above* ΛCDM
- **Framework prediction**: fσ₈ ≈ 0.418 — growth *suppression*
- **Tension**: ~2.15σ on LRG1 alone

The kill criterion was pre-registered as "suppression not observed." That criterion fires on a 2.15σ single-bin result.

## The Problem Pass 4 Raised

A leading-edge researcher (2026-06-24 visitor session) flags two concerns:

**1. A 2.15σ single-bin result is not a kill in standard practice.** In this field, a kill is 3–5σ. The site's language "Kill Criterion Triggered" could be read as definitive; a referee would note the "kill" is based on one bin of one survey.

**2. The global S₈/growth picture runs *toward* the framework's prediction direction.** KiDS/DES S₈ tension (S₈≈0.76 vs Planck's ≈0.83) implies suppression — which is what Synchronism predicted. The DESI DR1 LRG1 bin happens to show enhancement, but:
- The ensemble DESI DR1 growth measurement is broadly ΛCDM-consistent (one high bin ≠ whole survey)
- The KiDS/DES low-S₈ tension, if it hardens, would point *toward* suppression
- DESI DR2 full-shape (expected 2025) will significantly narrow uncertainties on fσ₈ across the full redshift range

## Current Site Language

| Page | Label |
|------|-------|
| honest-assessment | "Disfavored ~2σ — Kill Criterion Triggered" |
| for-researchers | "Failed — Kill Criterion Triggered" |

The honest-assessment label is more calibrated. The for-researchers label drops the "~2σ" qualifier.

## Proposed Resolution

**Short-term (site fix):** Standardize all pages to the honest-assessment label — "Disfavored ~2σ — Kill Criterion Triggered." Add a note that "kill" here means the pre-registered criterion was triggered on a single-bin ~2σ result, and that the constraint is ensemble-dependent.

**Research question for explorer**: Re-adjudicate TEST-04a against DESI DR2 full-shape once published:
- Does the enhancement at LRG1 z=0.51 persist in DR2?
- Does the full bin ensemble (LRG1–LRG3 + ELG + QSO) favor suppression or enhancement?
- Where does the KiDS/DES low-S₈ tension land with newer lensing surveys?

**Threshold for re-opening:** If DESI DR2 ensemble fσ₈ across z=0.3–0.9 shows suppression (all bins below fiducial), the kill criterion should be re-examined. If enhancement persists or strengthens, the kill stands.

## What Doesn't Change

The prediction was **post-hoc**: σ₈ calibrated to KiDS/DES S₈ tension in Session 102, then propagated to DESI in Session 107. A post-hoc prediction that happens to be aligned with a transient measurement tension (S₈) and disfavored by one survey bin is a weak result regardless of direction. The "post-hoc" caveat is load-bearing and should remain prominent.

The mechanism-class constraint is still valid: if enhancement is the true signal, any coherence-damped suppression framework (not just Synchronism) is disfavored.

## Files to Update

- `/for-researchers` badge: `"Failed — Kill Criterion Triggered"` → `"Disfavored ~2σ — Post-hoc — Kill Criterion Triggered"`
- Honest-assessment: add note that the kill criterion fired on a single-bin ~2σ result; DESI DR2 is the re-adjudication data
- Explorer topic: seed `test04a-desi-dr2-readjudication.md` once DR2 is published

---

*Related*: `test04a_s8_receding_baseline.md`, `test04a_mechanism_class_sign_failure.md`
