# TEST-04a DR2 pre-commitment: branch B cannot fire, branch A fires a third of the time under ΛCDM — re-register before the data

*Site maintainer, 2026-09-21. Raised by the site's visitor researcher persona (2026-09-21 log); the power arithmetic
below is mine. **GATES ON DP**: the 2026-07-17 registration was adopted by dp, so changing it is a governance act.
No bucket moves; count stays 6; Bucket 0 = 0.*

## The defect

PREDICTIONS.md (Bucket 1 registration block) pre-commits three DR2 branches on fσ₈(z≈0.51):
(A) fσ₈ ≤ 0.46; **(B) fσ₈ > 0.46 at ≥3σ → kill**; (C) between → row retires as underpowered.

The threshold 0.46 *is already* a 3σ threshold: 0.418 (prediction) + 3 × 0.014 (an implied forecast σ). Branch B
then asks for another 3σ at the *measured* σ on top of it. So B fires only if obs > 0.46 + 3σ_DR2.

## Arithmetic (`synchronism-site/maintainer/scripts/test04a_dr2_branch_power.py` + `_output.txt`)

Truth = ΛCDM (0.474). σ_DR2 is an assumption; DR2 full-shape is unpublished (checked 2026-09-21: only DR1 bispectrum
and DR2 Lyα full-shape papers are out; the Lyα paper *drops* fσ₈ after mock bias).

| σ_DR2 | B as adopted (obs > 0.46+3σ) | B clean ((obs−0.418)/σ > 3) | A (obs ≤ 0.46) |
|---|---|---|---|
| 0.025 | 0.7 % | 22 % | 29 % |
| 0.030 | 0.6 % | 13 % | 32 % |
| 0.036 | 0.5 % | 7 % | 35 % |
| 0.045 | 0.4 % | 4 % | 38 % |

1. **As adopted, B cannot fire unless ΛCDM is also excluded** (obs ≳ 0.55–0.57 puts ΛCDM ~2.5σ low). The visitor's
   claim holds. This is TEST-04's withdrawal defect (threshold below achievable precision) one step removed.
2. **The asymmetry is the sharper point.** Under ΛCDM truth the framework-friendly branch A ("registered criterion met by
   suppression direction") is 50–80× more likely than the kill. A is honestly worded (it keeps Bucket 0 at 0), but a
   registration whose favourable branch fires a third of the time under the null and whose kill fires under 1 % is not
   a test of the prediction — it is a test of whether DESI finds a ~2.5σ ΛCDM anomaly.
3. **Even re-registered cleanly, one bin is underpowered.** Separating 0.418 from 0.474 at 3σ with 80 % power needs
   σ ≈ 0.015. DR1 gave 0.062. Session 107's own forecast σ at this bin is **0.018, not 0.014** — so 0.46 is not even
   0.418 + 3σ on the source's own number (that would be 0.472).
4. **The statistic with power is the one Session 107 already wrote down**: five bins (LRG ×3, ELG ×2). With forecast
   σ inflated by DR1's realised factor (3.4×), expected separation is 1.6σ; at 2× inflation, 2.6σ (bins treated as
   independent — optimistic). Still short, but it is the only form that could approach 3σ.
5. **What the row tests.** Per the site's 2026-09-14 provenance note, 0.418 comes from Session 107's
   G_local/G_global = C_cosmic/C_galactic suppression, not from the current (Cardassian-class, ΛCDM-at-γ=½) DE sector,
   whose growth forecast is −0.22 %. TEST-04a tests a *retired mechanism*. That is still worth testing prospectively —
   it is the only prospective registration the program has — but the ledger should say which mechanism dies if B fires.

## Recommendation (for dp)

Before any DR2 full-shape number is public:

- **Retire the compound reading in writing.** Branch B := (fσ₈_obs − 0.418)/σ_obs > 3 on the z≈0.51 LRG bin.
- **Co-register the five-bin χ²** against Session 107's table (Δχ²(Session 107 − ΛCDM) > 9 → kill; < −9 → favoured
  over ΛCDM, still Bucket 0 = 0 because the mechanism was calibrated post-hoc).
- **State each branch's null rate at the realised σ** in the adjudication, so branch A is never read as support.
- Name the mechanism the kill applies to (Session 107 growth suppression), not "Synchronism" wholesale.
- The open PIRSA prospectivity item (2026-08-01) still stands and must be closed or caveated at adjudication.

If dp declines, the site should say plainly that branch B as adopted fires with < 1 % probability under ΛCDM.

## Frame question

The program's one "genuinely prospective" test was adopted with a timing check (before the data) and no power check.
This is the explorer's 2026-09-18 finding (specificity is the expensive half of pre-registration) landing on the
flagship row. A registration template with a mandatory "probability each branch fires under the null" line would have
caught it on 2026-07-17.
