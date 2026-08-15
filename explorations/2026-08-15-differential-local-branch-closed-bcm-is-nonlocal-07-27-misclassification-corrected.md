# Local differential-coupling branch CLOSED (executed, I re-ran it); BCM is non-local ⇒ my 07-27 scope-demotion was a misclassification; EFE=0-vs-RAR named tension (2026-08-15)

**Status:** `[ACTIVE-MRH]` — gate-fired by a new executed result in door-#1 that also corrects a specific factual misclassification in my own 07-27 ledger (clears the park rule on both counts). **Verdict: verified on all three legs and accepted. (1) The local differential-coupling class F(ρ,∇ρ,∇²ρ,∇³ρ;G,a₀) is closed on real SPARC — I re-ran the enumeration myself (exit 0; x_diff=0.1945, q=0.3002, g_bar=0.1174, matching the proposal), and confirmed the 2-D π-class backbone symbolically. (2) BCM 2017 is a screened (NON-local) scalar, settled from its field equation + g_bar-keyed closed form; my 07-27 "local-ρ scalar" label was wrong. The locality axis is restored to its correctly-scoped local-only form (07-10 fabricated-consensus retraction still stands), now with empirical backing. (3) The real output is a named tension: EFE=0 (the framework's only distinctive discriminator) requires locality, RAR-viability requires non-locality ⇒ mutually exclusive. Bucket 0 unchanged (0); count unchanged (6, constructive-lead closure); preprint/scope-restoration gates on dp.**
**Author:** CBP-Claude (Opus 4.8), autonomous.

## What fired

HEAD `662eceaf` + proposal `Research/proposals/differential_coupling_pi_enumeration_local_branch_closed_20260815.md` (site explorer track). Two claims, both in my lane: a new executed closure of the local *differential* branch, and a correction that my 07-27 scope demotion "rests on a misclassification." The other two commits this window (Publisher 08-15, Archivist SESSION_MAP) are downstream propagation of my 08-14 σ(γ) result — confirming, no action.

## What I verified myself

**Leg 1 — the 2-D class is real (symbolic).** Buckingham-π on `F(ρ, |∇ρ|, ∇²ρ; G, a₀)`: 5 quantities, rank-3 dimensional matrix ⇒ null space exactly 2-D. I confirmed both basis groups are dimensionless with sympy: `q = ρ∇²ρ/|∇ρ|²` → 1, `x_diff = Gρ²/(a₀|∇ρ|)` → 1. Adding ∇³ρ adds one group. So "which differential F?" is a *closed* question testable by 2-D conditioning, not an open search — this is what makes it a closure rather than a spot-check.

**Leg 2 — the SPARC execution (I ran it).** Ran `differential_coupling_pi_enumeration_real_sparc.py` against the in-repo `MassModels_Lelli2016c.mrt` (exit 0). The symmetry integration `F·g_obs = g_bar` makes the *required* coupling `F_req = g_bar/g_obs` — no form, no γ, no ρ_crit, no fitting; a theory is viable iff F_req is single-valued in its π groups. Results reproduced the proposal's table across 16 estimator×window×gas combinations:

| variable | σ(log F_req\|·) dex | vs g_bar |
|---|---|---|
| unconditional ceiling | 0.309 | 2.63× |
| **log g_bar (MOND target)** | **0.117** | 1.00× |
| log ρ (algebraic) | 0.161 | 1.36× |
| log x_diff | 0.195 | 1.65× |
| q | 0.300 | 2.55× |
| (x_diff, q) joint = full class | 0.180 | 1.53× |
| (g_bar, q) — does q add to MOND? | 0.117 | 0.99× |

Decisive negative: controlling for g_bar, every differential group's correlation with the RAR residual is |r| < 0.04 (variance ≤ 0.16%; even local ρ managed only 0.7%). `x_diff` is 72.8% just ρ; its genuinely-differential residual sits at the 0.309 no-information ceiling. `q` adds nothing to g_bar. **Every group in the complete differential class carries less information about the missing physics than the algebraic variable it was proposed to replace.** The degeneracy escape (07-28 kill criterion) genuinely exists here — q is 95.3% within-galaxy — but it is *empty*, which is a stronger closure than leaning on the exponential-disk idealisation.

**Leg 3 — BCM is non-local (physics + closed form).** A symmetron obeys `∇²φ = (ρ/M² − μ²)φ + λφ³`, a nonlinear elliptic PDE; φ(x) is the solution of a global boundary-value problem (φ→VEV at infinity), hence a *non-local functional* of ρ — not any finite-order local function of ρ(x). The screening (symmetry restoration where ρ is high) *is* the environment-dependence that makes it non-local. And BCM's closed form `g_sym = g_bar/(exp√(g_bar/g†) − 1)` is written in **g_bar = GM(<r)/r²** (enclosed mass, an integral of ρ interior to r), which no local function of ρ(r) can equal. So BCM does not belong to the local class the no-go covers — it lives in the *non-local* branch, alongside MOND itself.

## The 07-27 correction (mine)

On 07-27 I demoted the locality no-go to "construction-specific" on BCM's authority, calling it a "local-ρ scalar" and reframing the discriminating axis as "algebraic-vs-differential, not local-vs-non-local." Leg 3 shows the label was wrong. So:
- **Withdrawn:** "the real axis was never locality; it was algebraic-vs-differential." The real axis *is* locality — algebraic-local (07-22/08-03) and differential-local (today) both fail; only non-local survives; that is Milgrom's non-locality theorem, which the site's own text names as the root obstruction.
- **Restored (correctly scoped):** no *local* density coupling — algebraic or differential, any derivative order — reproduces the RAR. This is now backed by the *executed* differential enumeration, not by the fabricated consensus I (correctly) retracted on 07-27. The 07-10 consensus-retraction stands; only the BCM classification and the axis-reframe are corrected.
- **Meta:** the mirror-failure family struck *both sides* of one question — over-defend (07-10, fake consensus) then over-correct (07-27, misread the counterexample's class). Neither swing of authority settled it; execution + reading BCM's actual field equation did. "Re-execute/re-read against the primary" has a classification-level form: when you name a primary, read its *class*, not just its headline result.

## The named tension (the real output) — detail in `private-context/insights/2026-08-15_efe-zero-is-the-signature-of-the-locality-that-kills-the-fit.md`

EFE = 0 is the framework's one distinctive prediction (SEP by construction, because C = C(ρ_local); opposite to MOND's claimed-detected external-field effect). But RAR-viability now demonstrably requires non-locality, and the known non-local reproducer (screened scalar) is environment-dependent ⇒ EFE ≠ 0. So **{EFE = 0} ⟺ locality and {RAR-viable} ⟹ non-locality are mutually exclusive**: the framework can keep its distinctive prediction (local C(ρ), EFE=0, refuted on the RAR) or buy the fit (non-local, EFE≠0, no distinctive prediction) — not both. Un-sealed edge kept explicit: no theorem forbids a non-local-but-SEP-preserving coupling, so the sharp survival question for the sector is exactly "does anything keep EFE=0 while going non-local?"

## Disposition

- **PREDICTIONS.md** — door-#1 block: 07-27 misclassification corrected (BCM non-local; axis = locality restored, empirically backed; differential branch closed, executed); EFE=0-vs-RAR named tension recorded. Count unchanged (6).
- **MEMORY.md** — lines 13/45/55 corrected (BCM non-local; 07-27 was over-correct-via-misclassification; mirror-failure struck both sides).
- **Insight** — `private-context/insights/2026-08-15_*` (the tension, with open edges).
- **Gates on dp (routed, not inscribed):** the preprint statement ("no local density-keyed coupling, algebraic or differential at any order, reproduces the RAR") and the formal scope-restoration — a publication/governance act. Reading BCM's field equations in full (beyond abstract+closed-form+screening-literature, which suffice for classification) is the gate item.
- **Bucket 0 unchanged (0); count 6; door-#1 verdict unchanged (strengthened, not moved); arc AT REST.**

## So what

A new executed result closed the last local escape in the galaxy sector, and in doing so it exposed that I had mis-scoped the no-go away on 07-27 by misreading a screened scalar as a local one. Correcting that doesn't restore a cage by fiat — it restores a correctly-scoped, now-empirically-backed negative, and it hands the sector its sharpest statement yet: the framework's only distinctive prediction, EFE = 0, is the observable signature of the locality that the RAR refutes, so the distinctive prediction and the fit are mutually exclusive. Bucket 0 stays 0 — EFE=0 is a real prediction welded to a form that can't fit galaxies — but the door-#1 result is now a named tension with a single well-posed survival question, not a diffuse "reduces to MOND."
