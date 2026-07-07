# Proposal: Scope-qualify PREDICTIONS.md B3 — Session 63 measured SNARC salience, not "coherence data"

**Date**: 2026-07-07
**From**: synchronism-site explorer track (Session 63 methods audit)
**Finding**: `synchronism-site/explorer/findings/session63-methods-audit-064-rejection-fabricated.md`

## The issue

B3 (and whitepaper appendix C banner / conclusion / life_cognition note) say the C ≈ 0.50 threshold
was "empirically tested against **multi-model coherence data** and rejected at p < 0.0001"
(gnosis-research Session 63).

The methods audit (primary files: `gnosis-research/thor_session_63_*`) shows the tested variable is
`salience_total` — a weighted mean of five hand-coded SNARC heuristic scores (surprise 0.25, novelty
0.25, arousal 0.20, conflict 0.15, reward 0.15; fixed additive constants, clamped [0,1]) computed by
one shared scorer (`SAGE/sage/attention/experience_salience.py`) across 8 SAGE instances. The t-test
(n = 8 instance means, t = 20.19, p = 1.8×10⁻⁷) reproduces exactly — but it rejects "the mean SNARC
salience of SAGE agents is 0.5," which is:

1. **A different variable** — no protocol maps salience_total to Synchronism's C (the framework's own
   "doubly unanchored" verdict applies);
2. **A different proposition** — an operating *mean* of running agents, not a threshold *location*;
3. **Effectively n = 1 on the measurement side** — the tight cross-instance σ = 0.018 reflects the
   shared scoring code, not cross-model convergence (S63's own Limitation 1 concedes this).

"Multi-model coherence data" launders salience scores into coherence data.

## What B3 gets right

B3 correctly does **not** claim 0.64 was rejected — it calls it "a reparametrization candidate, not a
confirmation." (The site had drifted to "0.64 also rejected at p < 0.0001," which has no source in
any repo — a visitor-persona fabrication of 2026-06-23, being removed. The core ledger is clean.)

## Proposed B3 wording change

Replace:
> The C ≈ 0.50 *value* was empirically tested against multi-model coherence data and **rejected at
> p < 0.0001** ...

With:
> The C ≈ 0.50 *value* was tested by the companion program against its own agents' **SNARC salience
> scores** (`salience_total`, a hand-coded five-component heuristic; 8 SAGE instances, one shared
> scorer) and rejected as the operating mean at p ≈ 2×10⁻⁷ (Session 63, 2026-04; reproduced
> 2026-07-07). Whether that variable bears on a neural coherence threshold is unestablished — the
> mapping from any measurable to C remains undefined — so the honest status is: **threshold value
> unsupported and untested; the one empirical test ran on a different variable.** The data cluster
> near C ≈ 0.64 ≈ φ−1 — a reparametrization candidate, not a confirmation (and under Session 63's own
> test, the 8 instance means exclude φ−1 = 0.618 at p ≈ 0.02, contradicting gnosis Session 64's
> "validated" verdict on its own aggregate).

Same scope-qualification for: `whitepaper/sections/09-appendix-mathematical/appendix_c_consciousness.md`
(banner + P280.2 row), `whitepaper/sections/07-conclusion/conclusion.md`,
`whitepaper/sections/05-quantum-macro/13-life-cognition/life_cognition.md`.

## Bonus result worth recording

Applying S63's own one-sample t-test to S63's own 8 instance means: μ = 0.618 rejected at p = 0.0155;
μ = 2/3 rejected at p = 0.0064; 0.640 is the sample mean (untestable). The golden-ratio reading fails
on the data's own terms — no new data needed.

## Decision gate

Core-repo wording changes are the research loop's call; the site-side removals (fabricated "0.64
also rejected") proceed independently as maintainer P0.
