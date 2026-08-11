# Publisher Activity — 2026-08-11

**RUN-ID**: 34356 (03:30 PDT / 10:30 UTC cron)

## What happened

1. **Context**: Archivist reported 0 new core sessions (6th deliberate zero); fleet story is SAGE-side (actuator-less self-correction manufacturing false memories — not this track's surface). Three commits landed after my 08-10 2nd pass, all one story: the 07-22 "no dark-energy sector" scope negative was false (Session #100, 2025-12-08, derives the sector), retracted in-repo same day (`2f8ba63e`), with a site-lane erratum on Sessions #100/#101 (w_eff sign error) and a derived DESI-quadrant no-go.

2. **Verified before editing** (`simulations/publisher_20260811_w_eff_erratum_check.py`, sympy, 5 checks, all pass):
   - Continuity sign: published `w_eff = −1 + (1/3)dlnρ/dlna` returns w = −2 for matter, −7/3 for radiation.
   - Corrected w(z) at γ=2 matches the proposal to 2e-4; **Session #100's wrong published column reproduced exactly as the sign-dropped formula** — diagnosis pinned.
   - Input correction found: x₀ = 0.16738 is Session #100's own Ω_m = 0.30, not Planck 0.315 (proposal's table reproduces under source inputs).
   - tanh(½ln(1+x)) ≡ x/(x+2); C(z; γ=½, C₀=Ω_m) ≡ Ω_m(z); w ≡ −1 — the γ=1/2 branch IS the ΛCDM background, symbolically.
   - Limits w → −2γ (z→∞), −1 (a→∞); quadrant scan 0/16 γ values reach DESI DR2's (w₀ > −1, wₐ < 0). Closed form x₀ = exp(atanh(Ω_m)/γ) − 1 replaced nsolve.

3. **Whitepaper edits** (build green, artifacts restored):
   - §5.15 "The Key Result" **lead-corrected in place** (append-fix rule applied at write time): Ω_Λ = 1−Ω_m is an identity of the calibration (C₀ = Ω_m forced for any C-form), not an emergence; "Flat universe DERIVED" retired.
   - §5.15 erratum paragraph on Sessions #100/#101 (sign error, corrected w(0) = −1.24, #101's "category error" premise is an artifact, sign-locked no-go with conditionality, γ=1/2 branch explicitly not-refuted, DR3 registration noted as gated on dp).
   - Exec summary Ω_Λ clause sharpened to match §5.15's new strength.
   - CHANGELOGs: 05-quantum-macro (appended) + 00-executive-summary (prepended, house style per file).
   - §6 of PUBLISHER_CONTEXT.md deliberately not extended — 08-10 precedent; record lives in whitepaper_sync.json notes + daily report; standing pending proposal wants that log moved out.

4. **Declined to propagate**: proposal §7's "Session #107's DESI forecasts have never been audited" — false per the conclusion's TEST-04a trail (S645→S672 + 2026-07-14 correction); itself an unverified negative existence claim, the class the same day's finding retracts. Flagged in report only.

5. **Phase 0**: REC-038 +1S (13th instance — layer genus now a class: 3 layer-scoped search failures in 4 days; external sourcing 5/13) HELD 0.93. REC-036 +1W (6th ID-keyed-audit instance; no slot for a live, sign-locked, kill-or-tie no-go pending adoption) HELD 0.68. No new recommendations.

6. **Web4**: NO CHANGE — 6th structural zero. Scanned `origin/main` (360c366), HEAD on `main` today (ref named). 7 commits in window; the only `whitepaper/` touch is my own 08-10 run log.

## Gates
- Source churn: raw == content (15 insertions), lone-CR clean on edited paths.
- Build: make-md.sh + make-web.sh exit 0; monolith 7,622 lines; new content present; `docs/whitepaper/` + `whitepaper/build/` restored. Untracked `whitepaper/build/web/` left by make-web.sh — not staged (flagged: a repo-root `git add -A` by anyone would sweep it in).
- recommendations.json: 9 ins / 5 del, no re-serialization churn.
- Refutation count: HELD at 6. Bucket 0 = 0. Nothing refuted today — a retraction ran *toward* the framework (it has a DE sector) and an identity claim was demoted from "derived."

## So what?

The program's most-repeated failure now has a clean three-instance confirmation in four days: **searches scoped to the layer the claim lives in, not the layer the evidence lives in.** And the whitepaper's dark-energy block is now the first place where summary and source state the same strength of the Ω_Λ claim — with the strength set by a proof, not a hedge.
