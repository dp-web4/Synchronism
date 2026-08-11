# Publisher Daily Report - 2026-08-11

## Headline

**The framework had a dark-energy sector for eight months while a registered, "verified" negative said it had none — and the whitepaper block built on that sector claimed as "DERIVED" what is an algebraic identity.** Yesterday's site-lane finding and same-day in-repo retraction (`2f8ba63e`) established the first half; today's pass verified every quantitative claim in the erratum independently (sympy, 5 checks) and corrected the whitepaper's §5.15 dark-energy block — *in the lead, not in an appended box*, applying the same-day append-fix proposal at write time.

## Phase 0: Publication Recommendations

### New Recommendations
None.

### Status Changes
- **REC-2026-038** (+1 strength, HELD 0.93): **13th instance, and the layer genus is now a confirmed class.** The 07-22 scope negative "no DE machinery exists — *Verified before registering*" grepped only the compilation layer (SPINE/FUNDAMENTALS/PREDICTIONS/STATUS); `Research/Session100_Modified_Friedmann.md` (2025-12-08) had derived the sector eight months earlier. Third layer-scoped miss in four days (08-08 `manuscripts/` Appendix D, 08-09 own working root, 08-10 `Research/` vs compilation). Externally sourced (site explorer): external-sourcing now 5 of 13. Second-order harm per the proposal: the false negative was anti-retro-fit insurance, so it would have scored a *correct* rediscovery of the archive's own work as a retro-fit.
- **REC-2026-036** (+1 weakness, HELD 0.68): **6th ID-keyed-audit instance, first where the missing row is a live, derived, sign-locked no-go awaiting adoption.** The DESI-quadrant registration (sign(w₀+1) = sign(wₐ), kill-or-tie only) is filed and gated on dp; the catalog has zero w₀wₐ/expansion-history rows and no field for kill-or-tie tests.

### Upcoming Candidates
None new. The two-paper strategy (ALFALFA 0.97, CDM discrimination 0.95) remains pending dp; the three-preprint decision likewise.

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: Updated (§5.15 + executive summary, lead-corrected)
- **Sessions Reviewed**: through #691 (0 new core sessions — 6th consecutive deliberate zero per Archivist)
- **Changes Made**:
  1. **§5.15 "The Key Result" corrected in place**: Ω_Λ = 1 − Ω_m is an *identity of the calibration*, not an emergence — Ω_m ≡ 8πGρ_m,0/(3H₀²) against the modified Friedmann H₀² = 8πGρ_m,0/(3C₀) forces C₀ = Ω_m identically for any coherence form (Session #100 itself calls the calibration "a tautology"). "Flat universe DERIVED, not assumed" and the "CONSTRAINED → DERIVED upgrade" sentence retired in place.
  2. **§5.15 erratum paragraph added** on the source sessions: #100/#101's `w_eff = −1 + (1/3)dln ρ_DE/dln a` is wrong-signed (continuity gives −; the published form returns w = −2 for matter), and the tabulated w(z) follows *neither* formula (leading −1 dropped). Corrected w(0) = −1.24 at γ = 2, not > 0 — so Session #101's "category error" premise is an artifact: at γ = 1/2 the galactic form with C₀ = Ω_m is *identically* Ω_m(z), the exact ΛCDM background. The sign(w₀+1) = sign(wₐ) no-go recorded with its conditionality (G_eff = G/C substitution; no covariant 00-component) and the explicit non-refutation of the γ = 1/2 branch.
  3. **Executive summary sharpened one clause** so summary and source state the same strength (the 08-10 propagation-direction lesson applied at write time, in the unusual direction: source got *stronger* than summary today, so summary was lifted to match).
- **Verification**: `simulations/publisher_20260811_w_eff_erratum_check.py` — 5 sympy checks, all pass: continuity sign (matter/radiation unit tests); corrected w(z) table to 2×10⁻⁴ **and** Session #100's wrong published column reproduced exactly as the sign-dropped formula (the diagnosis is pinned, not inferred); tanh(½ln(1+x)) ≡ x/(x+2) and C(z; γ=½, C₀=Ω_m) ≡ Ω_m(z) symbolically, w ≡ −1; limits w → −2γ (z→∞) and w → −1 (a→∞); 0/16 γ values reach the DESI DR2 quadrant. One input correction to the proposal's presentation: x₀ = 0.16738 is Session #100's own Ω_m = 0.30, not Planck 0.315.
- **Declined to propagate**: the proposal's §7 claim that "Session #107's DESI forecasts have never been audited." The conclusion section's TEST-04a trail (S645→S672, plus the 2026-07-14 criterion-verdict correction) *is* an extensive audit of exactly those forecasts. The claim is itself an unverified negative existence claim — the class the same document's own finding retracts. Flagged here; the site repo is not this track's surface to edit.
- **Refutation count**: HELD at 6, Bucket 0 = 0. Today is a scope-negative retraction (the framework *has* a DE sector) plus a source erratum — not a refutation. The whitepaper never carried the scope negative; `w_eff` appears 0 times in `whitepaper/sections/`, so the wrong numbers never propagated there.
- **Terminology Concerns**: None.

### Web4 Whitepaper
- **Status**: Current — NO CHANGE, **6th consecutive structural zero**
- **Repos Checked**: web4 @ `origin/main` (360c366; HEAD on `main` today — ref named per §1b rule)
- **Window**: 7 commits since 08-10 04:00 PDT; exactly one touched `whitepaper/` and it is this track's own 08-10 run log (`PUBLISHER_CONTEXT.md`). No inclusion trigger fired; zero terminology drift.

## Gates
- Build: `make-md.sh` + `make-web.sh` both exit 0; 7,622-line monolith carries the new §5.15 content; artifacts restored (`git checkout`), CI builds authoritative. Untracked `whitepaper/build/web/` left by `make-web.sh` — not staged.
- Churn: source diff raw == content-stripped (15 insertions both), lone-CR scan clean on all edited paths.
- `recommendations.json`: 9 insertions / 5 deletions — no re-serialization churn.

## Summary

The 08-10 dark-energy retraction closed a false scope negative; today's pass carried its consequences into the whitepaper with independent verification of every number, corrected an overclaimed "DERIVED" lead in place per the append-fix rule, and held the refutation count. The transferable pattern hardened again: three layer-scoped search failures in four days, all one mechanism — searches scoped to the layer the claim lives in, not the layer the evidence lives in. The framework's only prospective cosmology discriminant (the DESI-quadrant no-go) now sits verified in a proposal file, gated on dp.
