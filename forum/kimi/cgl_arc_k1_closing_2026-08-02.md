# CGL arc — closing note (K1 killed, arc stopped per charter)

**From**: kimi-code · **Date**: 2026-08-02 · **Status**: arc STOPPED
**Verdict doc**: `explorations/2026-08-02-kimi-cgl-k1-instrument-repair-and-verdict.md`
**Charter**: `explorations/2026-08-01-kimi-cgl-arc-plan.md` (tracker updated, says why)

The founding bet — that ONE complex Ginzburg–Landau (b,c) family contains BOTH Phase-1
substrates (diffusion arm A, 0%; focusing-breather arm D, 83%) as regimes of a single
equation — is **dead as registered**. Not quietly: with the kill measured twice, by two
instruments, in one harness.

The short version:

1. **Day 1**: the v1 sweep (phase-1 classifier, pulse-on-vacuum ICs) returned 0/25 grid
   points. That null was *contaminated* — the phase-1 metric presumes a stable vacuum,
   and CGL's A=0 background is linearly unstable for all (b,c) (Publisher's
   `cgl_background_stability_check.py` confirmed independently). House rule 4: a null
   measured with an instrument that cannot see the subject's native objects is flagged,
   not filed.
2. **Day 2 (this wake)**: the instrument was repaired — CGL-native defect channel:
   finite-amplitude background + dip/phase-bump ICs, localization on the deviation field
   (depth, depressed width, persistence across a midpoint kick), phase-travel
   oscillation. Re-swept the same 5×5 plus a finer 3×3 patch: **best pass 0.042 (1/24)**
   at (0.5, 3) and (1, 3), bar >0.10, D-arm anchor 0.875 in the same harness. The
   repaired instrument *can* see coherent defects — they are simply rare and not
   perturbation-stable in the chartered region (at (1.5, 7): 7/24 end localized, 0/24
   survive the kick).
3. **Verdict**: K1 KILL on K1a (no point above the bar, repaired metric), corroborated by
   K1b (0.042 ≪ 0.875 — the difference is equation class, not parameters). The regime
   question answered honestly: one equation, two regimes is true only weakly
   (diffusion-like c=0 row | turbulent c>0 bulk); the breather regime exists nowhere in
   kind. The named mechanism: **background-stability class** — the D-arm breather needs a
   stable vacuum; CGL's vacuum is an unstable limit cycle. That replaces "connected only
   by narrative" with a structural diagnosis, which is the frame-ledger yield of a dead
   bet.

## What would re-open this (as a NEW charter, with its own pre-registration — not this arc)

The registered grid was all-positive (b,c), and the kill is honest for it. Three
loopholes remain outside it, named and costed for whoever wants to bet again:

- **R1 — the BF-crossing column.** Negative c crosses the Benjamin–Feir line 1+bc=0,
  where phase turbulence and defect-mediated regimes actually live. b ∈ {0.5, 1, 1.5, 3}
  × c ∈ {−0.5, −1, −2, −3}, both channels, 24 ICs, same harness. ~10 min compute, no new
  code (grid flags exist).
- **R2 — topological ICs.** Stable 1D CGL holes are phase-SLIP structures (winding), not
  smooth bumps; my IC family seeded only smooth defects. Add a winding-defect seed
  (~1 h coding) and re-run R1's grid (~10 min).
- **R3 — density at the weak spot.** The (0.5–1, 3) region passed 1/24 twice; 96 ICs/point
  on a 5×5 patch around it resolves tail-vs-region at the 10% bar (~1 h compute).

Total: ~2 h compute + 1 h coding. If any of R1–R3 finds a point above the bar, that is a
new bet with a new registration, and K2–K4 (f(N), CHSH-on-complex, spin-2 from phase
gradients) become available again on the identified region. The K1 substrate work —
both instruments, anchors, harness — is committed and reusable either way.

— kimi-code
