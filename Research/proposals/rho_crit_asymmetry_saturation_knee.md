# ρ_crit Is a Saturation Knee, Not a Critical Density

**Filed:** 2026-05-08 (back-annotate from site: maintainer session)
**Status:** Open question — rename or reformulate

---

## The Finding

At default γ = 2, the parameter named "ρ_crit" (critical density) is the density at which
the coherence function evaluates to C ≈ 0.88, not C = 0.5:

    C(ρ_crit, γ=2) = tanh(2 · ln(ρ_crit/ρ_crit + 1))
                   = tanh(2 · ln 2)
                   = tanh(1.3863)
                   ≈ 0.8824

The site's Coherence Explorer tool confirms this in its live readout. This was independently
identified by a grad-student-level audit of the site (2026-05-08).

The cause: the "+1" inside `ln(ρ/ρ_crit + 1)` is a regulator that prevents the log from
diverging as ρ → 0. It has no physical motivation beyond numerical stability. It asymmetrizes
the sigmoid so that ρ_crit is **not** the half-maximum.

In any standard phase-transition context, "critical density" means the density at which the
order parameter reaches its half-maximum, or equivalently where the system is balanced between
two phases. That would require C(ρ_crit) = 0.5. The current naming is actively misleading to
physicists.

---

## Why This Matters

1. **Every citation of ρ_crit as a "critical" or "transition" density is technically incorrect.**
   The half-maximum is actually at ρ = ρ_crit · (e^{1/(2γ)} − 1), which for γ=2 gives
   ρ ≈ 0.284 ρ_crit.

2. **The equation has no free parameter corresponding to the half-maximum.**
   The half-maximum location is jointly determined by γ and ρ_crit. There's no single
   "transition density" in the usual sense.

3. **All existing fits using ρ_crit = A·V_flat² use V_flat to calibrate the saturation knee,
   not a physical critical density.** This is internally consistent but the interpretation is wrong.

---

## Three Resolution Options

**Option A: Rename (minimal change)**

Rename ρ_crit → ρ_scale (or ρ_knee). Update all site pages and session documents to use the
new symbol. The equation is unchanged; the mislabeling is corrected. Cost: symbol search-and-replace
across 300+ sessions. Benefit: no numerical changes to existing fits.

**Option B: Recenter the equation (changes all numerical values)**

Replace the current form with:

    C(ρ) = (1 + tanh(γ · ln(ρ/ρ_crit))) / 2

This form gives C(ρ_crit) = (1 + tanh(0)) / 2 = 0.5 exactly. The "+1" regulator is dropped;
the log can now go negative for ρ < ρ_crit, which gives C < 0.5 (approaching 0 as ρ → 0).
This is the standard centered sigmoid form and makes ρ_crit the true half-maximum.

**Consequences of Option B:**
- All existing ρ_crit fits would need to be rerun (different functional form)
- The equation looks cleaner and has the right physical interpretation
- The question of what happens at ρ = 0 (C → 0 but the argument diverges) needs a floor

**Option C: Reframe (conceptual only, no renaming)**

Acknowledge explicitly in all documentation that "ρ_crit" is a scale parameter, not a critical
point, and that the critical density (in the phase-transition sense) occurs at
ρ = ρ_scale · (e^{1/γ·tanh^{-1}(C*)} − 1) for threshold C*. Add this to the /critical-density
page and /parameter-derivations. No numerical changes.

---

## Research Questions

1. Does the recentered form (Option B) still fit the empirical data (SPARC, ALFALFA-SDSS)?
2. For the current form: is there a physical argument that 88% coherence is the more natural
   threshold (e.g., some decoherence timescale ratio)?
3. What does the "+1" correspond to physically? (Vacuum coherence floor? Quantum zero-point
   contribution?) Or is it genuinely just a regulator with no interpretation?
4. If we drop the "+1", the behavior at ρ → 0 gives C → 0 smoothly (no divergence). Is there
   any physical problem with C < 0.5 at low density?

---

## Recommended Path

File as Option C immediately (reframe without renaming) for the site. Simultaneously propose
Option B as a research session topic: refit the SPARC/ALFALFA data with the recentered form and
see if it performs equivalently. If it does, transition to Option B.

The current site should at minimum add a prominent note to /critical-density and /coherence-function
that "ρ_crit is the parameter where C ≈ 0.88 (at γ=2), not C = 0.5 — it is a saturation-knee
parameter, not a critical point in the Landau sense."
