# Phase-2 — special relativity emerges; the model's only novel content is its discreteness (2026-06-22)

**Status:** `[ACTIVE-MRH]` — the most decisive *physics-axis* result of the substrate arc: it
locates, precisely, where any novel prediction must live, and bounds it.
**Sim:** [`simulations/phase2_emergent_dispersion_and_LIV.py`](../simulations/phase2_emergent_dispersion_and_LIV.py) · result: `simulations/results/phase2_dispersion_LIV_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous (dp: "follow what pulls you").

## What pulled me

Everything in the arc so far *reproduces* known physics (breathers, solitons, Madelung,
c-ceiling) — Bucket 0 (zero confirmed novel predictions) unchanged. The deepest question is the
program's central one: **can this model say anything continuum QFT doesn't?** A complex field on
a continuum *is* Klein-Gordon. The one place it can differ is its **discreteness** — so that's
the only place a novel prediction can come from. Measuring the emergent **dispersion relation
ω(k)** tests it directly, and ω(k) *is* the electron→photon curve (rest mass at k=0, light at
high k), so it also tests whether E²=(pc)²+(mc²)² emerges.

## Result

| k | ω_meas | ω_SR (continuum) | ω_lattice | v_g meas | SR dev | lattice match |
|---|--------|------------------|-----------|----------|--------|----------------|
| 0.05–0.2 | = | = | = | = | **<0.02%** | 0.00% |
| 0.6 | 0.774 | 0.781 | 0.774 | 0.74 | −0.9% | 0.00% |
| 1.2 | 1.235 | 1.300 | 1.235 | 0.76 | −5.0% | 0.00% |
| 2.0 | 1.756 | 2.062 | 1.756 | 0.52 | −14.8% | 0.00% |
| 3.0 (→zone edge) | 2.057 | 3.041 | 2.057 | 0.07 | −32% | 0.00% |

1. **Special relativity EMERGES at low energy.** For k≪1 the measured dispersion matches
   continuum SR (ω²=m²+c²k²) to <0.02% — E²=(pc)²+(mc²)² emerges from the substrate, and the
   *same curve* interpolates electron (k=0, ω=m = rest mass) and photon (high k, ω→ck). The
   electron/photon unification from Phase-1.6 is literally this dispersion curve.
2. **The substrate's exact dispersion is the discrete-lattice one** ω²=m²+2c²(1−cos k) (matched
   to 0.00% at every k), which deviates from continuum SR increasingly toward the Brillouin-zone
   edge.
3. **The discreteness signature — the model's only candidate novel content:**
   - a **UV cutoff**: a maximum frequency/energy ω_max ≈ 2.06 at the zone edge (no patterns
     above it — a natural UV completion, no infinities);
   - an **energy-dependent speed of light**: group velocity falls below the SR value as energy
     rises (v_g → 0 at the cutoff) — high-energy patterns propagate *slower*. That is
     **Lorentz-invariance violation (LIV)**, and it maps to a **real observable**:
     energy-dependent photon arrival times from distant gamma-ray bursts / blazars, which LIV
     searches actually measure.

## What this is, honestly

This is the closest the arc has come to a *novel prediction* — and the honest accounting is
that it is **a precisely-located, in-principle-testable prediction class that is neither unique
nor currently distinguishable:**

- **Not unique.** Lattice/discrete dispersion → UV cutoff + LIV is *generic* to discrete-
  substrate theories (loop quantum gravity, causal sets, doubly-special relativity). It is not
  a Synchronism fingerprint — it is the signature of *any* discrete substrate.
- **Planck-scale, hence tiny.** The lattice spacing is identified with the Planck length, so the
  physically relevant deviation scales as (k·a)² = (E/E_cutoff)² — at k=0.2 it is −0.02%, and
  reaches order-1 only at the Planck-energy zone edge. Current GRB photon-timing constraints
  already bound LIV at/near the Planck scale, so the model's predicted deviation sits **below
  current sensitivity.**
- **Bucket 0 unchanged.** Nothing is *confirmed*. This does not move the confirmed-novel count.

**But the advance is real, and it is SPECIFICATION, not confirmation.** The arc converted the
program's vague "underspecified, zero novel predictions" into a precise, defensible statement:

> Synchronism's novel physics, *if it has any*, lives **exactly in the discreteness of the
> substrate**, manifests as **a UV energy cutoff and energy-dependent light speed (LIV)**, sits
> at the **Planck scale**, and the **only experiment class that could ever test it is
> astrophysical photon timing** — where it is currently indistinguishable from generic
> discrete-gravity theories.

That is the most useful thing the physics axis can honestly say right now: not "we predicted
something," but "we now know precisely where our one possible prediction must be, what form it
takes, and what could falsify it." For a research program, knowing *where the novelty must
live* is a real coordinate, not a null.

## The arc, end-to-end (this session)

- **Stage-1:** monotonic saturation can't self-confine; a *focusing* nonlinearity makes a
  standing pattern (electron-like). [Foundation 3 refined.]
- **Phase-1.5:** momentum can't be *boosted* into a pattern; it must be intrinsic.
- **Phase-1.6:** a *complex* field carries momentum in its phase — the three momentum mechanisms
  unify (Madelung), wave/particle duality IS the mechanism, electron/photon from one rule, c as
  ceiling.
- **Phase-2 (this):** SR emerges; the model's only novel content is its discreteness (UV cutoff
  + Planck-scale LIV, GRB-testable, generic-to-discrete).

The substrate model is now concrete: **a complex Intent field with a focusing nonlinearity on a
discrete grid.** It reproduces a startling amount (rest mass, propagation, c-ceiling, relativistic
dispersion, electron/photon unification) — which is exactly why it's *not* novel physics — and
its single novel seam is now named and bounded.

## Honesty notes / caveats

- Agent-zero: a nonlinear mass-renormalization "bonus" measurement was attempted and
  **discarded** — the linear k-mode method is invalid at soliton amplitudes (it returned a
  meaningless −95%). A proper nonlinear-normal-mode method is deferred; no number reported.
- 1D; the dispersion is measured in the linear (small-amplitude) regime; the lattice
  spacing↔Planck identification is the framework's assumption (the spacing is a free parameter
  here). The −93% "deficit at high k" in the JSON is the zone-edge (Planck-energy) extreme, not
  the accessible-energy deviation (which is (E/E_Planck)²-suppressed).
