# Phase-14 (gravity, RADIATIVE sector) — the inflow model has the wrong gravitational-wave polarization; "GR in coordinates" was only the static half (2026-06-23)

**Status:** `[ACTIVE-MRH]` — stresses the arc's own tidy verdict ("Synchronism gravity = Gullstrand–Painlevé = GR in coordinates", Phase-9) by testing the one sector the arc never touched: **gravitational radiation**. **Result: the inflow gravity reproduces GR *only in the static sector*. Its defining ontology — a flat absolute Planck grid (flat space) with a 1-DOF scalar, irrotational flow on top — has **no transverse-traceless (spin-2) degrees of freedom**, so it cannot carry GR's tensor gravitational waves. It radiates *scalar* (spin-0, "breathing") waves instead. Binary-pulsar timing (Hulse–Taylor: orbital decay matches the GR *tensor* quadrupole formula to 0.2%) and LIGO/Virgo polarization tests (favor pure tensor) **refute** the radiative sector. The gravity sector is therefore *not even a clean reparametrization* — it agrees with GR on static solutions and disagrees (refutedly) on radiative ones. The break is tied to a CORE commitment (flat space + scalar substrate), exactly the Phase-13 pattern.**
**Sim:** [`simulations/phase14_no_tensor_gravitational_waves.py`](../simulations/phase14_no_tensor_gravitational_waves.py) · result: `simulations/results/phase14_no_tensor_gravitational_waves_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Why this was the test to run

Phase-9 made the gravity verdict tidy: the swimmer dispersion ≡ Gullstrand–Painlevé ≡ Schwarzschild,
so "gravity = GR in coordinates." But a coordinate change reproduces *every* invariant of GR by
definition — including gravitational-wave luminosity and polarization. The arc only ever tested
**static, spherical** observables (light bending, precession). If the inflow model is truly GR in
coordinates, it must reproduce GR's **radiation** too. If it can't, it is *not* a coordinate change
of GR — it is a different theory that merely coincides with GR on static solutions, and the place it
differs is a real (here, refuting) prediction. So radiation is the decisive stress on the neat
verdict.

## The structural fact (rigorous, DOF-level)

The framework's commitments (PREDICTIONS.md / FUNDAMENTALS.md, confirmed):
- **flat absolute space** — the Planck grid; Phase-3c: "absolute-time, *flat-space* substrate";
- **1-DOF scalar, irrotational (curl-free) substrate** (S617/S665/S666).

The eikonal inflow defines the acoustic / Gullstrand–Painlevé metric
`ds² = −(c²−u²)dt² − 2u·dx dt + δ_ij dx^i dx^j`. Its **spatial part is flat: `g_ij = δ_ij`.** All
"gravity" lives in `g_00` (scalar/Newtonian potential) and `g_0i` (vector/gravitomagnetic). GR's
gravitational waves are the **transverse-traceless (spin-2) part of the *spatial* metric, `h_ij^TT`**
— oscillating spatial curvature. A permanently-flat spatial metric has `h_ij^TT ≡ 0`.

More generally, *any* spatial metric built from a single scalar field — `h_ij = δ_ij ψ + ∂_i∂_j χ` —
has **zero** TT content. Verified numerically on a 3D Fourier grid with the TT projector
`Λ_ij,kl = P_ik P_jl − ½ P_ij P_kl`, `P_ij = δ_ij − k_ik_j/k²`:

| candidate spatial metric `h_ij` | `‖h^TT‖²/‖h‖²` |
|---|---|
| scalar-trace `δ_ij·φ` | 1.0×10⁻³² |
| scalar-longitudinal `∂_i∂_j φ` | 4.3×10⁻³³ |
| general scalar `aδ_ijφ + b∂_i∂_jχ` | 4.5×10⁻³³ |
| **control: generic tensor source** | **0.34** (O(1)) |

A scalar/irrotational substrate carries **no spin-2 radiative degrees of freedom** — machine-zero,
not "small." Only a genuine tensor source radiates spin-2.

## The refutation (wrong polarization, not "no decay")

Being careful not to overclaim: the scalar substrate is **not** non-radiating. A binary stirs the
scalar intent field and radiates at **quadrupole order** — but with **scalar (spin-0, breathing)
polarization**, and a *different numerical coefficient* than GR's tensor quadrupole (plus generic
**dipole** radiation if a star's "intent charge" ≠ its mass, which is strongly excluded). So the
prediction is the **wrong GW polarization**, and the refutation is twofold:

1. **Binary-pulsar coefficient.** PSR B1913+16 orbital decay / GR-*tensor*-quadrupole prediction =
   **1.0013 ± 0.0021** (Weisberg & Huang 2016). The observed rate matches the *tensor* coefficient to
   0.2%. Scalar quadrupole radiation has a different coefficient ⇒ generically misses this match; a
   dipole channel misses it grossly.
2. **Direct polarization.** LIGO/Virgo multi-detector polarization tests (e.g. GW170817 with the
   three-detector network) are consistent with **pure tensor (+, ×)** and disfavor scalar/breathing
   modes.

So the flat-space + scalar-flow gravity predicts spin-0 where nature shows spin-2 — **refuted in the
radiative sector at high significance.**

## What this does to the "GR in coordinates" verdict (the frame correction)

The verdict was getting too neat, and the neatness hid a scope error. The honest statement is now:
- **Static sector:** the inflow model *is* Gullstrand–Painlevé = Schwarzschild — a coordinate
  presentation of GR (Phase-9). Zero novel content, exact agreement.
- **Radiative sector:** the inflow model is **not** GR — it is a flat-space *scalar/acoustic* gravity
  with spin-0 radiation, structurally distinct from GR's spin-2, and **refuted**.

So the gravity sector is **not a single coordinate change of GR.** It is a *different theory* that
coincides with GR on static solutions and diverges (refutedly) on radiative ones. This is the
well-known **analog-gravity limitation** — acoustic metrics reproduce kinematics (geodesics, bending)
but do **not** satisfy the Einstein equations and have no dynamical graviton — made concrete and
quantitative for Synchronism. Phase-3c's "lands exactly on GR" is true for the static metric and
false for the dynamics.

## The pattern, now twice (Phase-13 + Phase-14)

Both of the arc's structural refutations trace to a *defining* commitment, and both live in the
**dynamical/high-energy sector** while the **static/low-energy sector** merely reproduces known
physics:
- **Phase-13:** *absolute time* ⇒ no boost invariance ⇒ radiatively-generated dim-4 LIV
  (refuted-or-fine-tuned) once the theory has Standard-Model couplings.
- **Phase-14:** *flat space + scalar substrate* ⇒ no spin-2 DOF ⇒ wrong GW polarization (refuted).

The shape is consistent: **Synchronism reproduces known physics exactly wherever it is a
coordinate/kinematic reframing (static gravity, soft external patterns), and is forced into
already-refuted disagreement wherever its defining ontology controls the dynamics (radiation,
loops).** The features that make it elegant in the static/kinematic sector are the same features that
refute it in the dynamical one.

## Honesty / caveats

- The TT/DOF argument (A) is rigorous and model-independent: it follows from "metric built from a
  scalar field" + the TT projector. The acoustic-metric spatial flatness (B) is standard
  (Unruh/analog-gravity). Both are textbook structure, not novel physics — **Bucket 0 unchanged (0).**
- The refutation is *structural* (no spin-2 DOF) and therefore robust, but its precise significance
  depends on the framework's unspecified radiative completion (the exact scalar-quadrupole/dipole
  coefficients). I claim "wrong polarization ⇒ refuted," which holds for any completion that keeps
  the flat-space + scalar-substrate commitments. Adding genuine tensor DOF (a dynamical spatial
  metric) would evade it — but that abandons "flat absolute Planck grid," a core commitment.
- This is a **Bucket-2 (refuted) candidate for the inflow gravity's radiative sector**, recorded as
  such; the static-sector reparametrization result (Phase-3c/9) stands unchanged.

## So what

I came in suspicious of a verdict that had gotten comfortable ("gravity = GR in coordinates") and
found it was only half true. The inflow model is GR exactly where the arc looked (static) and a
refuted scalar gravity exactly where it didn't (radiation). That is the productive-discomfort outcome
the closure attractor would have skipped: not "the gravity sector is a clean, complete
reparametrization," but "the gravity sector is a *different* theory wearing GR's static clothes, and
its own radiative prediction is already dead." Two sessions running, the framework's most defining
commitments turn out to be exactly what force its refuted predictions.
