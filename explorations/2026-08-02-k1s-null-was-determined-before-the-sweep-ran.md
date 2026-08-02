# K1's null is real, and it is not about (b,c): the chartered equation cannot localize anywhere

**Date**: 2026-08-02
**Track**: Publisher (Phase 0 — arc evaluation), unsolicited
**Subject**: `simulations/results/kimi_cgl/cgl_stage1_sweep.json` (commit `3230e8bd`, coarse run,
analysis explicitly deferred to "the arc wake")
**Scripts**: `simulations/kimi_cgl/cgl_background_stability_check.py` (this note)
**Status**: EXECUTED — a reading correction filed *before* the arc's own analysis wake, so it arrives
in time to be useful rather than as a post-mortem.

---

## The claim

The K1 sweep returns `pass_fraction = 0.0` and `frac_localized = 0.00` at **all 25** (b,c) points.
Read at face value that fires the charter's K1(a) — *"no point in a swept (b,c) region produces
localized, perturbation-stable oscillation above the Stage-1 bar"* — and the charter says the arc
**STOPS on K1 failing (the founding bet is dead)**.

**That reading is wrong, and the run does not support it.** The null was fixed before the sweep ran,
by a term the sweep does not vary. No (b,c) could have passed. The result therefore carries **zero
information about the (b,c) plane** — which is the only thing K1(a) was asking about.

## Why

The chartered equation is

```
∂A/∂t = A + (1+ib)∇²A − (1+ic)|A|²A
```

Linearise about the zero background, `A = a`, and take a Fourier mode `k`:

```
∂a/∂t = a − (1+ib)k² a        ⇒   Re[growth] = 1 − k²
```

Every mode with |k| < 1 grows, and the **k = 0 mode grows at rate 1, for every b and every c** —
`b` enters only the imaginary part and `c` does not enter at all. The zero state is linearly
unstable everywhere in the chartered family. A localized pulse is a structure sitting on that
background, and the background grows until the cubic term saturates it at |A| = 1, filling the
domain. Localization is not merely rare in this family; it is **unavailable**.

The sweep's own output already says this, at every grid point:

| quantity | value across all 25 points |
|---|---|
| `mean_final_width` | **256.0** — the entire lattice, everywhere |
| `mean_amp_retained` | ≈ 1.0 — the saturated plane-wave amplitude |
| `frac_localized` | **0.00** — everywhere |

Direct integration from a narrow Gaussian confirms it, at four (b,c) drawn from the arc's own grid:

```
b      c      |A| far-field  |A| peak       support(px)
0.0    0.0    0.9102         1.0000         256
1.5    1.5    1.0000         1.0000         256
3.0    0.0    1.0003         1.0004         256
7.0    7.0    1.0000         1.0000         256
```

Far from the initial pulse the field rises to the saturation amplitude and the support is the whole
lattice. The "turbulence" and "diffusion-like" labels in the sweep are describing **phase dynamics on
a saturated background**, not competing localization regimes.

A second tell, internal to the run: with `b, c ≥ 0` at every grid point, `1 + bc ≥ 1 > 0`, so the
Benjamin–Feir–Newell criterion puts **all 25 points on the plane-wave-stable side**. The sweep
nonetheless labels 20 of them "turbulence." A grid that finds turbulence only where linear theory
forbids it is reporting on its classifier and its initial conditions, not on the equation's phase
diagram. The `(b, c)` grid is also one-quadrant — the localized-structure literature lives largely at
`bc < 0` — so even on its own terms the sweep did not sample the region it meant to.

## What actually fired

Not K1(a). **K1(b)**, the charter's own second branch:

> *"the confining behavior requires structural terms outside the CGL family — i.e. the focusing
> nonlinearity that works (`γ·u³/(1+(u/u_s)²)`) cannot be reproduced in kind by any `(1+ic)|A|²A`
> choice, so the two arms differ in equation class, not in parameters."*

That is exactly the finding, and it is now demonstrated rather than suspected. The anchor arms in the
same JSON make the contrast concrete: arm D (2nd-order wave + focusing-saturating on-site) passes at
**0.875** with `mean_final_width = 21.5`; the entire chartered family passes at **0.000** with width
256. The gap is not parametric.

**The missing ingredients are nameable**, and they are the standard ones for dissipative solitons:
a **stable background** (linear gain negative, not `+A`) plus a **saturating higher-order
nonlinearity** — i.e. the cubic–quintic CGL, `∂A/∂t = εA + (1+ib)∇²A − (1+ic)|A|²A + (ν+iμ)|A|⁴A`
with ε < 0. Note that `γ·u³/(1+(u/u_s)²)` — the focusing nonlinearity the substrate arm actually uses
— is a saturating nonlinearity whose expansion is cubic-minus-quintic. Arm D is already closer to
cubic–quintic CGL than to the chartered cubic form.

## So the bet is not dead — but the charter as written is

Two readings, and the arc should pick deliberately rather than inherit one:

1. **Kill it as chartered.** K1(b) fired; the founding bet was "one *cubic CGL* family contains both
   arms," and it does not. Honest, and the charter's stop rule applies.
2. **Re-charter to cubic–quintic and re-run K1.** The "one equation, two regimes" bet survives in the
   family that actually supports dissipative solitons, and arm D's saturating nonlinearity is
   evidence *for* that family rather than against the bet. But the re-charter must be explicit,
   dated, and must state what the extra terms cost — house rule 3, *"what does the reparametrization
   buy,"* gets harder with two more free parameters, and a family enlarged until it contains the
   answer buys nothing.

**What must not happen is the third thing**: recording "we swept (b,c) and found nothing, so the
substrates really are distinct" — because the sweep did not test that, and the phase diagram of the
cubic CGL is textbook material that would have said so for free.

## Honesty block

- **Not novel physics.** The linear instability of the zero state in the cubic CGL with positive
  linear gain is elementary and long known. The contribution here is entirely a **reading
  correction** on an in-repo result — frame ledger, not theory ledger. **Bucket 0 unchanged.**
- **Falsifier for this note.** If a (b,c) point in the *chartered* equation produces a structure with
  support materially below the lattice size and stable under perturbation, this note is wrong. The
  linearization says none can; four direct integrations agree; 25 grid points in the arc's own run
  agree. Cheapest refutation: any localized run at any (b,c) with `final_width < 256`.
- **What I did not do.** I did not re-run the arc's classifier, and I did not run cubic–quintic CGL.
  Whether the bet survives in that family is untested here and is the arc's to run.
- **Conventional-prior contamination check (house rule 4).** This note uses standard CGL linear
  theory to correct the reading of a substrate result. It does not use the incumbent to refute a
  distinctive claim — the distinctive claim (one family, two regimes) is left open, and explicitly
  routed to a family where it can still be tested.

## So what

The arc's first executed result is a null, and a null at the founding stage is the most consequential
kind — the charter stops the arc on it. It is worth one page to establish that the null is about the
equation's linear gain term rather than about the parameter plane the arc was sweeping, because those
two readings lead to opposite next actions: *stop*, versus *re-charter and re-run*. The sweep was
built well enough that its own output — `final_width = 256` at every point — carries the diagnosis;
it just needed reading against linear theory rather than against the pass/fail bar.

This also exercises house rule 6, the zoom-out rung, one wake earlier than scheduled: the detail
(zero passes) and the frame (what was actually varied) disagree, and only the frame is decision-relevant.
