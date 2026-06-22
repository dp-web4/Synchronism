# Phase-3c — invert the frame: gravity as substrate inflow (2026-06-22)

**Status:** `[ACTIVE-MRH]` — the strongest result of the substrate arc, and a frame correction
that retracts conclusions from Phase-3 and Phase-3b. An absolute-time, flat-space substrate
**inflow** reproduces the *full* general-relativistic light deflection (`4GM/c²b`), dissolves the
Yukawa obstacle, derives the equivalence principle, and recasts gravitational time dilation as an
instrument effect — all from the framework's own ontology.
**Sim:** [`simulations/phase3c_inverted_frame_substrate_flow.py`](../simulations/phase3c_inverted_frame_substrate_flow.py) · result: `simulations/results/phase3c_inverted_frame_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous — prompted by dp's frame correction.

## The frame correction (dp)

Every prior probe was **pattern-centric** — the car in the wind tunnel. "A photon travels
*through* the substrate." "Mass *sources* a field." "What does the light-pattern *see*." That is
the observer-centric bias of classical science: computationally convenient, observationally
correct (lift/drag; deflection angles), but **inverted**. The fluid is *there*; the vehicle moves
through it. In race-car CFD, still air is **drawn toward and under** the approaching car —
understanding *that* (not "air moves around the car faster") produced real efficiency
breakthroughs. The car-centric description computes correctly and hides the mechanism. So did
epicycles.

Applied here: **mass is not a source emitting a field. It is a convergence the substrate flows
into.** The Intent substrate is drawn *in* toward a coherent region. Gravity is not the mass
reaching out; it is the substrate flowing in, and patterns are carried by the flow.

## What inverting dissolves and derives (before any number)

- **The Yukawa obstacle (Phase-3b) evaporates.** That obstacle came from modeling mass as an
  *emitter* — an emitted massive field decays as `e^{−w₀r}`, short-range. An **inflow** has
  nothing to emit and nothing to attenuate; its profile is set by flux through 3D shells, a
  **power law — long-range**. The short-range problem only ever existed in the pattern-centric
  frame. It was an epicycle.
- **The Gullstrand–Painlevé profile is derived, not fitted.** Requiring that a pattern carried by
  the flow self-advects at Newtonian `−GM/r²` (i.e. `(u·∇)u = −GM/r² r̂`) forces
  **`u(r) = √(2GM/r)` inward** — the escape-velocity / river profile. *(Verified: `(u·∇)u / (−GM/r²)
  = 1.0000` at r = 10, 40, 160.)*
- **The equivalence principle is free.** Free-fall = floating in the river = no force felt. GR's
  deepest feature; the pattern-centric frame has to impose it.

## The number (computed, not asserted)

Light swims at speed `c` **relative to the local (moving) substrate**. Eikonal Hamiltonian of
light in a moving medium (drag coefficient 1 — the substrate *is* the medium):

```
H(x,k) = c|k| + k·u(x);   dx/dt = c·k̂ + u;   dk/dt = −∇(k·u)
```

Deflection coefficient `Δθ·b/(GM/c²)`, measured from the **final wavevector** (clean; the
group-velocity estimate carries residual flow at finite box size), sweeping `b` and extrapolating
`b→∞` with a `1/√b` fit (the natural strong-field correction for a `√(2GM/r)` flow):

| b | coeff (from k) |
|---|---|
| 60 | 4.290 |
| 120 | 4.178 |
| 240 | 4.124 |
| 480 | 4.097 |

**Fit: `coeff = 3.979 + 2.34/√b` → `b→∞` limit = 3.98 ≈ 4.** (residuals < 0.015)

**The pure absolute-time inflow alone reproduces the full general-relativistic light bending
`4GM/c²b`** — not the Newtonian/Soldner half (2) I'd have bet on.

## Why — and how it unifies dp's two corrections

The eikonal `H = c|k| + k·u` defines an effective (acoustic / Gordon) metric whose time
component is

```
g_tt = −(c² − u²) = −(c² − 2GM/r) = −c²(1 − 2GM/c²r)   ← exactly Schwarzschild's g_tt
```

So the metric coefficient GR attributes to **gravitational time dilation** is here **manufactured
by the flow's `u²` term** — *not* by any clock slowing. The substrate's tick stays absolute; a
clock deep in the well *reads* slow because it is being **carried by the moving substrate**, an
**instrument/pattern effect**. This is exactly the framework's standing
"time-dilation-is-an-instrument-effect" bet (PREDICTIONS Bucket 1) — and it is now **vindicated as
the mechanism behind the full GR deflection**. Absolute time and substrate-centric flow are not
two separate corrections: **the flow generates the apparent time dilation, and full GR bending
follows.**

## Retraction (Phase-3)

My Phase-3 conclusion — "absolute time forbids the time-dilation half, so the substrate must
*also* slow `c` near mass (a separate 'spatial piece') to climb from 2 to 4" — was **itself
pattern-centric and is refuted**. Adding a `c`-gradient on top of the flow **overshoots**:
FLOW+slow gives coefficient ≈ 6.2 (`α=1`) and ≈ 8.4 (`α=2`). There is **no separate piece to
add**; the flow's `u²` already *is* the effective `g_tt`. The factor-of-2 does not "split into
space + something" — it is one mechanism, the inflow.

## Honesty

**Not novel physics.** The river model / Gullstrand–Painlevé coordinates / analog ("acoustic")
gravity are a known way to write GR as a flat-space flowing medium; that this flow reproduces GR
light bending is established. **Bucket 0 unchanged — this reproduces GR, it does not beat it.**

**What the result *is*:** the framework's own ontology — absolute time, substrate-centric, Intent
drawn toward coherence — lands *exactly* on the flat-space, absolute-time, flowing-medium
formulation of GR. It dissolves the Yukawa obstacle, derives the equivalence principle and the GP
profile, and recasts gravitational time dilation as the instrument effect Synchronism already
claimed. That is a strong **internal-consistency** result and a vindication of the inverted frame,
even though the physics is not new.

**The substrate-native open question (the real frontier now):** the `√(2GM/r)` inflow was put in
by hand here (justified by self-advection = Newtonian). The open problem is whether the
*substrate's own rules* — Intent flux converging on a coherent region — **force** that profile.
If a coherent Intent convergence necessarily draws substrate at the escape-velocity rate, gravity
is *derived*; if the profile is free, it is *fit*. That is the next thing to actually compute, and
it is now a clean, substrate-internal question — not a borrowed mechanism.

## Caveats

Eikonal / geometric-optics (a photon as a ray; appropriate for deflection). Drag coefficient 1
(the substrate is the propagation medium itself). Planar — a single mass + one ray lies in an
exact symmetry plane of the 3D problem, **not** a dimensional-reduction artifact (no
impenetrable walls). The `√(2GM/r)` profile is imposed (its substrate-derivation is the open
question above). `1/√b` extrapolation, intercept 3.98 with sub-0.02 residuals.
