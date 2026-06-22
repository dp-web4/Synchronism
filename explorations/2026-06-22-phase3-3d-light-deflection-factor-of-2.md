# Phase-3 — 3D light propagation & gravity: the factor-of-2 discriminator (2026-06-22)

**Status:** `[ACTIVE-MRH]` — the substrate arc resumes in 3D (retracting a premature "plateau"
read; 1D/2D were artifact-laden napkin sketches). Lands the most famous test in the history of
gravity and extracts a sharp substrate requirement + the real open question.
**Sim:** [`simulations/phase3_3d_light_deflection_factor_of_2.py`](../simulations/phase3_3d_light_deflection_factor_of_2.py) · result: `simulations/results/phase3_light_deflection_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous (dp: "1D/2D are napkin sketches with built-in
impenetrable walls; reality is 3D; degrees of freedom matter. What keeps photons linear in 3D,
and why does gravity bend them, weakly but nonzero? Keep proposing.").

## Retraction first

The prior turn called the physics arc a "plateau." That was wrong — it generalized from
**dimensionally crippled** 1D/2D models, where missing dimensions act as impenetrable walls
(the Peierls–Nabarro pinning that froze the 1D soliton is a 1D artifact; true topological charge
and transverse stability are 3D-only). dp's correction stands: degrees of freedom change the
physics, not just the resolution. The arc is not plateaued; it moves to 3D.

## The two questions, and the proposal

- **Q1 — what keeps a photon straight in 3D?** Nothing *active*. In a homogeneous, isotropic
  substrate a propagating mode goes straight by **symmetry** — translational invariance →
  momentum conservation (Noether); the eikonal/ray limit of a wave in a uniform medium is a
  straight line (Fermat). Straightness is the default — the *absence* of a gradient.
- **Q2 — why does gravity bend it?** Mass = a localized high-Intent core that makes the
  substrate locally inhomogeneous: it lowers the local reconstruction rate *c* (an effective
  refractive index n > 1). A gradient in *c* bends the ray toward the slower region
  (gradient-index optics). This reuses the framework's own "c = pattern-reconstruction rate"
  and "mass = high-Intent core."

## The sharp test — the factor of 2 (history's most famous gravity discriminator)

Effective index n(r) = 1 + α·(GM/c²)/r. Weak-field deflection Δθ = 2α·(GM/c²)/b.
- **α = 1** (scalar, "variable-c only") → Δθ = 2GM/(c²b) = **Newton/Soldner**.
- **α = 2** (GR effective index) → Δθ = 4GM/(c²b) = **General Relativity** (Eddington 1919, observed).

A scalar substrate gives **half** the observed bending — the value that killed scalar gravity.

## Result (3D eikonal ray-trace, geometric units c=GM=1, r_g=1, weak field b≫1)

| α | meaning | Δθ·b (→ 2α expected) | reading |
|---|---------|----------------------|---------|
| 0 | homogeneous (no mass) | **0.000** at every b | **straight line** — Q1 confirmed |
| 1 | scalar / variable-c only | 2.08 → 2.03 → 2.01 → **2.005** | = 2GM/(c²b), Newton/Soldner |
| 2 | GR effective index | 4.34 → 4.13 → 4.06 → **4.025** | = 4GM/(c²b), General Relativity |

**Factor of 2 measured = 2.035** (α=2 / α=1), across impact parameters (converging as b grows
into the deep weak field).

## What it means

- **Q1 ✓:** a homogeneous substrate propagates the ray dead straight; linearity needs no
  mechanism, only symmetry.
- **Q2 ✓:** a local c-gradient bends the ray, weakly, toward the mass.
- **The substrate requirement (the real content), stated faithfully to absolute time:** in
  Synchronism **time is absolute** — the tick is the universal advance of grid state; a pattern
  *in* the grid cannot perturb the substrate clock. So Synchronism is **denied the time half**
  that GR uses. GR builds its factor-of-2 as *half from spatial curvature + half from
  gravitational time dilation*; an absolute-time substrate must obtain the **whole** 2× from
  **spatial/structural** degrees of freedom alone. A mass modeled as *only* a slowdown of Intent
  reconstruction (one structural effect) gives `n−1 = GM/c²r` → **half** the observed bending →
  **falsified**, exactly as scalar gravity was in 1919. To reach GR's `n−1 = 2GM/c²r` the grid
  must supply both halves spatially: the **reconstruction rate** (cells a pattern advances per
  absolute tick — kinematic, *not* a clock rate) **and** the **grid spacing** contracting near a
  high-Intent core, each contributing `GM/c²r`. Two spatial effects, **no clock change**. *(The
  earlier draft of this doc said "perturb time AND space" — that was a category error imported
  from GR's bookkeeping; corrected per dp 2026-06-22.)*

## The open question — the live frontier (keep proposing)

Absolute time splits this into **two** sharp, falsifiable sub-questions:

**(a) Can grid structure alone supply the full 2×?** Both candidate mechanisms are spatial —
reconstruction rate (cells per absolute tick) and grid spacing. If a high-Intent core
*necessarily* contracts the local spacing by the same fraction it slows reconstruction (each
`GM/c²r`), the substrate **derives** gravitational lensing (α=2 falls out, with *no* clock
term). If only one mechanism is forced, a spatial-only substrate tops out at the **half** →
falsified like scalar gravity. Model the core's effect on *both* local cells-per-tick and local
spacing from the substrate rules; check whether the ratio is forced to 1:1 or free.

**(b) Can the framework explain *measured* gravitational time dilation without a clock?**
Pound–Rebka, GPS, and optical clocks at centimeter height differences directly measure clocks
running slow deeper in a potential — and GR attributes *half* the light-bending to exactly this.
An absolute-time substrate must re-derive that as an **instrument/pattern effect**: a pattern
deeper in the well reconstructs through more grid-traversal per absolute tick, so an *embedded*
clock *reads* slow while the universal tick is unchanged. This is the framework's existing
"time-dilation-as-instrument-effect" bet (PREDICTIONS Bucket 1) — now made **load-bearing for
gravity**. If that account reproduces redshift + Shapiro delay + the time-half of deflection,
Synchronism is a genuine same-prediction / different-ontology alternative; if it cannot, absolute
time is **falsified by gravity**. Neither (a) nor (b) is answered — this is where to keep proposing.

## Honesty

NOT novel physics: the "GR as an effective medium with n = 1 + 2GM/c²r" picture, and the
factor-of-2 itself, are textbook. The value is (1) answering Q1/Q2 with a runnable genuinely-3D
model (no missing-dimension walls), (2) extracting the exact substrate requirement *as
constrained by absolute time* — the full 2× must come from grid structure alone, with no clock
term, which is a **stronger** demand than GR's space+time split, and (3) naming the two precise
open questions (can spatial structure alone give the 2×; can measured time dilation be an
instrument effect) that decide whether the substrate *explains* or merely *fits* gravity. Bucket
0 unchanged (nothing confirmed).

Caveats: eikonal / geometric-optics limit (a photon as a ray — appropriate for deflection; no
diffraction modeled); weak field (b ≫ r_g); α is imposed here by hand (resolving whether the
substrate *forces* α=2 is the open question above), and the spacing↔Planck / Intent↔mass
identifications are the framework's.
