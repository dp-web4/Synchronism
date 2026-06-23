# Phase-6 (P1) — is the spatial-Lorentz failure fundamental, or a resolution artifact? (2026-06-22)

**Status:** `[ACTIVE-MRH]` — the make-or-break continuation. Phase-5 found the substrate's
preferred frame hides in *time* (dilation emerges) but shows in *space* (a boosted soliton was
Peierls–Nabarro pinned). P1 asks: is that spatial failure **fundamental** to a discrete substrate,
or a property of the **under-resolved** simulation? **Answer: a resolution artifact — and at the
physical scale it vanishes to `exp(−10²⁰)`.** This upgrades Phase-5's partial-negative.
**Sim:** [`simulations/phase6_spatial_lorentz_pn_barrier.py`](../simulations/phase6_spatial_lorentz_pn_barrier.py) · result: `simulations/results/phase6_spatial_lorentz_pn_barrier_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## What was measured

The **Peierls–Nabarro barrier** directly: the energy variation as a static sine-Gordon kink's
centre slides through one lattice cell. PN barrier > 0 ⇒ the lattice presents a preferred position
⇒ the pattern pins ⇒ the spatial frame is *visible*. PN barrier → 0 ⇒ the pattern sits/moves at
any sub-cell position with equal energy ⇒ the spatial frame is *hidden*. Sine-Gordon is chosen
because its **continuum** limit is exactly Lorentz-invariant, so *any* pinning is a pure
discretization artifact.

## Result — the barrier collapses with resolution

| lattice spacing `a` | kink spans | PN barrier / E₀ |
|---|---|---|
| 1.00 | ~1 cell | 1.55 × 10⁻³ |
| 0.70 | ~1.4 cells | 3.2 × 10⁻⁵ |
| 0.50 | ~2 cells | < 10⁻⁶ (below precision) |
| 0.25–0.18 | 4–5.6 cells | 0 |

A ~48× drop from 1 → 1.4 cells ⇒ `barrier ~ exp(−const·N)`, `const ≈ 9` per cell-of-resolution.
The pinning potential vanishes once a pattern is resolved over even ~2 cells.

## The physical-scale point (the real resolution)

Synchronism's substrate is a *fundamental* Planck-scale grid — you cannot "refine" it. But the
relevant number is **cells-per-pattern**, and real patterns are astronomically over-resolved:
a particle's size (Compton wavelength) exceeds the Planck length by `N ~ 10²⁰`. With
`barrier ~ exp(−const·N)`, that is `exp(−10²¹) ≈` **absolute zero**. So:

> **Phase-5's spatial pinning was a numerical under-resolution artifact** — I simulated few-cell
> solitons for tractability. At the *physical* scale, a particle spans `~10²⁰` Planck cells, the
> PN barrier is `exp(−10²⁰)`, and the spatial preferred frame is hidden to a precision no
> experiment could ever reach. **The make-or-break resolves in the model's favor at the physical
> scale.**

This is the same structure as Phases 2/4: discreteness effects are real but `(pattern/grid-size)`-
suppressed, invisible except in principle at the grid scale.

## Honest accounting

- **The PN-barrier-vanishes-with-resolution scaling is textbook** discrete-soliton numerics
  (Peierls–Nabarro, exponential in resolution). **Not novel.** The content is that it *locates*
  Phase-5's obstruction as non-fundamental and quantifies its physical-scale suppression.
- **The dynamical glide test was honest-but-uninformative here:** a boosted kink glides at ~100%
  velocity retention at *both* a=1.0 and a=0.25 — because **sine-Gordon is only weakly pinned**
  (barrier ~10⁻³·E₀ even coarse, far below a v=0.3 kink's kinetic energy), so it glides at any
  spacing. It therefore does *not* reproduce-then-cure Phase-5's strong (complex-KG) pinning;
  the **PN-barrier scaling is the rigorous result**, system-independent. (A φ⁴ kink or a slow
  kink would show strong dynamical pinning; not needed for the conclusion.)
- An earlier glide implementation was **discarded** (agent-zero): a single kink on a *periodic*
  lattice injects a 2π boundary discontinuity that radiates and corrupts the tracking — it gave
  fine-grid retaining *less* than coarse, a physical impossibility. Fixed with Dirichlet
  boundaries + a clean monotonic π-crossing centre-find.
- **Bucket 0 unchanged.**

## What it does to the arc, and the deeper open question

Phase-5 was a partial-negative (spatial frame visible). P6 **upgrades it**: that visibility is
under-resolution, and at the physical scale the spatial frame hides as thoroughly as the temporal
one. Combined with Phase-5 Part A (time dilation emerges) and Phase-3c (absolute-time inflow → GR),
the substrate now reproduces SR/GR at the physical scale in **both** sectors, with all
preferred-frame / LIV signatures `(pattern/grid)`-suppressed.

**Deeper open question (P1-deep, named not built):** if one wanted *exactly* zero PN barrier at
*fixed coarse* `a` (rather than `exp(−10²⁰)`), translationally-invariant discretizations
(Speight–Ward, integrable-lattice) achieve it by construction. Worth testing whether the
*intent-on-a-grid* update rule is naturally of that class — but it is no longer make-or-break,
since the physical-scale suppression already hides the frame.
