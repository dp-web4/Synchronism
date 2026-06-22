# Phase-3b — can the Intent substrate host long-range gravity at all? (2026-06-22)

**Status:** `[ACTIVE-MRH]` — a **productive negative** that sits *upstream* of the factor-of-2:
before asking whether the substrate gets gravity's 2× right, ask whether it can produce a
long-range `1/r` field at all. It can't — not from its particle sector.
**Sim:** [`simulations/phase3b_intent_field_range_yukawa_vs_gravity.py`](../simulations/phase3b_intent_field_range_yukawa_vs_gravity.py) · result: `simulations/results/phase3b_intent_field_range_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## A discarded near-circular first attempt (kept visible on purpose)

The first version of this probe imposed a `sech`-shaped core in a massive field and "found" the
induced index well was short-range — but that just echoed the shape I typed in, and its FDTD
cross-check silently failed (returned `None`). **Discarded** rather than reported: a measurement
whose dummy didn't run is not a measurement (agent-zero). The honest question is non-circular and
sharper — *solve* for the field instead of imposing it.

## The question, stated so it can't be circular

Phase-3 imposed `n(r)=1+α/r` by hand and watched rays bend — that tests geometry, not the
framework. Here: **a static field sourced by a point mass has a spatial range set by the field's
own rest mass.** Don't impose the shape; solve the genuinely-3D radial static field equation
(`−∇²φ + w₀²φ = source`, the static Klein-Gordon / screened-Poisson) and read off the range.

Faithful to absolute time (dp's correction): this is a static **spatial** field. No tick, no
clock varies.

## Result (3D radial solve, Thomas tridiagonal)

| w₀ (field rest mass) | fitted κ | κ/w₀ | range (1/e) | expected Yukawa range 1/w₀ |
|---|---|---|---|---|
| 1.0 | 0.9996 | **1.00** | 1.0 | 1.0 |
| 0.3 | 0.300 | **1.00** | 3.1 | 3.3 |
| 0.1 | 0.100 | **1.00** | 7.8 | 10.0 |
| 0.0 (massless) | 0.004 ≈ 0 | — | 42.3 | ∞ (long-range `1/r`) |

The static field is **Yukawa** `e^{−w₀r}/r` whenever the field is massive: the decay constant κ
*equals* the field rest-mass w₀ exactly, so the range is `1/w₀` — **short**. Only the **massless**
field (w₀=0) gives κ≈0, the pure `1/r` Coulomb/Newton form — **long-range**.

## What it means — the real gap, upstream of the factor-of-2

The framework's Intent field carries a **rest mass w₀ by design** — that w₀ *is* the particle's
mass (the entire point of Phase-1/2: a standing pattern's rest oscillation = mass). A massive
mediator forces a **short-range** force — this is Yukawa's 1935 result, and it's exactly why the
weak and strong forces are short-range while EM and gravity are long-range.

So the substrate's **particle sector cannot be the source of gravity.** Its field dies off as
`e^{−w₀r}` over a distance `1/w₀` (a particle's own Compton wavelength) — it cannot reach across
a solar system. And `1/r` exists only in **3D** — the dimension mattered here too.

This sits **upstream** of everything in Phase-3: there is no point asking "does the substrate get
gravity's factor-of-2 from space alone, or is it falsified" until a `1/r` field exists to bend
light *at range*. The factor-of-2 analysis assumed the long-range field; this shows the framework
doesn't yet produce one.

## The honest open problem (sharper than before)

To host gravity, Synchronism needs a **separate massless Intent mode** — a graviton-analog, or a
collective long-wavelength excitation of the field, or a genuine `1/r` far-tail — that the current
specification does **not** provide. Three candidate directions, none yet modeled:
1. **A massless sector**: a second field component with w₀=0 (like the photon is the massless EM
   mode). What in the Intent ontology would be massless, and why?
2. **A collective/entropic effect**: gravity as a long-wavelength, many-pattern statistical
   effect (à la entropic-gravity proposals) rather than a fundamental mediator — emerging from the
   substrate's coherence dynamics, not a single field mode.
3. **The grid itself as the mediator**: if mass perturbs the *lattice geometry* (spacing) and that
   perturbation is long-range, gravity is a property of the grid, not a field on it — which would
   reconnect to the absolute-time "the 2× must be spatial/structural" requirement from Phase-3.

## Honesty

Not novel physics (Yukawa range, massless⇒long-range are textbook). The value: a **non-circular,
derived** demonstration that the framework's particle-giving massive field is the *wrong kind of
object* to be gravity, which **relocates the open problem** from "get the factor-of-2 right" to
"produce a long-range field at all." That's a more honest and more fundamental statement of where
the physics stands. Bucket 0 unchanged; this sharpens Bucket 1's gravity line into a stated
open problem with three candidate attacks.
