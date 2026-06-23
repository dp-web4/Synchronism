# Phase-7 (P2) — does a transport DOF derive gravity, or is the profile still fit? (2026-06-22)

**Status:** `[ACTIVE-MRH]` — follows the temporal/spatial research clue (the spatial sector is the
wall) into the unpaid part of Phase-3c (the GP inflow was *imposed*). **Result: the transport DOF
*derives* long-range gravity — the Yukawa obstacle dissolves from the substrate's own transport —
but it does *not force* the GR-matching profile. "Gravity is fit" shrinks from "fit the whole
inflow" to "fit one equation-of-state." Progress; loan still unpaid.**
**Sim:** [`simulations/phase7_transport_inflow_profile_forced_or_fit.py`](../simulations/phase7_transport_inflow_profile_forced_or_fit.py) · result: `simulations/results/phase7_transport_inflow_profile_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## The setup

Phase-3b: the complex-KG substrate's *massive* field can't host long-range gravity (Yukawa,
short-range). Phase-3c's fix — gravity as substrate **inflow** — reproduced full GR but *imposed*
the GP profile `u=√(2GM/r)`. P2 adds the transport degree of freedom dp specified (intent has
momentum / "where it wants to go", conserved under continuity) and asks: does a coherent high-Intent
core (a *sink* toward which intent reconstructs — a "mass") **derive** a long-range inflow, and do
the rules **force** the GP profile?

**Derivation (steady radial inflow, intent conserved):** continuity gives `ρ·v_r·r² = const`. So:
- **Incompressible** (`ρ=const`, the minimal/default EoS): `v_r ∝ 1/r²` — a power law, **long-range**.
- **GP / GR-matching:** `v_r ∝ 1/√r` requires a *specific* EoS, `ρ ∝ r^(−3/2)`.

Either way the inflow is a **power law, not Yukawa** — so the transport DOF dissolves the
short-range obstacle *for any EoS*, derived rather than imposed. That is the win.

## The test and result

Trace light (eikonal swimmer `H=c|k|+k·u`, from Phase-3c) through each inflow and measure the
deflection-vs-impact-parameter scaling `Δθ ∝ b^(−q)`. GR (observed) needs `q=1`.

| inflow profile | `Δθ ∝ b^(−q)` | matches GR (`q=1`)? |
|---|---|---|
| GP `v ∝ 1/√r` | **q = 1.02** | ✅ (= `4GM/c²b`, Phase-3c) |
| incompressible `v ∝ 1/r²` (minimal continuity) | **q = 4.04** | ❌ (`Δθ ∝ b⁻⁴` — falls far too fast) |

The minimal/default transport profile gives the *wrong* gravity: `b⁻⁴` deflection, negligible and
nothing like the observed `1/b`. The GR-matching profile exists but requires choosing the EoS
`ρ ∝ r^(−3/2)`.

## What it means — the freedom is located, not eliminated

**P2 answer: long-range is DERIVED; the GP profile is NOT forced.** The transport DOF is real
progress — it removes the Yukawa obstacle from the substrate's own rules, not by hand. But the
*profile* (hence whether light bends like GR) is set by the **equation-of-state** — how intent
density responds to a convergence — and minimal continuity (incompressible) gives the wrong answer
(`b⁻⁴`). Matching GR still requires *choosing* `ρ ∝ r^(−3/2)`.

So "gravity is fit" has **shrunk**: from Phase-3c's "fit the entire inflow profile" to "fit one
equation-of-state." That is genuine narrowing — but the loan is still unpaid. The open question is
now sharp and *spatial-emergence* in character (the research clue holds): **what, in the substrate's
rules, forces `ρ ∝ r^(−3/2)`?** That is a statement about how the substrate's density responds to a
coherence sink — exactly the kind of spatial-structure question the whole arc keeps bottoming out
on.

## Honesty

Not novel physics — effective-medium / fluid-analog gravity and the continuity derivation are
textbook; Bucket 0 unchanged. The content is: (a) the transport DOF **derives** the long-range
property (a real upgrade over Phase-3c's imposition), and (b) it **fails to force** the GR profile,
**locating the remaining freedom precisely as the equation-of-state** rather than the whole inflow.
Caveats: planar eikonal (exact orbital plane), weak field, the sink/EoS are minimal choices, and the
incompressible deflection underflowed past `b=80` (so its `q=4.04` is fit from the two resolved
points — enough to show "far steeper than `1/b`", which is the only claim made).
