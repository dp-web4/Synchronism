# Phase-8 (P2-deep) — multi-faith on capacity rules: which derives GR? (2026-06-22)

**Status:** `[ACTIVE-MRH]` — applies dp's **multi-faith** method (hypothesize promising paths, test
with feedback, refine/abandon/redirect) to the question Phase-7 left: *which capacity rule produces
the GR density profile `ρ ∝ r^(−3/2)`?* **Result: the GR-selecting rule is located precisely and
uniquely — `ρ ∝ v³` (held intent scales as the cube of flow speed) — derived from continuity with
no gravitational force postulated. "Gravity is fit" has shrunk to a single number.**
**Sim:** [`simulations/phase8_capacity_rule_multifaith.py`](../simulations/phase8_capacity_rule_multifaith.py) · result: `simulations/results/phase8_capacity_rule_multifaith_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Method (dp 2026-06-22)

We don't know the capacity rule — we *discover* it. Hypothesize a few promising candidates, test
which lands on GR, refine/abandon/redirect. The capacity rule (how much intent a cell holds /
receives, gated on momentum) *is* the equation-of-state, so each candidate predicts a density
profile and a light deflection. Steady radial inflow with intent conserved: continuity
`ρ·v·r² = const` is the backbone; the capacity rule is the closure.

## Faith A — momentum-coupled capacity `ρ ∝ vⁿ`

"A cell's held intent scales with its flow speed as `ρ ∝ vⁿ`" (`n` = how strongly capacity is gated
on momentum). With continuity: `v ∝ r^(−2/(n+1))`, `ρ ∝ r^(−2n/(n+1))`. Sweep `n`, measure the
light-deflection scaling `Δθ ∝ b^(−q)` by eikonal ray-trace:

| `n` (momentum coupling) | density exp | deflection `q` | GR? |
|---|---|---|---|
| 0 (incompressible) | 0.00 | 4.00 | — (Phase-7's wrong `b⁻⁴`) |
| 1 | −1.00 | 2.00 | — |
| 2 | −1.33 | 1.34 | — |
| **3** | **−1.50** | **1.02** | ✅ **GR** |
| 4 | −1.60 | 0.83 | — |
| 5 | −1.67 | 0.71 | — |

**`n=3` uniquely lands on GR** — `ρ ∝ r^(−3/2)` *and* `Δθ ∝ 1/b` — and it does so **from continuity
alone, with no gravitational force postulated** (the inflow is driven by the sink consuming intent,
not by an assumed `1/r²` pull). The deflection `q` crosses 1 exactly at `n=3`; the others are
eliminated (too steep below, too shallow above).

**So the capacity rule GR requires is sharp and falsifiable: `ρ ∝ v³`.** Gravity-is-fit has narrowed
across the arc:
- Phase-3c: fit the *entire* inflow profile (imposed `√(2GM/r)`);
- Phase-7: fit *one equation-of-state* (`ρ ∝ r^(−3/2)`);
- Phase-8: **the momentum-coupling exponent is exactly `3`** — a single number.

## Faith B — saturating capacity (redirect, not abandon)

Clip the held density at a ceiling `ρ ≤ I_max` (a cell can only hold so much). This only affects the
**centre** (`r < r_sat`, where `ρ` would exceed `I_max`); the **tail stays `ρ ∝ r^(−3/2)`** (GR
preserved outside `r_sat`). So a saturating capacity gives a **cored centre + GR tail** — which is
*exactly the structure of galactic dark-matter halos* (the core–cusp profile). **Redirect:** keep
the saturating rule for the galactic frontier (P4 / B2 RAR), where it would alter rotation curves
**without** changing solar-system lensing. That is a genuine lead toward the loan-repayment target.

## Honest accounting

- **Not novel physics**, and **`n=3` is *located*, not *derived*.** Power-law EoS → power-law
  profile is textbook fluid-analog gravity; finding the exponent that matches GR is curve-matching
  until the cell's *actual receive-dynamics* are shown to **force** `n=3`. That derivation is the
  open question — but it is now **a single number**, which is a far sharper target than "what is the
  inflow profile." **Bucket 0 unchanged.**
- The multi-faith did its job: one path (`n=3`) is the GR target, one (saturation) redirects to the
  galactic frontier, the rest are eliminated — precisely the refine/abandon/redirect outcome.
- Caveats: steady-state, weak-field, planar eikonal; the `ρ ∝ vⁿ` family is one closure family (a
  physically-motivated one — momentum-gated capacity — but not the only conceivable one); the
  incompressible `b⁻⁴` underflows past `b=80` (exponent from resolved points).

## The next question, now a single number

**Why `n=3`?** What, in a cell's receive-dynamics ("how much more intent it can take, gated on its
momentum"), would make held density scale as the *cube* of flow speed? That is the derive-vs-fit
line for gravity, reduced from a whole profile to one exponent. If a natural cell rule forces `n=3`,
gravity is derived and the loan moves; if `n=3` has to be hand-set, gravity stays fit — but either
way the question is now crisp enough to actually answer.
