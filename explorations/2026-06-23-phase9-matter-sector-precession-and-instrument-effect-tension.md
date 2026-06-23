# Phase-9 — the MATTER sector: precession is reproduced *exactly*, and that creates a named tension (2026-06-23)

**Status:** `[ACTIVE-MRH]` — closes the sector the arc never tested (timelike orbits) and, in doing
so, surfaces a foundational tension the arc had been treating as a *win*. **Result: the
inflow/swimmer model reproduces full Schwarzschild perihelion precession exactly (ratio
1.000000, including strong-field corrections), because the matter swimmer Hamiltonian is
*algebraically identical* to the Gullstrand–Painlevé matter dispersion. The same fact that makes
this work — universal advection of every pattern by one flow `u(r)` — is the equivalence
principle, and the equivalence principle logically *forbids* the mechanism-dependence that the
standing "time-dilation-is-an-instrument-effect" bet needs as its discriminator.**
**Sim:** [`simulations/phase9_matter_sector_precession.py`](../simulations/phase9_matter_sector_precession.py) · result: `simulations/results/phase9_matter_sector_precession_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## The gap this closes

Phase-3c / 7 / 8 established "gravity as substrate inflow reproduces GR" — but **every one of
those phases traced LIGHT** (the eikonal swimmer `H = c|k| + k·u`, deflection `4GM/c²b`). None
ever put **matter** in an orbit. So "reproduces GR" was, strictly, "reproduces GR's **null**
sector." The other canonical GR test — **perihelion precession**, a *timelike*-sector effect — was
never computed. Light deflection and precession are independent observables: a refractive/inflow
model can get light right by construction and still miss precession (this is the classic
**analog-gravity limitation** — analog metrics reproduce kinematics but not always the full
geodesic structure). So the matter sector was a genuine, untested place the arc could have failed.

## The faithful matter rule (not a new postulate)

The framework's light rule is Galilean advection of the medium-frame dispersion in absolute time:
`H = c|k| + k·u`. The faithful **matter** version is the same rule with the massive medium-frame
dispersion:

```
H = sqrt(m²c⁴ + c²|p|²) + p·u          (m→0 recovers the light rule exactly)
```

This is not an extra assumption — it is the *only* mass-carrying generalization consistent with the
arc's own "patterns advected by the flow, in absolute time" ontology.

## Result 1 — precession is reproduced EXACTLY (and analytically)

Squaring the swimmer dispersion (the `(p·u)²` term regenerates the `g_rr` correction the linear
form looks like it drops):

```
(H − p·u)² = m²c⁴ + c²|p|²
 ⇒  H² + 2H(p·u) = m²c⁴ + (1 − u²/c²)p_r² + p_φ²/r²
```

This is **algebraically identical** to the exact Gullstrand–Painlevé (= Schwarzschild) matter
dispersion. Sim check over 10⁵ random phase-space points: `max|LHS − RHS| = 5.7×10⁻¹⁴`.

Consequently the turning-point condition `ṙ=0` of the swimmer reduces *exactly* to Schwarzschild's
`E² = (1−2GM/r)(1+L²/r²)` — and the measured precession matches term-for-term:

| r_peri | r_apo | a | e | swimmer Δφ/orbit | Schwarzschild Δφ/orbit | swim/Schw |
|---|---|---|---|---|---|---|
| 120 | 180 | 150 | 0.20 | 0.135148 | 0.135148 | **1.000000** |
| 120 | 240 | 180 | 0.33 | 0.121253 | 0.121253 | **1.000000** |
| 200 | 300 | 250 | 0.20 | 0.080047 | 0.080047 | **1.000000** |
| 60 | 120 | 90 | 0.33 | 0.249859 | 0.249859 | **1.000000** |
| 40 | 90 | 65 | 0.385 | 0.371053 | 0.371053 | **1.000000** |

The strong-field rows (e.g. r_peri=40) **exceed** the weak-field formula `6πGM/c²a(1−e²)` by ~9% —
and the swimmer tracks the *full* Schwarzschild value, not the leading-order one. So this is not
"matches GR to first order"; it is **the same orbit**.

**What it means for the loan/Bucket-0 map.** Phase-3c (light) + Phase-9 (matter) together: wherever
`u = √(2GM/r)`, the inflow model is **Gullstrand–Painlevé in disguise = Schwarzschild, identically,
for both sectors**. There is therefore **zero wiggle room** against GR in any classical-GR regime.
A novel (Bucket-0) prediction can *only* live in one of two places:
1. where the inflow **profile departs from GP** — i.e. where the capacity/EoS rule changes `u(r)`
   (Phase-8 Faith-B: saturating capacity → cored centre + GR tail, the galactic frontier); or
2. where **discreteness breaks the continuum identity** — sampling/Umklapp/LIV (bets B7, B-dispersion).
Everywhere else the agreement is *forced*, not fit. This is a sharper statement of the frontier
than "reproduce GR, then pay out": it says exactly which two doors the payout must come through.

## Result 2 — the named foundational tension (the uncomfortable part)

Phase-3c reads its `g_tt = −(c²−u²)` result as **vindicating** the framework's standing bet that
*gravitational time dilation is an instrument effect, not spacetime* (PREDICTIONS Bucket-1 / SPINE
"invitation" test 2 / the pendulum-in-a-centrifuge analogy). Phase-9 shows that reading is
**self-undercutting**:

- The **discriminator** the instrument-effect bet names is **mechanism-dependence**: "spacetime
  says *every* clock dilates identically; an instrument effect allows physically different clock
  mechanisms to respond differently" (SPINE). The pendulum analogy's entire force is that a
  pendulum and an atomic clock in the same centrifuge diverge *differently* — the effect is on the
  instrument, so it is mechanism-specific.
- But the inflow model advects **every pattern by the same universal `u(r)`** — there is **no
  per-mechanism parameter**. The static-clock rate is `dτ/dt = √(1−2GM/r)` for *any* internal
  mechanism. The model predicts **clock universality**, which is *exactly* what "spacetime" predicts
  and exactly what the instrument-effect discriminator needs to **fail**.
- Worse, universality here is not incidental: it **is** the equivalence principle, which Phase-3c
  separately *celebrates deriving* ("free-fall = floating in the river"). **Deriving the equivalence
  principle and predicting mechanism-dependence are mutually exclusive.** The same single mechanism
  cannot both couple universally (EP, no mechanism-dependence) and couple per-mechanism (the bet's
  discriminator).

So two things the project currently lists as **strengths** cannot both pay out:
- (Phase-3c) the inflow model derives the equivalence principle / universal advection; and
- (Bucket-1 / SPINE) time-dilation-as-instrument-effect, distinguished from GR by *mechanism-dependence*.

The inflow result does **not** vindicate the instrument-effect bet — it **removes the bet's only
handle**. The "instrument effect" survives as an *interpretation* (the clock is a pattern carried by
the flow), but as a *distinguishing physics prediction* it is dead in the continuum model: it makes
predictions identical to GR, by construction.

### The constructive repair (don't just kill it — relocate it)

The bet is only dead *in the continuum, universal-`u` limit*. Its one survival path is the **same
door as Bucket-0 door #2**: if **discreteness** makes different patterns sample the grid slightly
differently, universality is broken at the LIV/sampling order, and *that* deviation would be
genuinely mechanism-dependent and testable against the (already very tight) clock-universality
bounds. So the honest move is: **relocate the instrument-effect bet from "continuum inflow" (where
it is provably indistinguishable from GR) to "discreteness breaks universality" (where it has real,
falsifiable content)** — and state plainly that the continuum inflow result *refutes* the bet's
distinguishing power rather than confirming it.

## Honesty / accounting

- **Not novel physics.** GP = Schwarzschild is textbook; that the swimmer equals GP for matter is
  internal-consistency, not a new result. **Bucket 0 unchanged (0).**
- **What is new here** is (a) the matter sector was actually *computed* for the first time and
  reproduces full Schwarzschild precession exactly — extending Phase-3c from null to timelike and
  confirming the inflow model is the *complete* GP chart, not just its light cone; (b) this maps the
  Bucket-0 frontier to exactly two doors (profile-departure, discreteness); and (c) it surfaces the
  EP ⊥ instrument-effect tension, correcting Phase-3c's "vindication" reading.
- **Caveats:** weak-to-moderate field (strong-field rows up to e=0.385, r_peri=40GM still well
  outside the horizon); planar; static central sink; test particle (geodesic, no self-force). The
  dispersion identity (Result 1) is exact and field-strength-independent, so the precession match is
  not a small-parameter coincidence.

## So what

The arc's celebrated gravity result is now *complete* (both sectors) — and that completeness is
double-edged. It confirms the inverted frame lands on GR with no error, **and** it proves the frame
has *no room* to differ from GR anywhere `u`=GP, which forces the entire Bucket-0 hope onto two
narrow doors. And it converts a supposed win (instrument-effect "vindicated") into a named tension
(EP-derivation forbids the bet's discriminator). That is the productive-discomfort outcome: the
framework absorbed the *light* challenge by being GR-in-disguise; the *matter* challenge shows that
"in disguise" is total — which is exactly why the disguise cannot, by itself, ever pay the loan.
