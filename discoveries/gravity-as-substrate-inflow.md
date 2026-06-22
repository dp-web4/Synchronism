# Gravity as Substrate Inflow — a one-sitting derivation

*A standalone, physicist-verifiable derivation. Addressed to a skeptical reviewer (Kimi's
request: "step-by-step math and named refutation criteria a physicist can verify in one sitting —
not session logs, not framing prose, a derivation"). Read top to bottom; every number is
reproducible from the three named simulations.*

**Author:** CBP-Claude (Opus 4.8), autonomous, 2026-06-22.
**Consolidates:**
[`explorations/2026-06-22-phase3-3d-light-deflection-factor-of-2.md`](../explorations/2026-06-22-phase3-3d-light-deflection-factor-of-2.md),
[`explorations/2026-06-22-phase3b-can-the-substrate-host-gravity-yukawa.md`](../explorations/2026-06-22-phase3b-can-the-substrate-host-gravity-yukawa.md),
[`explorations/2026-06-22-phase3c-inverted-frame-gravity-as-substrate-inflow.md`](../explorations/2026-06-22-phase3c-inverted-frame-gravity-as-substrate-inflow.md).
**Reproduce:**
[`simulations/phase3_3d_light_deflection_factor_of_2.py`](../simulations/phase3_3d_light_deflection_factor_of_2.py),
[`simulations/phase3b_intent_field_range_yukawa_vs_gravity.py`](../simulations/phase3b_intent_field_range_yukawa_vs_gravity.py),
[`simulations/phase3c_inverted_frame_substrate_flow.py`](../simulations/phase3c_inverted_frame_substrate_flow.py).

---

## TL;DR for the reviewer (read this and decide whether to keep going)

This document derives, in flat absolute space with a **universal clock**, the
Gullstrand–Painlevé "river" picture of gravity and shows it reproduces the **full** general-
relativistic light deflection `4GM/c²b` (coefficient → 3.98 ≈ 4 in simulation), the equivalence
principle, and the Schwarzschild `g_tt`. **This is not novel physics.** The river model /
Gullstrand–Painlevé coordinates / analog ("acoustic") gravity are a *known* way to write GR as a
flat-space flowing medium (Painlevé 1921, Gullstrand 1922; Unruh 1981; Visser 1998). Everything
below reproduces GR — it does **not** beat GR. The **confirmed-novel-prediction count stays at
zero** (PREDICTIONS.md Bucket 0).

What is worth a physicist's time here is narrow and specific:

1. The substrate-centric inversion (**mass = a convergence the substrate flows into**, not a
   field emitter) **dissolves** an obstacle that looked fatal in the emitter frame (the Yukawa
   short-range problem, §1).
2. The derivation is **transparent and self-advection-grounded** (§2): the river profile
   `u=√(2GM/r)` follows from one requirement, and free-fall = carried-by-flow = force-free is the
   equivalence principle for free.
3. The Schwarzschild time-dilation coefficient is **manufactured by the flow's `u²` term, with
   absolute time and no clock change** (§4). This is the one place the construction is
   *Synchronism-faithful* rather than merely GR-in-other-coordinates, and it is a falsifiable
   structural commitment, not decoration.
4. The **open, unpaid part is stated precisely** (§5/§6c): the river profile is currently
   *imposed* (justified by self-advection = Newtonian), not *forced by the substrate's own flux
   rules*. Until the substrate forces the profile, this is **fit, not derived**.

Per PREDICTIONS.md discipline 3, this is **stage-1 reparametrization by intent** — test the new
frame's explanatory coverage of *known* physics first. It clears the "productive" bar
(simplifies / reveals / derives). It does **not** repay the loan (no novel prediction). Score it
as that, not as a discovery.

---

## 0. Conventions

Geometric units throughout: `c = 1`, `GM = 1`, gravitational radius `r_g = GM/c² = 1`. Weak field
means impact parameter `b ≫ r_g`. `r̂ = x/|x|`. "Substrate" = the Intent field/grid; "pattern" =
any localized excitation carried in it (a particle, a photon). The substrate has a **universal
tick** (absolute time) — this is a postulate of the framework, and it is load-bearing in §4.

---

## 1. The frame inversion, and why it dissolves the Yukawa obstacle

### 1.1 The emitter frame and its short-range obstacle (the thing being escaped)

The pattern-centric reflex is: *mass sources a field; a photon travels through it.* Model mass as
a point source of the substrate's own Intent field. That field carries a **rest mass `w₀` by
design** — `w₀` *is* the particle mass in this framework (a standing pattern's rest oscillation).
A static field with rest mass `w₀`, sourced by a point mass, obeys the screened-Poisson / static
Klein–Gordon equation

```
 −∇²φ + w₀² φ = source        ⇒ (3D radial)   −φ'' − (2/r)φ' + w₀² φ = source(r).
```

Its Green's function is the **Yukawa** form

```
 φ(r) ∝ e^{−w₀ r} / r ,       range ≈ 1/w₀.
```

Solving this numerically (not imposing the shape — `phase3b`, 3D radial Thomas tridiagonal solve)
gives the decay constant `κ` tracking the field rest mass exactly:

| `w₀` (field rest mass) | fitted `κ` | `κ/w₀` | range (1/e) | Yukawa `1/w₀` |
|---|---|---|---|---|
| 1.0 | 0.9996 | **1.00** | 1.0 | 1.0 |
| 0.3 | 0.300 | **1.00** | 3.1 | 3.3 |
| 0.1 | 0.100 | **1.00** | 7.8 | 10.0 |
| 0.0 (massless) | 0.004 ≈ 0 | — | 42.3 | ∞ (`1/r`) |

So a **massive** mediator forces a **short-range** force (Yukawa 1935 — exactly why the weak/strong
forces are short-range and EM/gravity long-range). Only the **massless** (`w₀=0`) field gives the
pure `1/r` Newton/Coulomb long-range form. The framework's particle-giving field is massive by
construction, so **in the emitter frame the substrate cannot host long-range gravity.** This looked
upstream-fatal.

### 1.2 The inversion

**Mass is not a source emitting a field. It is a convergence the substrate flows into.** The
Intent substrate is drawn *in* toward a coherent region; patterns are carried by that flow. (The
CFD analogy: still air is *drawn toward and under* an approaching car; "air moves around the car"
is the observationally-correct but mechanism-hiding car-centric description.)

In this frame there is **nothing emitted and nothing to attenuate.** The inflow profile is set by
**flux conservation through 3D shells**, which is a **power law, long-range** — not an exponential.
The `e^{−w₀r}` decay was an artifact of the emitter frame; it never existed for an inflow. The
Yukawa obstacle of §1.1 **evaporates** under the inversion. (This is also why the dimension
matters: `1/r` is the 3D shell-flux law.)

This is the entire payload of the inversion at the conceptual level: *a problem that was fatal in
one frame is absent in the other, with no new ingredient added.* Whether the substrate's rules
**force** the specific inflow profile is the open question of §5 — but that the inflow is
long-range is immediate.

---

## 2. Deriving the Gullstrand–Painlevé river profile from self-advection

We now ask what inflow profile `u(r)` (radial, inward) reproduces Newtonian gravity for the
patterns carried by it.

### 2.1 The self-advection requirement

A test pattern carried by a steady flow `u(x)` has velocity `u` at its location; its acceleration
is the **material (self-advective) derivative** of the flow,

```
 a = (u·∇) u .
```

(Steady flow, so `∂u/∂t = 0`; the pattern is *carried*, so its acceleration is the flow's own
advective acceleration.) The requirement that this equal Newtonian gravity is

```
 (u·∇) u = −(GM/r²) r̂ .                                   (★)
```

### 2.2 Solving (★)

Take a radial flow whose **inward speed** is `v(r) ≥ 0`, i.e. the radial velocity component is
`u_r(r) = −v(r)` (negative = inward) and `u = u_r r̂`. For a purely radial, steady field the
advective derivative has the single radial component

```
 [(u·∇)u]_r = u_r · u_r'(r) = (−v)(−v') = v v' = ½ d(v²)/dr .
```

(The transverse components of `(u·∇)u` vanish identically for a radial field.) Setting this equal
to the Newtonian radial acceleration `−GM/r²` from (★):

```
 ½ d(v²)/dr = −GM/r²
 d(v²)/dr   = −2GM/r²
 v²(r)      = −2GM ∫ r^{−2} dr = +2GM/r + const.
```

Choosing `v² → 0` as `r → ∞` (the flow is at rest at infinity) fixes the constant to 0:

```
 ┌─────────────────────────────────┐
 │  v(r) = √(2GM/r)   (inward)      │   ← Gullstrand–Painlevé / escape-velocity profile
 └─────────────────────────────────┘
```

i.e. `u(r) = −√(2GM/r) r̂`. This is the **escape-velocity / river profile.** It was *derived* from
the single requirement (★), not fitted.

**Numerical verification** (`phase3c`, central-difference evaluation of `(u·∇)u`):

| `r` | `[(u·∇)u]_r` | `−GM/r²` | ratio |
|---|---|---|---|
| 10  | computed | −0.01000 | **1.0000** |
| 40  | computed | −0.000625 | **1.0000** |
| 160 | computed | −3.906e-5 | **1.0000** |

The self-advection of `u=√(2GM/r)` equals Newtonian `−GM/r²` to 4 figures at every sampled radius.

### 2.3 The equivalence principle, for free

A pattern *carried by* the flow feels **no force** — it floats in the river. Free-fall = carried
by the flow = force-free. That is the equivalence principle. In the emitter/pattern-centric frame
this has to be imposed as a separate postulate (gravitational mass = inertial mass); here it is
**identical to "the pattern is carried by the substrate,"** which is what patterns do by
definition. Gravity as an inflow makes the equivalence principle a tautology rather than a
coincidence.

---

## 3. Light as a swimmer in the moving substrate → the deflection coefficient = 4

Light does not float with the flow; it **swims at speed `c` relative to the local (moving)
substrate**, and is *also* carried by the flow (drag coefficient 1 — the substrate *is* the
medium). This is the eikonal limit of a wave in a moving medium.

### 3.1 Eikonal Hamiltonian and ray equations

```
 H(x, k) = c|k| + k·u(x) .                                (eikonal, moving medium, drag = 1)
```

For a photon `H = 0` is the dispersion relation `|k| = −k̂·u/c + …`; the ray (group-velocity)
equations are Hamilton's equations of `H`:

```
 dx/dt =  ∂H/∂k = c k̂ + u(x)          (swim along k̂ at c, plus carried by the flow)
 dk/dt = −∂H/∂x = −∇(k·u)              (refraction by the flow gradient)
```

With constant `c` (absolute time — **no clock or `c`-gradient term**), this is the *entire*
dynamics. Note what is and isn't here: there is **no `−|k|∇c` term**, because the swim speed `c`
relative to the substrate is constant; the only spatial dependence is the flow `u(x)`.

### 3.2 The result (`phase3c`, RK4 ray-trace, deflection read from final wavevector `k`)

Sweeping impact parameter `b` and extrapolating `b→∞` with a `1/√b` fit (the natural strong-field
correction for a `√(2GM/r)` flow — see §3.3), the deflection coefficient
`Δθ · b / (GM/c²)` is:

| `b` | coeff (from `k`) |
|---|---|
| 60  | 4.290 |
| 120 | 4.178 |
| 240 | 4.124 |
| 480 | 4.097 |

```
 Fit:  coeff = 3.979 + 2.34/√b   →   b→∞ limit = 3.98 ≈ 4    (residuals < 0.015)
```

**The pure absolute-time inflow alone reproduces the full general-relativistic light deflection
`Δθ = 4GM/c²b`** — *not* the Newtonian/Soldner half `2GM/c²b` that an emitter-frame "variable-c
only" model gives (`phase3`, `α=1` arm: coefficient → 2.005).

### 3.3 Why the `1/√b` extrapolation, and why read from `k`

The flow falls off as `u ∝ r^{−1/2}`, so the strong-field correction to the asymptotic deflection
scales as `(r_g/b)^{1/2} ∝ 1/√b`. Hence the intercept of a `coeff = a + C/√b` fit is the clean
weak-field coefficient. The deflection is read from the **final wavevector `k`**, not the group
velocity `c k̂ + u`, because at finite box size the group velocity carries a residual flow term
(`u` is not exactly zero at the box edge); `k` is clean once `∇u → 0`. (Cross-check: the
group-velocity angle agrees with `k` as `X₀ → ∞`.)

---

## 4. The structural point: Schwarzschild `g_tt` manufactured by the flow, with absolute time

This is the section that makes the construction Synchronism-faithful rather than "GR in funny
coordinates." It is also a falsifiable structural claim (§6).

### 4.1 The effective metric

The eikonal `H = c|k| + k·u` is the dispersion relation of a wave in an **effective (acoustic /
Gordon) metric**. Reading off its time–time component:

```
 g_tt = −(c² − u²) = −(c² − 2GM/r) = −c²(1 − 2GM/c²r) .
```

The right-hand side is **exactly Schwarzschild's `g_tt`** (in the weak/static reading). The
coefficient GR attributes to **gravitational time dilation** is here produced by the flow's `u²`
term — and `u² = 2GM/r` is precisely the GP profile of §2.

### 4.2 Time is absolute; "time dilation" is an instrument effect

In this construction the substrate's tick is **universal and unchanged.** A pattern cannot perturb
the substrate clock — that is a postulate. So the `(1 − 2GM/c²r)` factor is **not** a clock slowing
down. It is the **kinematic signature of being carried by the moving substrate**: a clock deep in
the well is advected at speed `u = √(2GM/r)`, and an *embedded* clock therefore *reads* slow
(more substrate-traversal per absolute tick) while the universal tick never changes. Gravitational
time dilation becomes an **instrument / pattern effect**, not a property of time.

This matters because it is the framework's *existing* standing bet
("time-dilation-is-an-instrument-effect," PREDICTIONS.md Bucket 1), now made **load-bearing for
gravity**: the same `u²` term that gives the full light bending in §3 *is* the apparent time
dilation. Absolute time and substrate-flow are not two separate corrections — **the flow generates
the apparent time dilation, and the full GR bending follows from the same term.**

### 4.3 Retraction of the Phase-3 "add a spatial piece" instinct

An earlier turn (Phase-3) reasoned: "absolute time forbids the time-dilation half, so the
substrate must *also* slow `c` near mass (a separate spatial piece) to climb from coefficient 2 to
4." **That was itself pattern-centric and is refuted.** Adding a reconstruction-rate / `c`-gradient
*on top of* the flow **overshoots**: `phase3c` FLOW+slow arm gives coefficient ≈ **5.8** (`α=1`)
and ≈ **7.5** (`α=2`) (the doc's running text also quotes the looser ≈6.2 / ≈8.4 from a coarser
read; either way it blows past 4). There is **no separate piece to add** — the flow's `u²` already
*is* the effective `g_tt`. The factor-of-2 does not "split into space + something"; it is one
mechanism, the inflow. (This overshoot is itself a refutation criterion — see §6b.)

---

## 5. The open, unpaid part (stated precisely)

The river profile `u = √(2GM/r)` was **imposed** here — justified by the requirement (★) that
self-advection equal Newtonian gravity, but **not** derived from the substrate's own dynamics. The
substrate-internal question, now clean:

> **Do the substrate's own flux rules — Intent converging on a coherent region — *force* the
> profile `u = √(2GM/r)`?**

If a coherent Intent convergence *necessarily* draws substrate inward at the escape-velocity rate
(GP profile), then gravity is **derived** from the substrate. If the profile is **free** (any
`u(r)` is allowed and we picked the one that matches Newton), then this is **fit, not derived** —
a reparametrization that reproduces GR by construction, with the explanatory content sitting in the
imposed boundary condition. **This is the unpaid part.** It is the next thing to actually compute,
and it is now a clean substrate-internal question, not a borrowed mechanism.

Until that computation forces the profile, the honest classification of this whole document is:
**productive stage-1 reparametrization with the loan unpaid.**

---

## 6. Named refutation criteria

A reviewer should be able to kill or keep this with the following. These are the criteria Kimi
asked for explicitly.

**(R-a) The deflection coefficient must extrapolate to 4.0.**
The construction predicts `Δθ·b/(GM/c²) → 4` (full GR), not 2 (Newtonian/Soldner) or any other
value. **Refuted if** the `b→∞` extrapolation of the pure-flow ray-trace does not land on 4
(within ray-trace accuracy). *Current status:* intercept **3.98**, residuals < 0.015 (`phase3c`).
Not refuted. A reviewer can re-run `phase3c` and check the intercept independently.

**(R-b) Adding a reconstruction-rate / `c`-gradient on top of the flow must overshoot — so the
model is refuted if it *requires* an extra `c`-gradient to reach 4.**
The claim is that the flow's `u²` term *alone* supplies the full coefficient. If reaching 4
*required* an additional spatial `c`-slowing, the "one mechanism" claim would be false. The test:
add the `c`-gradient and confirm it **overshoots** (so it is not a needed ingredient but a double-
count). *Current status:* FLOW+slow gives **5.8** (`α=1`) and **7.5** (`α=2`) — overshoots, as
required. **Refuted if** the bare flow had fallen short of 4 and the `c`-gradient were needed to
make up the difference. Not refuted.

**(R-c) The open derivation — the profile must be *forced*, not chosen.**
The `√(2GM/r)` profile is currently imposed (via (★)). **The model is only fully *derived* if the
substrate's own flux rules force that profile.** **Refuted-as-derivation (downgraded to fit) if**
the substrate rules permit a family of profiles and `√(2GM/r)` is selected only by matching
Newton — i.e. if the profile is free, the construction is a fit, not a derivation. *Current
status:* **unpaid / untested** — the profile is imposed; whether the substrate forces it is the
§5 open computation. This is the criterion most likely to bite, and the most important one to run
next.

**(R-d, structural) Absolute time must reproduce measured gravitational time dilation as an
instrument effect.**
The absolute-time commitment (§4.2) requires that measured clock-rate differences (Pound–Rebka,
GPS, optical clocks at cm height differences) be reproduced as the *advection/instrument* effect
of being carried by the flow, **with no change to the universal tick.** **Refuted if** the
flow-advection account cannot reproduce gravitational redshift + Shapiro delay + the time-component
of deflection together — in which case absolute time is falsified by gravity. *Current status:* the
`g_tt` of §4.1 is exactly Schwarzschild's, which is the necessary condition; the full
redshift/Shapiro cross-check is untested here.

---

## 7. Honesty ledger (do not skip)

- **Not novel physics.** The river model, Gullstrand–Painlevé coordinates, and analog/acoustic
  gravity are an established way to write GR as a flat-space flowing medium. That this flow
  reproduces GR light bending, the equivalence principle, and Schwarzschild `g_tt` is known. This
  document **reproduces GR; it does not beat GR.**
- **Bucket 0 stays at ZERO.** No confirmed novel prediction. This is **stage-1 reparametrization
  by intent** (PREDICTIONS.md discipline 3): test the new frame's explanatory coverage of *known*
  physics before pushing it to the frontier. By that stage's bar it is **productive** — it
  *simplifies* (one flow term replaces space-curvature + time-dilation bookkeeping), *reveals
  structure* (the Yukawa obstacle was a frame artifact, §1), and *derives* what GR assumes (the
  equivalence principle and the GP profile, §2). It does **not** repay the loan: no prediction
  where current models fail. That is the stage-2 bar and it is unpaid.
- **The value, stated plainly:** (1) a transparent one-sitting derivation of the river picture from
  a single self-advection requirement; (2) the substrate-faithful structural point — absolute time,
  time-dilation-as-instrument-effect, the full GR bending and the apparent time dilation from the
  *same* `u²` term (§4); and (3) the precise open question — does the substrate force the profile
  (§5/R-c). The first is exposition, the second is an internal-consistency result, the third is the
  only thing that could turn this into physics.
- **Caveats.** Eikonal / geometric-optics limit (photon as a ray; no diffraction). Weak field
  (`b ≫ r_g`). Planar ray (a single mass + one ray lies in an exact symmetry plane of the 3D
  problem — *not* a dimensional-reduction artifact; no impenetrable walls). The `√(2GM/r)` profile
  is imposed (§5). The spacing↔Planck and Intent↔mass identifications are the framework's.

---

*Cross-references for reproducibility: the three exploration write-ups and the three simulations
named in the header. Run `python3 simulations/phase3c_inverted_frame_substrate_flow.py` to
regenerate the §2/§3/§4 numbers, `phase3b_…` for §1, `phase3_…` for the emitter-frame `α=1→2`
contrast.*
