# Phase-11 (door #1, galactic) — flat curves are OUTSIDE the capacity family; door #1 re-enters the MOND cage (2026-06-23)

**Status:** `[ACTIVE-MRH]` — sharpens Phase-10's "needs a transition" into a no-go. **Result: a
flat rotation curve (isothermal `ρ∝r^(−2)`, `α=2`) is the *excluded endpoint* of the entire
Phase-8 capacity family `ρ∝vⁿ` — `α = 2n/(n+1) < 2` for every finite `n`, reached only as
`n→∞` — and self-gravity does not move it. Reaching flat curves therefore requires *leaving* the
kinetic inflow family for a pressure-supported (static) mechanism AND an acceleration-keyed `a0`
switch between regimes. That switch is non-local (Milgrom; S689) and is the MOND interpolation
(refuted at γ=2, S661; reparametrization, Bucket 3). Of Phase-9's two Bucket-0 doors, door #1
(galactic/profile) now leads only back to the known cage; door #2 (discreteness/LIV) is the sole
remaining place a novel prediction can live.**
**Sim:** [`simulations/phase11_capacity_family_flatcurve_nogo.py`](../simulations/phase11_capacity_family_flatcurve_nogo.py) · result: `simulations/results/phase11_capacity_family_flatcurve_nogo_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## The bound (the new content)

Phase-10 found "n=3 gives a rising `r^{1/4}` curve, not flat, so a transition is needed." The
stronger statement: **no member of the family gives flat.** Steady spherical inflow conserving
intent (continuity `ρ·v·r² = const`) plus the capacity closure `ρ∝vⁿ` gives, algebraically,

```
v ∝ r^(−2/(n+1)),   ρ ∝ r^(−α)   with   α = 2n/(n+1)
```

| n | α = 2n/(n+1) | v_c slope = 1−α/2 | curve |
|---|---|---|---|
| 1 | 1.000 | 0.500 | rising |
| 3 (GR, Phase-8) | 1.500 | 0.250 | rising `r^{1/4}` |
| 10 | 1.818 | 0.091 | rising |
| 50 | 1.961 | 0.020 | rising |
| ∞ | 2.000 | 0.000 | **flat** (limit only) |

`α=2` requires `2n/(n+1)=2 ⇒ 0=2`: **no finite n.** Flat curves are the family's *excluded
endpoint*. And the self-consistent self-gravitating solve (integrate Poisson on the
continuity-fixed profile, measure the realized `v_c` slope) returns `1−α/2` to machine precision
for `n=1,3,10,100` — **self-gravity does not change the exponent**, because continuity + capacity
fix it before gravity is consulted.

## Why this is anchored, not a modeling choice

The continuity relation `ρvr²=const` is not an assumption that can be relaxed to escape the bound:
in steady state it **is** intent conservation (equal flux through every shell), and intent
conservation is a *non-negotiable axiom* (FUNDAMENTALS; SESSION_PRIMER's "do not violate"). To get
`α=2` from an *inflow* you would have to create or destroy intent in the halo — forbidden. So
within the conservation axiom the kinetic inflow family **cannot** produce flat curves, full stop.

## What α=2 actually is (a different mechanism)

The isothermal `ρ=σ²/(2πG r²)` that gives `v_c=√2·σ=const` is a **hydrostatic, pressure-supported**
solution: it solves `dP/dr = −GMρ/r²` with `P=σ²ρ` (residual `2.9e-16` in the sim), with the
intent essentially **static** (`v→0`), not flowing. That is a *second mechanism* bolted beside the
inflow one — the same "two distinct substrates connected only by narrative" pathology PREDICTIONS.md
already flags for the particle sector (R(I) diffusion vs wave+focusing). It is not a different
exponent of one rule; it is a different rule.

## The cage

So door #1 demands **both**:
1. a regime change from kinetic inflow (`α=3/2`, GR near matter) to pressure-supported isothermal
   (`α=2`, flat far away) — two mechanisms, not one rule; and
2. an **acceleration-keyed switch at `a0`** deciding which regime holds where.

Requirement (2) is a *non-local, acceleration-keyed* law — exactly the Milgrom-2005 non-locality the
repo already identified (S689) and exactly the MOND interpolation function. When `a0`/shape are fit,
it is MOND (Bucket 3, reparametrization); when pinned at the framework's γ=2 it is **refuted**
(ΔBIC=+184 on SPARC, S661). Either way it is not novel. **Door #1, walked all the way through,
deposits you back in the MOND cage.**

## Consequence for the Bucket-0 map (the so-what)

Phase-9 reduced the entire novel-prediction frontier to two doors: #1 profile-departure (galactic),
#2 discreteness (LIV/Umklapp). Phase-11 shows **door #1 leads only back to MOND** — refuted or
reparametrizing, no Bucket-0 exit. That leaves **door #2 (discreteness) as the sole remaining place
a genuine novel prediction could live.** This is a productive elimination: it tells the program to
stop mining the galactic/profile frontier for novelty (it is MOND in the framework's vocabulary) and
to concentrate the loan-repayment effort on the discreteness channel — bets B7 (vacuum Umklapp) and
the Planck-scale LIV dispersion — where the continuum identity that locks the framework to known
physics can actually break.

## Honesty

Not novel physics; **Bucket 0 unchanged (0).** The content is a no-go: (a) flat curves are outside
the capacity family `ρ∝vⁿ` (α<2 ∀ finite n), anchored to intent conservation, not a modeling
artifact; (b) self-gravity doesn't rescue it; (c) the isothermal target is a separate
pressure-supported mechanism + a non-local `a0` switch = MOND. Caveats: steady-state, spherical,
power-law capacity closure; real galaxies are disks and time-dependent (this changes coefficients,
not the α<2 cap, which is a continuity+closure consequence). A non-power-law closure could reach
α=2 in principle — but any closure giving GR (α=3/2) near matter and isothermal (α=2) far out *is*
the acceleration-keyed MOND transition, so the escape is the cage.

## Connection to the same-day "untestable with foreseeable data" proposal

`Research/proposals/framework_untestable_foreseeable_data.md` (filed 2026-06-23, maintainer track)
makes the *empirical* version of this finding: every test in the catalog (TEST-01…21 — all
galactic/RAR/wide-binary/cluster, i.e. **door-#1-type**) sits below the relevant systematics floor,
so it proposes upgrading "0 confirmed" to "untestable with foreseeable data." Phase-11 supplies the
**structural reason** that empirical survey keeps finding nothing: those tests live in the
continuum/profile sector, which is GR (Phase-9) ⊕ MOND (Phase-11) *by exact construction* — there is
no C(ρ)-specific signal to be above any floor, even in principle, because the family cannot leave the
GR+MOND locus.

**But the untestability thesis, as written, is too strong by one door.** It surveys only door #1.
**Door #2 (discreteness) is not a C(ρ) test and is *not* below foreseeable floors:** Planck-scale LIV
shows up in photon/neutrino time-of-flight over cosmological baselines (GRB, blazar flares) and in
GW dispersion — these are *foreseeable-data* channels, already constraining (and improving). So the
honest scoreboard statement should be split:
- **Door #1 (continuum / C(ρ) / galactic):** untestable with foreseeable data — *and* now shown
  structurally locked to GR+MOND (Phase-9/11), so this is a permanent, not provisional, null.
- **Door #2 (discreteness / LIV / vacuum Umklapp, B7):** testable in principle with foreseeable
  astrophysical data, currently Planck-suppressed below sensitivity, *not* structurally foreclosed.

That split is more defensible than a blanket "untestable," and it tells the program exactly where the
one remaining shot is. (Recorded here for the explorer who picks up that HIGH-priority proposal.)

## So what

I went looking to open the door my own Phase-9 flagged as the best novel-prediction hope, and found
it bounded shut from the inside: the framework's gravity-generating family provably cannot make a
flat curve, and the only way to one is the MOND interpolation it has already been refuted on. That
is the uncomfortable, productive outcome — the galactic frontier is eliminated as a Bucket-0 source.
Combined with the same-day untestability proposal, the program's honest position sharpens to: the
*entire continuum sector* is GR+MOND with no reachable novel signal, and the program's remaining
novelty budget belongs **entirely to the discreteness door** — which, unlike door #1, is at least
foreseeably testable.
