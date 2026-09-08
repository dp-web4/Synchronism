# MRH-validity, the sub-arc map: dp's epicycle diagnostic, and who can pick up what

**From**: kimi-code · **Date**: 2026-09-08 · **Acting on**:
`2026-09-08-kimi-mrh-validity-charter.md`, `2026-09-08-kimi-mrh-validity-first-results.md`
**Prompt** (dp, same day): *"the 'mo' in MOND, and CDM, are both un-instrumented
epicycle-like explanations where the standard model no longer matches observed
behavior (either there is matter we can't see, or matter behaves differently at
that scale). Map out next steps, sub-arcs if need be (that others could pick up),
then if something pulls, follow it."*

## dp's observation, formalized — and it is a diagnostic, not an analogy

In the program's terms the two camps are the two *kinds of response* to the same
signature — a residual that outgrew its equation's registered validity bound:

- **CDM reifies the excluded variable.** The horizon's included set (baryons,
  Newton/Einstein dynamics) is kept; the loss is declared to be a *thing* the
  horizon does not include. The epicycle question: **what instrument closes the
  horizon?** Direct detection, indirect detection, colliders — three decades of
  null results. In program language: *the posited excluded variable resists
  instrumentation.* That does not refute it (the variable may couple only
  gravitationally), but it prices it: every null year leaves the anomaly's
  entire information content living in its correlation with the included set.
- **MOND re-registers the equation's validity bound.** The included set is kept;
  the equation is declared to have been applied outside its registered horizon
  (strong-field, solar-system-tested gravity) and modified at the boundary. The
  program's honest counter-charge: **MOND commits the same violation it accuses
  GR of** — it has no registered validity on the cosmological horizon (CMB
  acoustic peaks, BAO, structure formation), and its community applies it there
  anyway. Under K1-style discipline, MOND's own bound statement is missing.
- **The epicycle parallel, made precise:** epicycles were a closed-horizon
  equation — exact until the instrument term changed (parallax). Both modern
  responses are *un-instrumented* in dp's sense: neither has produced the
  observation that changes closure for the other. The program's contribution is
  to make "closure" measurable instead of rhetorical: **the epicycle index** —
  the fraction of the anomaly's variance that is a function of the included set
  alone, instrument floor removed. If E → 1 under honest cuts and LOO
  discipline, the excluded sector is epicycle-class *regardless of which camp
  is right*: the instrument sees only included-set structure.

### Day-zero baseline (already run — the pull, followed)

`simulations/sparc_real_data/sa2_epicycle_index_dayzero.py` on the repo's local
SPARC copy (175 galaxies, 3317 points after stated cuts; fiducial M/L=1, no He
correction — simplifications recorded in the script):

```
anomaly   log10(g_obs/g_bar):        std 0.326 dex
residual  log10(g_obs/(g_bar*nu)):   std 0.181 dex   (McGaugh nu, a0 fiducial)
instrument floor (from errV):        std 0.069 dex
E_raw  = 1 - Var(residual)/Var(anomaly)          = 0.69
E_corr = 1 - (Var(res)-Var(inst))/Var(anomaly)   = 0.74
```

**~70–74% of the galactic anomaly's variance is priced by a one-parameter
function of an included variable, before any excluded sector is posited.** The
remainder above the instrument floor (~0.167 dex) is what either the
included-set correction series (the repo's 6-var model priced its version of
this at LOO R²=0.885 — different cuts, do not chain the numbers) or a genuine
excluded sector must account for. The structure by acceleration decade is
registered: E_local ≈ 0.25 in the deep-low-g regime, ≈ 0 in the Newtonian
regime (where the anomaly itself is mostly scatter).

Also woven in from the repo's own closure record
(`Research/proposals/rar_shape_test_closure_galaxy_program.md`): the γ=2
compander is refuted on this data (ΔBIC +184; ≈33 under conservative
correlation correction) and free-γ converges to 0.49 ≈ MOND's function. Under
the program that result reads: *the one-parameter included-set function is
measured, and its shape parameter prefers the MOND form — a statement about the
included set's structure, silent about the excluded sector.*

## The sub-arc map (each packaged for pickup: question, instrument, falsifier, first step)

**SA-1 — the exponent census** *(the program's empirical front; feeds/kills K2)*
For each domain where the repo measured γ: extract the fluctuation exponent and
test it against the CLT value 1/2. Four separable work packages:
- **SA-1a phonon** (Debye-class solids; literature data) — small, bounded.
- **SA-1b galactic** (SPARC residual structure vs N_corr scaling) — pairs with SA-2.
- **SA-1c quantum** (counting statistics cases from the quantum track).
- **SA-1d chemistry** (the ~2671-session dataset — does the repo's own largest
  corpus carry extractable exponents? The tautology-audit method transfers).
Falsifier (all): every domain reads 1/2 → K2 fires, γ dissolves into statistics,
and that is the reportable answer. *Pickup profile: any seat; SA-1a and SA-1d
are the most self-contained.*

**SA-2 — the closure test for the galactic horizon** *(dp's diagnostic; the pull
I am following)*
Registered rungs: (1) day-zero E — **done, this doc**; (2) E under honest
cuts + per-galaxy Υ* as an included-set parameter, with LOO discipline (the
published 0.1437-dex scatter with quality cuts should raise E — measure, don't
assume); (3) recompute the 6-var correction series on the SAME point cloud as
stage 2, so the two stages of the variance budget chain legitimately — the
deliverable is a single audited budget: anomaly = instrument + included-set
function + included-set series + remainder; (4) the remainder's structure: does
it correlate with any included-set variable left, or is it sector-demanding?
(5) the honest MOND-side rung: write MOND's own validity-bound statement
(register where its equation is measured vs where it is merely assumed — CMB,
BAO, clusters) so both camps stand priced. Falsifier: if the remainder at rung
(4) is structureless noise, the excluded sector is epicycle-class on current
instruments; if it carries included-set-invisible structure, the closure test
FAILS and the excluded-variable posit earns its rent. *Pickup profile: this is
the arc's heaviest computation; shared with the cosmology-track seats.*

**SA-3 — non-Gaussian loss mapping** *(theory)*
One worked non-Gaussian system where the Fano-direction bound is loose vs the
true price of exclusion; calibrates how far the Gaussian anchor (e^{2·CMI})
overstates precision. *Pickup profile: a theory-leaning seat; small, clean.*

**SA-4 — the epicycle case made quantitative** *(history-of-science instrument
case study)*
Actual naked-eye-era positional precision vs the information content of the
excluded heliocentric frame: compute the closure claim currently argued in
prose. Bounded, charming, self-contained. *Pickup profile: ideal first task for
a new seat; produces a citable number.*

**SA-5 — the web4 MRH bridge (K3)**
One definition of MRH that compiles in both worlds: context-scoping graph ↔
restriction map/σ-subalgebra. Deliberately after SA-1..3 — the physics earns
the bridge, not vice versa. *Pickup profile: a web4-native seat.*

**SA-6 — the two-level time resonance** *(long shot, parked)*
Tomita-Takesaki modular time from restriction vs the substrate track's
two-level time ontology. Opens only if SA-1..3 land; touching it earlier is
contamination risk.

## What pulls next (my own queue, in order)

1. **SA-2 rung 2** — E with quality cuts + LOO Υ*, the honest version of the
   day-zero number.
2. **SA-1b** — the galactic exponent measurement, sharing SA-2's pipeline.
3. Writing both up as the program's second results doc.

The rest is genuinely up for grabs — any seat picking up an SA should file a
dated doc in this directory naming the sub-arc and its registered falsifier
before running anything, per the program's house rules.

## Provenance

kimi seat, interactive session with dp, 2026-09-08. Instrument:
`simulations/sparc_real_data/sa2_epicycle_index_dayzero.py` (deterministic;
rerun for every number quoted). Sources read: the two program docs,
`Research/proposals/rar_shape_test_closure_galaxy_program.md`, SPARC data
format. No holdout analogue exists in this program; cuts are stated and the
script is the audit.
