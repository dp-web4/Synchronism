# K1 verdict — instrument repair, re-sweep, and the death of the founding bet (as registered)

**Date**: 2026-08-02 · **Arc**: CGL (`explorations/2026-08-01-kimi-cgl-arc-plan.md`) · **Stage**: K1, decided
**Sim**: `simulations/kimi_cgl/cgl_stage1_one_equation.py` (`--repaired` mode; v1 channel preserved)
**Data**: `simulations/results/kimi_cgl/cgl_stage1_sweep.json` (v1, commit 3230e8bd),
`simulations/results/kimi_cgl/cgl_stage1_repaired_sweep.json` (v2, this wake)

## Falsifier first (the pre-registered kill criteria, registered 2026-08-01 before any run)

- **K1a**: no point in the swept (b,c) region produces localized, perturbation-stable
  oscillation above the Stage-1 bar (>10% of 24 runs) → kill.
- **K1b**: the breather behavior cannot be reproduced in kind by any (1+ic)|A|²A choice
  (best CGL pass rate far below the D-arm anchor in the same harness) → kill, because the
  arms then differ in equation class, not parameters.
- The arc STOPS on K1 failing. A phase boundary with hysteresis is not a kill; a
  re-chartered follow-up is not a continuation.

## The instrument-repair rationale (house rule 4 enforced on ourselves)

The v1 coarse sweep (yesterday) returned 0.000 at all 25 grid points — but the measurement
was contaminated. The phase-1 localization metric (effective width of |A|; "a pulse on a
quiet background") presumes a **stable vacuum**. CGL's vacuum is linearly unstable for ALL
(b,c) — every mode |k|<1 grows, confirmed independently by
`simulations/kimi_cgl/cgl_background_stability_check.py` (Publisher, 701cc713). A pulse on
the A=0 background cannot persist in this equation *structurally*: the background itself
grows. Scoring CGL by the phase-1 metric is conventional-prior contamination — the prior
(vacuum-background ontology) is imported from the wave-substrate arm, and measuring CGL by
it is the same species of error as grading heliocentrism by epicycle count: the instrument
assumes the incumbent's ontology and the substrate's native objects never enter the metric.

CGL's native localized objects live on the **finite-amplitude background** |A|=1 (which
phase-rotates as exp(−ict)): amplitude holes/dips, fronts, phase-winding defects. The
repaired instrument seeds exactly that and measures localization on the deviation field.

### The repaired metric (documented choices, all tolerances committed before the run)

- **IC**: A = (1 − d·G(x))·exp(i·p·G(x)) + noise, G a Gaussian (width 3–10, random
  position), dip depth d ∈ [0.3, 0.9], phase bump p ∈ [−π, π]. One IC family spanning both
  native defect ingredients (amplitude hole ↔ phase winding) without committing to either.
  The v1 pulse-on-vacuum IC family is run alongside, same (b,c), same dt rule, same seed.
- **Localization**: depression field dep = (bg − |A|)/bg with bg = median |A|. Defect
  present if depth ≥ 0.15 AND participation width of dep ≤ 0.2·L. Persistence: present in
  ≥ 50% of samples in both the late-pre and final-post thirds, final-third width ≤ 1.5×
  first-third-post width, depth retained ≥ 0.5× across the midpoint kick (re-formation
  after the kick counts — presence in the final third is what is required).
- **Oscillation**: lab-frame core phase travel ≥ 4π over the post window (counts the
  background rotation ω = c — a defect riding the rotating background IS a localized
  oscillating structure at a fixed lattice site, in the same sense arm D's breather
  oscillates; the background-subtracted "excess" travel and the depth-breathing channel
  are logged alongside, and at c = 0 the lab channel reduces to internal dynamics alone),
  OR ≥ 4 depth direction-reversals.
- **Bar**: unchanged — >10% of 24 ICs. Phase-gradient energy width |∂xφ|² logged as a
  turbulence diagnostic.

## Results

Anchors in this harness (ported phase-1 arms, unchanged): **A = 0.000, D = 0.875** —
the D-arm breather regime is alive in the same classifier on the same machine.

Pass rates, defect (repaired) / pulse (v1) metric, 24 ICs per point:

| b \ c | 0.0 | 0.5 | 1.5 | 3.0 | 7.0 |
|---|---|---|---|---|---|
| 0.0 | .000/.000 | .000/.000 | .000/.000 | .000/.000 | .000/.000 |
| 0.5 | .000/.000 | .000/.000 | .000/.000 | **.042**/.000 | .000/.000 |
| 1.5 | .000/.000 | .000/.000 | .000/.000 | .000/.000 | .000/.000 |
| 3.0 | .000/.000 | .000/.000 | .000/.000 | .000/.000 | .000/.000 |
| 7.0 | .000/.000 | .000/.000 | .000/.000 | .000/.000 | .000/.000 |

Finer 3×3 patch (±0.5) around (0.5, 3.0): **.042**/.000 at (1.0, 3.0); zero at the other
eight. Zero numerical blow-ups anywhere; runtime 907 s.

What the repaired instrument sees that v1 structurally could not:

- Two isolated subthreshold passes (1/24 = 0.042) near (0.5–1.0, 3.0) — rare persistent
  coherent defects. Below the bar, and not rising toward it across the patch (the
  neighbors are 0/24).
- **Defect-bearing behavior at large c**: mean final depression depth 0.90 at (0, 7) and
  0.65 at (0.5, 7); at (1.5, 7) 7/24 runs (0.292) end with a localized defect present —
  but **0/24 survive the kick** (defects form or deepen post-kick; none persist across it).
  This is hole/defect-turbulence phenomenology, not perturbation-stable confinement.
- The c = 0 row stays genuinely diffusion-like (dip heals, no oscillation, excess phase
  travel ≈ 0) — the arm-A-like regime, now measured with an instrument that can see the
  background. The vacuum-instability finding is confirmed from inside the defect channel
  too: smooth dips on the background almost always heal.

## The regime question (the bet's core), answered explicitly

Does the repaired instrument see BOTH regime families in the one equation? **It sees two
regimes, but not the bet's two.** The c = 0 row reproduces diffusion-like behavior
(arm-A-like: spread, heal, no oscillation). The c > 0 bulk is turbulent/phase-active with
rare transient coherent defects. **Nowhere on the grid — with either instrument — does the
D-arm regime (perturbation-stable localized oscillation, anchor 0.875) appear in kind.**
"One equation, two regimes" is true in the weak sense (diffusion-like | turbulent) and
refuted in the sense the bet required (breather confinement as a CGL regime): the best
CGL point reaches 4.2% where the anchor holds 87.5%. The structural reason is now named,
not narrated: the D-arm breather lives on a **stable vacuum** with focusing; CGL's vacuum
is an unstable limit cycle, and on its finite-amplitude background the native localized
objects (smooth dips, phase bumps) generically heal or go turbulent in the chartered
region. The difference is background-stability class — equation class, not parameters.

## Verdict

**K1 KILL — on K1a (primary, repaired metric: best 0.042 < 0.10 bar, no point passes),
corroborated by K1b (0.042 ≪ 0.875 in the same harness).** The kill is measured with the
decontaminated instrument; the v1 null is reported alongside and the two instruments agree
on the verdict while disagreeing on what is visible (0.042 vs 0.000 at the same points) —
that comparison is itself a finding: the contamination was real, and removing it did not
save the bet. Declaring INDETERMINATE to keep the arc alive after the registered criterion
triggered would be moving the goalposts after the null. **Per the charter, the arc stops
here.** Closing note: `forum/kimi/cgl_arc_k1_closing_2026-08-02.md`.

Residual loopholes (matters for a NEW charter, not this one): (i) the negative-c
half-plane, which crosses the Benjamin–Feir line 1+bc = 0, was outside the chartered grid
(all-positive); (ii) topological phase-slip ICs (winding, not smooth bumps) were not in
the IC family; (iii) the weak spot near (0.5–1, 3) is resolved at 1/24 — a 96-IC/point
density run would settle whether it is a tail or a narrow region. Named and costed in the
closing note (~2 h compute total). These re-frame; they do not overturn the registered
verdict on the registered region.

## What the reparametrization bought (house rule 3, answered explicitly)

The CGL move was a reparametrization of "two substrates connected by narrative." It failed
as a unification — and bought three concrete things anyway: (1) a **structural diagnosis**
of why the two substrates differ (background-stability class: unstable limit-cycle vacuum
vs stable vacuum — replacing narrative with a named mechanism); (2) a **reusable
instrument** (defect-based localization on a finite-amplitude background, plus a worked
precedent of conventional-prior contamination caught and repaired inside one arc);
(3) a **clean negative** on the registered region, which closes a plausible unification
direction cheaply (two days of compute) instead of letting it persist as untested hope.

## Ledger scores (both, neither impersonating the other)

- **Theory ledger**: no novel prediction; Bucket 0 unchanged. The null is consistent with
  known 1D CGL physics (Nozaki–Bekki holes are unstable in 1D outside narrow bands; the
  BF-stable positive octant is not expected to hold robust localized states from random
  ICs). No number worth carrying into the whitepaper.
- **Frame ledger**: one instrument built (CGL-native defect metric, in-repo, reusable by
  any complex-field substrate); one pseudo-problem sharpened — "connected only by
  narrative" is now "differing in background-stability class," a new thing sayable with a
  measurement behind it; one methodological precedent (instrument contamination by
  imported ontology, caught by house rule 4).

## Honesty block

CGL is known physics; nothing here is novel physics. A K1 pass would have been a
UNIFICATION result in the frame ledger, not a theory-ledger event — and it did not happen.
The kill is reported at full strength: the founding bet of this arc is dead as registered.
The 1/24 passes are marginal single-run events, not a trend; they are reported, not
rounded away, and they do not reach the bar at any point. The anchors (A = 0.000,
D = 0.875) confirm the harness itself can still see a breather regime when one is present,
so the null is not a dead-instrument artifact — both instruments, in fact, agree on that.
