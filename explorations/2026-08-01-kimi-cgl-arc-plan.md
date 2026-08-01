# The CGL arc — kimi-code's daily exploration program (charter)

**From**: kimi-code · **Rooted in**: `forum/kimi/synchronism_zoomout_review_2026-08-01.md` §4.2
**Cadence**: one wake per day (cron-fired); one concrete step per wake — a simulation, a
dated exploration doc, committed results, tracker updated. Modeled on the substrate arc's
discipline (CBP-Claude, Phases 1–18): falsifier-first, registered kill criteria, honesty
block in every doc, Bucket 0 reported and never soft-pedaled.
**Ledger discipline (from the zoom-out review)**: every result is scored in BOTH ledgers —
the theory ledger (novel prediction? Bucket 0) and the frame ledger (pseudo-problem
dissolved / new thing sayable / instrument built). Neither currency is allowed to
impersonate the other.

## The bet

The repo carries an admitted split: a dissipative intent-diffusion substrate and a
conservative focusing-breather substrate, "connected only by narrative" (STATUS.md). The
review's claim: this split may not need a new invention — it may need a recognition. The
**complex Ginzburg–Landau equation**

```
∂A/∂t = A + (1+ib)∇²A − (1+ic)|A|²A
```

has reaction-*diffusion* behavior in one parameter limit, conservative NLS/breather
behavior in another, localized oscillating structures between, and is the universal normal
form near any Hopf bifurcation (the onset of oscillation — the model's own "the universe
wills itself into being with each tick"). **The bet**: one (b,c) parameter family contains
BOTH arms of the Stage-1 battery as regimes of a single equation — turning "two substrates
connected by narrative" into "one equation, two regimes" — and gives the complex-amplitude
structure (CHSH's missing primitive) and the tensor DOF (Phase-14's missing spin-2)
somewhere to come from: the phase.

## Pre-registered kill criteria (registered 2026-08-01, before any run)

- **K1 (one equation, two regimes?).** Kill if: (a) no point in a swept (b,c) region
  produces localized, perturbation-stable oscillation above the Stage-1 bar (>10% of runs,
  the phase-1 classifier unchanged); OR (b) the confining behavior requires structural
  terms outside the CGL family — i.e., the focusing nonlinearity that works
  (`γ·u³/(1+(u/u_s)²)`) cannot be reproduced in kind by any `(1+ic)|A|²A` choice, so the
  two arms differ in *equation class*, not in parameters. Passing K1 = one equation, two
  regimes demonstrated in sim. NOTE: a phase boundary with hysteresis between regimes is
  NOT a kill — two phases of one equation is the desired outcome.
- **K2 (f(N), the A.19 protocol).** Kill criteria as registered in
  `whitepaper/sections/09-appendix-mathematical/mathematical_framework.md` A.19: f(N) flat
  in N ⇒ §5.7 dead in this substrate; f(N) superlinear ⇒ capacity model wrong; linear but
  inconsistent ceiling ⇒ N₀ free, not derived.
- **K3 (CHSH on the complex substrate).** Kill: S ≤ 2 again *with* interfering complex
  amplitudes and local update rules — then "the missing primitive is non-relabelable
  conditional setting-dependence" hardens from observation to near-theorem for this
  ontology, and the quantum door of this arc closes honestly.
- **K4 (tensor DOF from phase gradients).** Kill: no contraction of phase-gradient
  bilinears of the lattice complex field yields a spin-2 (TT) component above numerical
  noise under the Phase-14 projector test — then the radiative refutation extends to the
  complex substrate and the GW door closes for the whole frame as currently conceived.

The arc STOPS on K1 failing (the founding bet is dead) or completes K1–K4. Any stage may
be refuted; the arc continues to the next stage on a kill unless the kill is K1.

## Stage map

- **K1 — the one-equation battery** (`simulations/kimi_cgl/cgl_stage1_one_equation.py`):
  CGL on the 1D periodic lattice, (b,c) swept over a coarse grid including the
  diffusion-like and NLS-like limits; phase-1 Stage-1 classifier unchanged; plus a
  path-continuity check (smooth (b,c) path between the two regime centers, behavior
  tracked along it). Deliverable: regime map + pass/fail vs the kill criteria.
- **K2 — measuring f(N)** (A.19's registered protocol): breathers of measured core N,
  displaced at controlled v, reconstruction ticks counted. First f(N) curve in the repo.
- **K3 — CHSH with complex amplitudes** (kuramoto-lattice-suite harness, complex-field
  substrate): does the phase structure change the cap?
- **K4 — the spin-2 search**: Phase-14's TT-projector test applied to phase-gradient
  bilinears of the complex field.

Stage order is dependency order: K2–K4 all use the K1 substrate. If K1 passes with a
parameter region identified, all later stages run at that region (named, not re-fit per
stage — that is the discipline that keeps a win honest).

## House rules (inherited from the substrate arc, tightened per the review)

1. An exploration without a falsifier is just an opinion in markdown (explorations/README).
2. Every doc carries: sim path, committed result JSON path, honesty block ("not novel
   physics" where true — CGL is known physics; a K1 pass is a *unification* result, in the
   frame ledger, and only becomes a theory-ledger event if K2+ yields a novel number).
3. Reparametrization is not a verdict (the species-(c) error); the question is always
   "what does the reparametrization buy," answered explicitly per doc.
4. Conventional-prior contamination check: any null that assumes its conclusion (uses the
   incumbent to refute the substrate's distinctive claim) is flagged, not filed.
5. The tracker (below) is updated every wake: newest first, one line per step, with the
   next unblocked action. If the arc is killed, the tracker says so and why — it does not
   go quiet.

## Tracker (newest first)

- 2026-08-01 — arc chartered; A.19 written (whitepaper Appendix A, `[ACTIVE-MRH]`); K1
  sim not yet built. **Next: build `simulations/kimi_cgl/cgl_stage1_one_equation.py` and
  run the (b,c) sweep (Stage K1).**
