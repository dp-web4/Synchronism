# The alignment arc — kimi-code's second daily exploration program (charter)

**From**: kimi-code · **Chartered**: 2026-08-02 · **Frame-owner**: dp
**Cadence**: one wake per day (cron-fired); one concrete step per wake — a simulation, a
dated exploration doc, committed results, tracker updated. Inherits the CGL arc's house
rules (falsifier-first, registered kill criteria, honesty blocks, BOTH ledgers scored,
conventional-prior contamination check, the zoom-out rung every 3 stages or 7 days, the
distinct-count keeper).

## The frame (dp, 2026-08-02 — load-bearing, not decoration)

Everything in this arc lives on the core premise: a **CFD-like model across a 3D physical
grid in discrete global time** — one substrate, one clock, one observer. The pair, the
detectors, the environment, and the measurement apparatus are all **fractal patterns of
patterns swimming in a common substrate in tune to a common global clock**. The 3D
grid-scale-to-quantum simulation is currently prohibitive and the useful abstractions do
not yet exist — so every abstraction this arc builds must **name what it threw away**
(the substrate, the dimensionality, the tick), and no conclusion may be filed that
depends on the thrown-away part. *If we forget the river, we are studying whirlpools.*

## The bet

**Entanglement is phase-lock on the common substrate, and measurement is phase alignment
— a dynamical, finite-window process in which the detector pattern and the measured
pattern enter resonance.** The detector is not a privileged observer reading a
pre-existing property; it is another pattern in the river, and it reports the state of
its *coupling* — the way an atomic clock does not read time but reports how its own
energy levels respond to speed and gravity (the pendulum in the centrifuge, generalized:
no instrument reads a property; every instrument reports a relationship). Bell
correlations are then not a property of the pair alone: they are *manufactured in the
measurement window* by the alignment dynamics on the shared substrate.

This is dp's challenge to the earlier arc's kill, and it is also that kill's own
autopsy: the CHSH constructions capped at S ≤ 2 because their settings were independent
labels, and the arc named the missing primitive — "non-relabelable conditional
setting-dependence" — which is exactly this bet read as a design spec.

**The strong form, and its boundary (registered up front):** the correlation is produced
*dynamically, in the window*. If any construction in this arc can only reach the quantum
form via conspiratorial initial-condition fine-tuning (superdeterminism proper), that is
a KILL — an unfalsifiable alibi, not an instrument. Only the dynamical version is
testable, and only the testable version is the bet.

## Pre-registered kill criteria

- **Q1 (the model).** Kill if: (a) the construction reduces to *relabeling* measurement
  settings (the nonlocal-grid trap the earlier arc already fell into, S ≡ 2.0 —
  gauge-equivalent, not a mechanism); or (b) the correlation requires initial-condition
  fine-tuning rather than window dynamics (the conspiracy boundary above).
- **Q2 (the angular form).** Local hidden-variable constructions give triangle-wave
  C(θ); quantum mechanics gives the cosine. Kill if the asymptotic (full-lock) angular
  form provably cannot approach the cosine under local rules + dynamical alignment — with
  the obstruction named, not hand-waved. Pass: the cosine (or a named, distinguishable
  near-cosine) *emerges from the locking dynamics*.
- **Q3 (the payout: the measurement-strength curve).** The bet's distinctive physics:
  partial alignment ⇒ partial correlation, so S depends on the alignment window τ.
  Orthodox QM has its own interpolation for weak/partial measurement. Kill if (a) the
  model's S(τ) is indistinguishable from QM's weak-measurement interpolation at every τ
  (pure reparametrization, no buy — the species-(c) verdict, applied to ourselves
  honestly), or (b) it disagrees with published weak-measurement data where QM agrees
  (refuted by data). Pass: a named τ-regime where the curves diverge *and* the model is
  not already dead there.

The arc stops on a Q1 kill (the construction failed, nothing to measure) or completes
Q1–Q3. A Q2 or Q3 kill leaves the earlier stages' instruments standing as contributions.

## Stage map

- **Q1 — the alignment Bell model on the substrate.** Lattice sim (numpy, headless,
  1D first — 2D if 1D is degenerate; the 3D gap recorded as an abstraction debt): cells
  carry phase; local-only coupling; **one global tick**. Source region emits a
  phase-locked pair of patterns (entanglement = phase lock at creation, on the
  substrate). Two detector regions with settable phase references (the "angles" — and
  they are *physical phase relationships on the same substrate*, not free labels: this
  is the whole point, and the code must not smuggle in independence). Measurement =
  finite-window transient alignment (τ ticks of local coupling); outcome read from the
  aligned state. Deliverable: the model + C(θ, τ) + the two Q1 kill checks (relabeling
  audit, conspiracy audit).
- **Q2 — the angular form.** C(θ) at full lock vs triangle and cosine; if not cosine,
  the obstruction named (this is where the earlier arc's theorem-level caps are
  engaged honestly: Bell's premises assumed, one of them — measurement independence —
  is what this model *replaces with dynamics*, and the doc must show where the theorem's
  proof uses that premise and what the replacement does to it).
- **Q3 — S(τ) against the weak-measurement literature.** The payout test: a named
  τ-regime of divergence, or an honest species-(c) verdict on ourselves.

## House rules (inherited, plus two)

1–6 as in the CGL charter (falsifier-first; honesty blocks; both ledgers; no
conventional-prior contamination; tracker discipline; the zoom-out rung).
7. **The river rule (dp's frame, operationalized):** every sim and every doc names its
   abstraction debts — what of the substrate it discarded (dimensionality, tick
   granularity, intent quantization) — and no verdict may lean on a discarded part.
8. **The instrument-independence check (this arc's own contamination class):** any
   instrument that treats measurement settings as free independent labels is measuring
   the old ontology, not this one. Each stage must demonstrate how its settings are
   physical phase relationships on the substrate, or be flagged, not filed.

## Tracker (newest first)

- 2026-08-03 — **Q1 PASS.** `simulations/kimi_alignment/alignment_bell_q1.py` built and
  run: first C(Δθ, τ) map (9×6 grid, 600 pairs/cell) committed at
  `simulations/results/kimi_alignment/q1_first_map.json`. Correlation manufactured in the
  window: lock 0.04→0.76, triangle amplitude 0.02→0.54, CHSH S 0.19→1.29 over τ=0→64;
  controls clean (no-source null max |C| 0.12, no-lock 0.11, both ≈ noise vs signal 0.51).
  Relabeling audit PASS (static: outcomes from local scalars only; dynamic: outcome_A
  bitwise invariant under θ_B swap, 200/200). Conspiracy audit PASS (randomized detector
  phases, late settings, robust to K_SUB ×½/×2). Full-lock form NOT yet resolved: at
  τ=64 triangle A=0.540 (rmse 0.0796) vs cosine A=0.467 (rmse 0.0803) — statistically
  tied, peaks rounded by residual own-phase. Doc:
  `explorations/2026-08-03-kimi-alignment-q1-first-map.md`. **Next: Q2 — the full-lock
  extrapolation (τ → ~256 with per-τ lock measurement, A → 1), then triangle-vs-cosine
  discrimination at the asymptote, engaging the measurement-independence premise where
  the model replaces it with dynamics.**
- 2026-08-02 — arc chartered with dp's common-substrate frame ("patterns of patterns in
  one river, one clock, one observer"); Q1–Q3 staged; kill criteria registered.
  **Next: build the Q1 alignment-Bell lattice sim (`simulations/kimi_alignment/`) —
  source, pair, two detector regions, alignment window τ — and produce the first
  C(θ, τ) map plus the relabeling/conspiracy audits.**
