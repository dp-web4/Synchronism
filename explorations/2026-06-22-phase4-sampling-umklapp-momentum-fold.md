# Phase-4 — the universe as a sampler: discreteness → Umklapp (momentum mod G) (2026-06-22)

**Status:** `[ACTIVE-MRH]` — tests dp's hypothesis that a base sampling/tick rate could underlie
quantization and high-energy effects, by treating the substrate as a discrete-time sampler (DSP)
and asking what it must then do. Yields a clean, controlled, **testable** prediction.
**Sim:** [`simulations/phase4_sampling_umklapp_momentum_fold.py`](../simulations/phase4_sampling_umklapp_momentum_fold.py) · result: `simulations/results/phase4_sampling_umklapp_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous — from the dp↔GPT Nyquist thread.

## The frame (dp)

Not cell-properties yet — the stage before: **a grid of cells advancing state on global clock
ticks.** The question: *if the substrate is a discrete-time state machine sampled at a base rate,
are quantization / Planck-scale / high-energy effects what you'd expect to fall out?* This rests
on a **fractal-leverage bet**: a Planck-scale discrete-time grid, if it is a state machine, should
behave *meaningfully like a digital signal processor* — the similarity holds *within the
discrete-sampling MRH*, which is the MRH under test. Nyquist and Umklapp are **invited**
(MRH-appropriate analogies), not smuggled. The one move *not* allowed: "the cutoff co-dilates with
patterns" — the sampling cutoff is the grid's property, **pattern-unaware**, by the model's
current statement. (That statement is about the model as specified, no more, no less.)

## The test

A sampler has a band limit: spatial `|k| ≤ π/a` (the first Brillouin zone = the Nyquist
wavenumber). A *linear* signal below it samples faithfully. A *nonlinear interaction* **mixes** —
patterns at `k1, k2` generate sum content at `k1+k2`. If that exceeds the band limit the sampler
**cannot represent it** — it **aliases**, and on a grid the aliasing of crystal-momentum *is*
**Umklapp**: the daughter appears at `k1+k2 − G`, `G = 2π/a`. Directional momentum needs a
**complex** field `e^{ikx}` (a real `cos` has a sign-symmetric spectrum and can't show direction)
— the same Phase-1.6 lesson, a third time.

## Result (complex field, χ² nonlinearity)

| | |
|---|---|
| **Static fold** | daughter's signed `k` = `wrap(2k)` **exactly** at every `k`; for `k/π > 0.5` it **reverses to negative `k`** (forward pump → backward daughter). Directional Umklapp. |
| **Two-pump collision** | `q1+q2 = 3.5 > π` → daughter at **−2.79 = wrap(3.5)**, not the naive 3.5. Momentum folded. |
| **Refinement control (agent-zero)** | 4× finer grid (band `π → 4π`) → the *same* interaction does **not** fold (daughter = 3.5, unaliased). **The fold vanishes as the Nyquist limit rises** — so it's a genuine sampling-cutoff effect, not coded-in physics. |

(A dynamical-evolution confirmation was attempted and **discarded** — the explicit complex
integrator went unstable and didn't grow the daughter; agent-zero says don't ship a check that
didn't work. It's redundant anyway: aliasing is a property of the *sampled representation*, so the
static χ² result *is* the exact mechanism, not "just arithmetic.")

## What the analogy buys — the testable prediction

If the universe is a discrete-time grid (the fractal Planck-DSP bet), then continuum QFT's *exact*
momentum conservation (Noether, from continuous translation symmetry) is replaced by conservation
**mod `G ~ Planck momentum`** (only discrete translation symmetry survives). **Vacuum Umklapp:**
extreme-momentum-transfer interactions should show apparent momentum non-conservation by `G` /
anomalous back-scatter. This is:
- a **distinct LIV channel** from Phase-2's dispersion-LIV (that one modifies free `ω(k)` near the
  zone edge; this one folds *interaction* momentum), and
- (dp's point) reachable at **lower per-quantum energy**, because nonlinear mixing **stacks**
  momentum toward the cutoff even when each quantum sits well below it.

So Phase-2 and Phase-4 are two faces of one discreteness: free-flight dispersion-LIV, and
interaction Umklapp-LIV.

## Honesty

Umklapp is **textbook condensed-matter** (lattice momentum mod `G`) — **not novel**. The novel,
falsifiable claim is its application to the **vacuum substrate** (momentum mod a *Planck-scale*
`G`). It is currently **unobserved** — momentum conservation is exact to high precision in every
tested interaction — so it is a **risky Bucket-1 bet**, Planck-suppressed for single interactions,
**not a confirmation**. Registered as **B7**. **Bucket 0 unchanged.** And the whole result is
**conditional on the fractal bet** that a Planck grid behaves like a sampler — that is the
hypothesis under test, stated as such, not assumed.

## Open tension (carried, not resolved)

The cutoff is a substrate property and (per the model) does **not** co-dilate with patterns — so
"how close is an interaction to the cutoff" is **frame-dependent** (it depends on velocity
relative to the substrate rest frame, the universal-clock frame). Standard physics makes that
frame-*invariant*. This is the **same preferred-frame / Lorentz-invariance tension** as the
relativity thread: the sampling cutoff's frame-(in)variance *is* the Lorentz-invariance question
for the discrete substrate. Unresolved, and not papered over.
