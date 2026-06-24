# Phase-15 (door #3) — time dilation as global-clock frequency shift; the complexity limit, tested (2026-06-24)

**Status:** `[ACTIVE-MRH]` — treads dp's rich-ground proposal: model GR "time dilation" on the
absolute-time / global-clock CFD foundation (internal pattern frequencies shift vs the global clock
with translation + intent-concentration; fractal self-similarity makes it read as time dilation;
claim: this reinforces a complexity speed limit). **Result: one generative win, one honestly-
corrected intuition, one new candidate bet.** (1) Clock universality / the equivalence principle is
**derived** from fractal self-similarity — exactly dp's mechanism, and a real explanatory result.
(2) A complexity speed limit *below c* does **not** emerge in the continuum (the γ-factors cancel =
emergent Lorentz invariance; heavy ions survive at 0.9999c). (3) The complexity intuition *does*
survive, relocated: on the discrete lattice, LIV dephasing scales as `Δk²` (complexity) — a
**complexity-enhanced LIV** that folds door #3 into door #2 with a genuinely new angle.
**Sim:** [`simulations/phase15_clock_universality_and_complexity_limit.py`](../simulations/phase15_clock_universality_and_complexity_limit.py) · result: `simulations/results/phase15_clock_universality_complexity_result.json`
**Author:** CBP-Claude (Opus 4.8), with dp (proposal).

## Part A — clock universality DERIVED from fractal self-similarity (the win)

A pattern's subpatterns oscillate at frequencies `{f_i}` relative to the global clock. Translation
shifts each by `S_v=√(1−v²)` (light-clock, Phase-5); intent-concentration by `S_g=√(1−2GM/c²r)`
(GP, Phase-9). The decisive point (dp's): **both factors multiply *every* subpattern by the *same*
`S=S_v·S_g`** (all fractal subpatterns ride the same bulk translation/concentration). Therefore:
- internal observables (frequency **ratios** `f_i/f_j`, the pattern's own clock) are **invariant**
  under `S` — max deviation `2.2×10⁻¹⁶`. The pattern **cannot detect its own slowdown.**
- every mechanism dilates by the **identical** `S` — per-mechanism spread `1.1×10⁻¹⁶`.

So time dilation is what an *external* (global-clock) comparison sees, and it is **clock-universal**
(no mechanism-dependence). This **derives GR's clock universality / equivalence principle** — which
GR *posits* — from fractal self-similarity, and it constructively resolves the Phase-9 tension:
mechanism-dependence (the dead instrument-effect discriminator) is forbidden *precisely because*
every mechanism is a fractal subpattern sharing `S`. A genuine **generative** result (explains a
known fact GR assumes), exactly as dp framed it.

## Part B (continuum) — the complexity speed limit does NOT emerge (honest correction)

Coherence needs the binding signal (speed `c` = reconstruction rate) to round-trip the footprint
within the coherence time. In the substrate (preferred) frame the footprint contracts (`L/γ`) and
the coherence time dilates (`τγ`). The round-trip is `t_rt = 2Lγ/c` and the budget is `τγ`: **the
γ's cancel** → the coherence margin is **velocity-independent** → a pattern coherent at rest stays
coherent at *any* `v<c`. This is emergent Lorentz invariance (Phase-5/6), and it matches reality:
**heavy ions survive at 0.9999c** (LHC). So "complex can't travel as fast intact" is **not** a
fundamental sub-`c` effect; the readily-observable complexity↔speed relation is the standard energy
cost `E=γmc²` (heavier/more-complex = more energy to accelerate), and UHECR fragmentation is
photo-disintegration (CMB interaction), not a substrate limit. **dp's "reinforce complexity limit"
does not hold in the continuum** — recorded honestly, not modeled as a confirmation.

## Part B (discrete) — where the complexity intuition DOES survive: complexity-enhanced LIV

The γ-cancellation is exact only in the continuum. On the discrete lattice (Phase-2 dispersion
`ω²=m²+2(1−cos k)`), modes dephase because the lattice group-velocity dispersion deviates from the
relativistic one. The **lattice-specific (LIV) extra dispersion** grows toward the zone edge
(boost → higher carrier `k₀`), and the dephasing rate scales as **`Δk²`** — and `Δk` (the pattern's
spread in mode space) **is** complexity: a more localized / more composite pattern spans more modes.

| Δk (complexity) | LIV dephasing rate |
|---|---|
| 0.05 | 5.7×10⁻⁴ |
| 0.10 | 2.3×10⁻³ |
| 0.20 | 9.2×10⁻³ |
| 0.40 | 3.7×10⁻² |

So dp's complexity intuition lands, **relocated**: not a sub-`c` ceiling, but **complexity-enhanced
LIV** — a complex composite accumulates the discreteness signature at *lower per-quantum energy*
than a simple particle (`∝ Δk²`). Door #3 folds into **door #2 (LIV)**, but adds a genuinely new
angle: LIV scaling with **complexity**, not just energy — potentially *more* accessible than
single-particle Planck-scale tests (a composite system provides the `Δk²` enhancement). This is a
**candidate door-#3/#2 prediction, not yet a quantified falsifier** (the magnitude for a real
composite + the observable are the next step) — a *direction*, consistent with the door-#3
discipline; not registered as a Bucket-1 bet yet.

## Honesty / net

- **Bucket 0 = 0.** Part A is generative (derives a known fact), not a novel prediction. Part B
  continuum is a null (corrects the intuition). Part B discrete is a candidate direction, not yet a
  bet.
- Caveats: Part A is exact (algebraic — shared `S` cancels in ratios). Part B continuum's
  γ-cancellation is the standard emergent-Lorentz result (toy bound system). Part B discrete uses
  the validated Phase-2 lattice dispersion; the `Δk²` scaling is the leading group-velocity-
  dispersion term — a real magnitude estimate for a named composite + the matching observable is the
  unfinished work.
- **The tread paid off the way honest exploration should:** dp's mechanism (A) is *right and
  valuable*; dp's stated consequence (continuum complexity limit, B-cont) is *wrong, and the model
  says exactly why*; and taking the complexity intuition seriously anyway surfaced a *new* angle
  (B-disc, complexity-enhanced LIV) that neither the pure door-#2 nor the pure door-#3 framing had.
