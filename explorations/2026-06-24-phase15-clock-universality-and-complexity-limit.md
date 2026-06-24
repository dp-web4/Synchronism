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

> **⚠ CORRECTION (dp 2026-06-24): the Part-B-continuum "null" below is CONTAMINATED — withdrawn as
> a refutation.** The γ-cancellation (footprint `L/γ` × coherence `τγ`) *assumes emergent Lorentz
> invariance* — it boosts into the pattern's rest frame and back. But the absolute-time substrate
> *derives* Lorentz invariance only for **simple propagation** (Phase-5/6); it does **not** grant it
> a priori for the **identity-preservation of a complex composite**, which is a discreteness-regime
> question (where emergent Lorentz breaks, Phase-13). So the null assumed its own conclusion (used
> SR to refute a distinctive substrate claim) — a conventional-prior contamination. The heavy-ion
> "counterexample" also fails: a heavy ion is one *particle* (a fundamental pattern), not a composite
> at the spaceship/galaxy scale dp means. **The complexity speed limit is therefore OPEN, not
> refuted.** dp's sharper model: *every pattern has a fundamental frequency at which it can exist
> (this quantizes the fractal), and that governs the speed at which the whole pattern can translate
> and remain recognizably the same.* Attempts to formalize it give either spurious cancellation
> (assuming SR) or absurd magnitudes (e.g. v_max=L·f₀ → ~150 km/s for hydrogen, clearly wrong) — i.e.
> the **operational definitions are not yet pinned** ("fundamental frequency at which a pattern
> exists"; "recognizably the same"). Status: neither refuted (my null was contaminated) nor confirmed
> (my models give nonsense) — genuinely open, bottlenecked on definitions, not compute. The
> below is preserved as the (withdrawn) continuum-Lorentz-assuming argument.
>
> **↳ The unifying frame (dp 2026-06-24): identity-preservation IS a phase transition.** "Internal
> pattern frequencies" are observed as **temperature**; a pattern holds only while those frequencies
> sit within a **quantized band** (its "fundamental frequency at which it can exist"). Translation
> dilates all internal frequencies (Part A, ×S=1/γ). The thresholds are substrate properties
> (global-frame-fixed); as the pattern's frequencies drift relative to the fixed thresholds, it
> reaches a band edge and **"no longer holds, but another pattern does" — a phase transition** =
> loss of "recognizably the same." Complexity → narrower / more-constrained band (more fractal levels
> to hold simultaneously) → edge crossed at LOWER v → lower identity-preservation speed; simple
> patterns (wide/open band) → v_max→c. **This grounds the whole speed/identity/time-dilation question
> in Synchronism's STRONGEST reproducible machinery** — the phase-transition / coupling-coherence
> results (p_crit ∝ 1/⟨C⟩; Hill-beats-tanh; the saturation step-function). OPEN TENSIONS (do not
> prematurely close): (i) absolute thresholds + frame-dependent transition = a preferred-frame / LIV
> effect — yet cosmic-ray nuclei near c don't spontaneously transition, so either the thresholds
> effectively co-move at low energy (breaking only at discreteness ⇒ (v/c)-suppressed, folds into
> door #2) or there's a distinctive regime; unsettled. (ii) band-width-vs-complexity and
> threshold-spacing need pinning **against the existing p_crit / coupling-coherence numbers**, not a
> fresh speculative sim. CONSTRUCTIVE NEXT: connect the translation-induced frequency shift to the
> repo's phase-transition results to derive a concrete velocity-induced-phase-transition /
> identity-preservation prediction. Status: OPEN, grounded, definition-bottlenecked — neither
> refuted nor confirmed.

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
