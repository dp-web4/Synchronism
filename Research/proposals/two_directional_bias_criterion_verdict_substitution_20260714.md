# Two New Failure-Mode Instances Sharpen the Reflexivity Finding: Criterion-Verdict Substitution and Explanatory-Prose Over-Promotion

**Date**: 2026-07-14
**Source**: Site visitor track, four-persona pass 2026-07-14 (Pass 3 Grad Student, Pass 4
  Leading-Edge Researcher), verified against site source before this write-up
**Status**: Proposal — extends
  `directional_bias_law_fails_null_reflexivity_is_predictor.md` (2026-07-09), does not
  supersede it
**Reproduces via**: `synchronism-site/visitor/logs/2026-07-14.md` (Pass 3, Pass 4); site
  source at `src/app/tier-1-existing/page.tsx:69-72`, `src/app/honest-assessment/page.tsx:460-509`,
  `src/app/core-idea/page.tsx:70,196-197`, `src/app/research-philosophy/page.tsx:230-231`

## What the 2026-07-09 proposal established

Self-referential statistics (session counts, contribution counts, demotion rates) break
100% of the time (3/3), because they have **no primary source outside the loop to walk
to** — the mechanism is absence of ground truth, not motivated reasoning. Physics
statistics (CHSH, RAR, DESI, LIV…) break only 15% (4/27), because they ship with a fixed
script or a citable dataset.

## Today's finding: a break *inside* the low-rate physics class, with a different mechanism

TEST-04a (DESI RSD growth) has a primary source — two DESI papers, both cited on the site
(arXiv:2411.12021 for σ₈, and the framework's own registered criterion). It should be a
low-risk statistic under the reflexivity model. It broke anyway, and in the same
over-refuting direction as every other self-referential break:

- **Registered kill criterion** (`tier-1-existing.ts:70`): `fσ₈(z=0.51) > 0.46 rules out at
  >3σ; >0.45 disfavors at >2σ`.
- **Verdict delivered** (`tier-1-existing.ts:72`, `honest-assessment.tsx:460-473`):
  `σ₈ = 0.841 ± 0.034` — a **different statistic** — badged "Kill Criterion Triggered."
- **On the actually-registered statistic**: `honest-assessment.tsx:171` gives LRG1
  `fσ₈/(fσ₈)_fid = 1.16 ± 0.13` at the fiducial 0.474 stated on the same page → `fσ₈ =
  0.550 ± 0.062`. Distance from the 0.46 threshold is **~1.5σ**, not the >3σ the criterion
  demanded to "rule out." The badge's own >2σ clause (0.45) is met, but the >3σ clause the
  "Kill Criterion Triggered" label implies is not.
- **σ₈ is GR-conditioned**: it is inferred from a full-shape EFTofLSS fit whose
  perturbation kernels assume GR growth — the site's own cited counterterm literature
  (Cabass, Simonović, Zaldarriaga) already concedes this moves fσ₈ by 1–2σ within ΛCDM.
  Using a GR-conditioned amplitude to falsify a modified-growth model is circular in a way
  the site does not flag.
- **DESI's own purpose-built modified-gravity analysis is uncited**: Ishak et al.,
  "Modified Gravity Constraints from the Full-Shape Modeling of Clustering Measurements
  from DESI 2024," arXiv:2411.12026 (JCAP 09 (2025) 053), constrains exactly the parameter
  a growth-suppression model lives on: μ₀ = 0.11 (+0.45/−0.54) from DESI FS+BAO+BBN+nₛ,
  tightening to μ₀ = 0.05 ± 0.22 with CMB+DES-SN. A ~12% fσ₈ suppression is μ₀ of order
  −0.2 to −0.4 (exact mapping needs the time-dependence assumption worked out) — **inside
  DESI-alone's 1σ band**. This is the analysis a referee asks for first, and it gives a
  *weaker* verdict than the site's σ₈ shortcut.

This is **not** the reflexivity mechanism (absence of a primary source). Both statistics
here are real, cited, and independently checkable. The failure is a **protocol
mismatch**: the criterion was registered on one variable, the verdict was adjudicated on
a different variable, and nobody checked that the substitution was licensed. Name this a
third failure mode, distinct from *fabricated-refutation-laundering* (S63, invented data)
and *retraction-survival* (CDM verdict, citing a source that later reversed itself):
**criterion-verdict substitution** — the registered falsification statistic and the
delivered verdict statistic silently diverge, and the divergence survives every prior
citation-walk because both individual numbers are independently real.

## A second, mirror-image instance: over-promotion lives in explanatory prose, not verdicts

Pass 3 (Grad Student) independently found the opposite bias in a different textual
register. `core-idea/page.tsx:70,196-197`: *"Critical exponents (β, ν) are off by ~2×
— the diagnostic result that rules out C(ρ) as a Landau-theory continuum order
parameter."*

C(ρ) = [(1+x)^{2γ}−1]/[(1+x)^{2γ}+1] is real-analytic on x ∈ (−1, ∞): no singularity, no
free-energy functional, and **no length scale anywhere in the equation**, so a
correlation-length exponent ν is not "off by 2×" — it is undefined by construction. The
"2×" traces to β_eff = 1 (from the small-x expansion tanh(γ ln(1+x)) ≈ γx) versus
mean-field β = 1/2 — but β_eff = 1 is a tautology for *any* analytic function vanishing at
the origin, not a measurement. **"Off by 2×" implies the framework got within a factor of
two of being a Landau theory; it is not that kind of object at all — the gap is not 2×,
it is infinite.** The correct sentence ("no critical point, therefore no critical
exponents; asking for β and ν is a category error, not a near miss") is *more* damning
and *shorter* than what is on the page — over-promotion here does not even buy comfort, it
buys a weaker sentence.

The same register produced `/key-claims:196-197`, where a Session #581 quote uses γ for
what `/honest-assessment` correctly calls **B** (the boost ratio g_obs/g_bar). γ is
defined site-wide as 2/√N_corr, bounded above by 2; "⟨γ⟩ = 10.82" is impossible under the
site's own definition. This one doesn't make the framework sound more sophisticated — but
it **buries** the site's cleanest parameter-independent result (C ≤ 1 bounds the
quadrature boost; the framework cannot reach the observed deep-MOND discrepancies for
*any* parameters) under a symbol collision, which has the same net effect: a duller
sentence standing in for a sharper one.

## Why this is a distinct mechanism from reflexivity, worth tracking separately

Reflexivity (07-09) explains *verdict* statistics that lack a primary source. It does not
predict or explain:

1. A verdict statistic *with* a primary source breaking because the criterion and the
   delivered number are different variables (TEST-04a).
2. Explanatory/motivational prose — not verdicts at all — softening a result into a
   near-miss framing that isn't checked because it reads as scene-setting, not a claim
   (critical exponents; the Hill-function negative-cooperativity framing Pass 3 proposes
   as the correct replacement).

Put together with 07-09's reflexivity result, the site now has **three** independently
confirmed asymmetries, not one:

| mechanism | direction | why it survives citation-walks |
|---|---|---|
| Reflexivity (07-09) | over-refute / over-claim self | no primary source to walk to |
| Criterion-verdict substitution (new) | over-refute | both numbers real; nobody checks *which* number the criterion named |
| Explanatory-prose promotion (new) | over-promote (sounds more like physics) | reads as motivation/analogy, not a checkable claim |

The unifying thread Pass 4 names directly: **the self-critical statements on this site
are the least-audited statements on this site.** Nobody checks a refutation, because it
sounds humble, and humility reads as safe. It isn't — a false refutation is a false claim,
and (new this session) a *softened* explanatory sentence is a false claim in the other
direction, hiding in the part of the page nobody treats as a claim at all.

## Proposed mechanical rule (extends the 07-09 rule, does not replace it)

07-09 proposed: every self-referential number ships its regenerating script or is marked
self-reported. That rule doesn't catch either of today's instances — TEST-04a's two
numbers both already have scripts/citations; the critical-exponents sentence isn't a
number at all. Add two checks:

1. **Criterion-verdict variable match.** Every kill-criterion badge must state the
   variable the criterion was registered on *and* confirm the delivered verdict is
   adjudicated on that same variable, not a related-but-different one. If the site
   switches variables (as here, fσ₈ → σ₈), that switch must be justified in the same
   sentence, not silently absorbed into "Kill Criterion Triggered."
2. **Analogy/motivation sentences get the same citation-walk as verdicts.** Any sentence
   that compares the framework's math to a named physics structure (Landau theory,
   Ising model, mean-field scaling) is a checkable claim, not scene-setting, and should be
   walked exactly like a p-value: does the equation actually have the structure being
   compared to (a length scale, a free energy, a singularity), or does it only share a
   qualitative shape?

## Scope

This does not change any bucket in PREDICTIONS.md — TEST-04a's underlying disfavor and
the framework's overall "no field equation, no confirmed novel prediction" status are
unaffected; the local-density no-go and the absence of a growth-sourcing mechanism still
kill the framework's cosmology on independent grounds. What changes is *which* sentence
carries the kill, and how confidently. The honest verdict, once the substitution is
fixed, is closer to "the test as registered never had the statistical power to
discriminate" than to "kill criterion triggered" — which is a *more* useful negative
result for a referee, not a weaker one, for the same reason the critical-exponents fix
above is stronger once correctly stated.
