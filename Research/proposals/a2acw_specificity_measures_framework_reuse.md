# Back-annotate from site: A2ACW Specificity Is Measurable-But-Mislabelled — It Scores Framework-Reuse, Not Novelty

**Filed:** 2026-06-27 (explorer track, synchronism-site)
**Source finding:** `synchronism-site/explorer/findings/a2acw-specificity-measures-framework-reuse-not-novelty.md`
**Builds on:** `a2acw_specificity_null_baseline.md` (2026-05-22, the R1/R2 result),
`a2acw_temporal_asymmetry_redesign.md`, `a2acw_vs_ai_discovery_landscape.md`
**Trigger:** 2026-06-27 visitor Pass 4 (leading-edge researcher): *"the specificity denominator
is empty"* + maintainer-seeded topic `a2acw-specificity-null-denominator.md`.

## What the predecessor established (and stands)

R1 (literal "names prior art") → 0% specificity; R2 (steelmanned "reduces to prior art,
nothing added") → 100% specificity but all discrimination supplied by an unautomated novelty
judgment. Diagnosis: retrieval ≠ discrimination. **This is not overturned.**

## What this adds

1. **"Empty denominator in principle" is false.** The contamination literature's standard
   time-based partition (items dated after the evaluator's training cutoff are provably
   uncontaminated; arXiv:2406.04244 and the Min-K%/ConStat/CoDeC family) constructs a clean
   specificity corpus. The 2026-05-22 corpus (Dirac…Noether) was merely the *maximally
   confounded* cell — canonical **and** memorized — so it could not separate "flagged for
   conceptual priority" from "flagged for memorization." That is a corpus defect, not an
   impossibility.

2. **A clean post-cutoff specimen still gets flagged by R1 — so 0% specificity is structural.**
   AlphaEvolve's 2025 result (4×4 complex matmul in 48 scalar multiplications, beating
   Strassen 1969, externally/automatically verified) is a genuine, provably-uncontaminated,
   post-cutoff discovery. Run through the protocol: **R1 flags it** (it names Strassen /
   tensor-decomposition one round in), **R2 passes it** (it *extends* — adds a new, verified,
   better object). AlphaEvolve is the clean specimen that **splits R1 from R2**, and shows R1's
   0% is not an artifact of using famous in-training examples.

3. **R2's smuggled judgment is the test-suite question.** R2 reaches opposite verdicts on
   AlphaEvolve (extends → pass) and Synchronism (reduces to MOND → flag). The discriminating
   question — "does it add a prediction/object the prior framework didn't have?" — is, for a
   physics theory, *exactly* what ΔBIC=+184 / DESI TEST-04a / the wide-binary feasibility kill
   answer. R2's 100% specificity is the test suite's authority wearing the protocol's costume.

4. **Therefore the metric measures the wrong axis.** A2ACW's automatable rule (R1) classifies
   **framework-reuse** ("new object inside an existing framework"), which describes the
   overwhelming majority of valuable science (Dirac, BCS, Higgs, AlphaEvolve). The true-negative
   class it would need to *pass* — framework-founding work — is near-empty across all of science
   and is precisely the class an in-distribution model cannot recognize anyway. So a cleanly
   measured specificity would still read ≈0%, not from method-blindness but from the construct
   being mislabelled. Fixing the corpus does not rescue the metric.

## The clean separation (use this in any write-up)

Synchronism's conviction rests on two properties the site/archive currently conflate:

| Property | Established by | Indictment strength |
|---|---|---|
| Framework-reuse (C(ρ) ∈ MOND/compander family) | A2ACW (R1) + human audit | **Weak** — shared with AlphaEvolve/Dirac/BCS; the norm for good science |
| Zero new predictions | the physics test suite | **Strong** — the actual indictment |

A2ACW supplies only row 1. Row 1 alone does not convict. The conviction is row 2, which A2ACW
had no part in producing. **Do not let the A2ACW null borrow the test suite's authority.**

## Restated defensible claim

> A2ACW's literal adversarial rule is a *framework-reuse classifier*: it flags any ansatz that
> sits inside an existing formal framework, which includes nearly all genuine discoveries
> (demonstrated on AlphaEvolve's verified post-cutoff matmul record). Its 0% "specificity" is
> structural and corpus-independent. What separates a worthless reparametrization from a
> valuable framework-reuse is whether the new object adds predictive content — an unautomated
> judgment that, for Synchronism, is delivered by the test suite, not the protocol.

## Decisive experiment (one session, paper-able)

Run AlphaEvolve's 48-mult (and 2–3 other post-cutoff verified results) through the documented
PRIMARY/CHALLENGER protocol under both R1 and R2, evaluator cutoff predating the result.
**Predicted: R1 flags, R2 passes (R2 pass requires a human "extends" assertion).** Disconfirming
outcome (R1 passes a genuine post-cutoff discovery) would partially rehabilitate automated
specificity and is the more publishable surprise. Second corpus: historical framework-founding
work (Noether, path integral, Galois) — prediction: R1 flags those too, confirming the
true-negative class is empty in practice.

## Actions

1. **Gate any A2ACW preprint on the separation.** The honest contribution is a *retrieval-recall*
   characterization (well-posed against ground-truth prior-art links), not a sensitivity/
   specificity classifier (ill-posed — "pass" is absence of retrieval from a non-self-auditable
   index). Replace the specificity row entirely with recall.
2. **Update `a2acw_v2_three_axis_protocol.md`:** the "specificity unmeasured" caveat upgrades to
   "specificity measurable via temporal holdout but mislabelled; R1 ≈0% is structural; the axis
   scored is framework-reuse."
3. **Site (`/research-philosophy`, `/honest-assessment`, `/for-researchers`):** keep the
   framework-reuse row (A2ACW) and the predictive-emptiness row (tests) separate; add the
   AlphaEvolve one-liner so the null is exactly as strong as it is and no stronger.

## Status

Productive sharpening. Converts "specificity is unmeasurable / denominator empty" (a dead end)
into "specificity is measurable and measures framework-reuse, the wrong axis" (a decidable,
publishable claim with a one-session confirming experiment). Consistent with — and tightening
of — the 2026-05-22 retrieval-not-detector result.
