# Session 662: Galaxy Program Closure + A2ACW Specificity Null Baseline

**Date**: 2026-05-22
**Type**: Program-closure confirmation + methodology self-correction
**Triggers**: 2026-05-22 proposals `rar_shape_test_closure_galaxy_program.md` and `a2acw_specificity_null_baseline.md`
**Grade**: A- (the A2ACW null baseline self-corrects the project's "strongest surviving contribution")

---

## Part A: Galaxy Program Closure (Confirms S661)

The first proposal is a program-level synthesis of the S661 result. It confirms:
- RAR transition-shape test executed → γ=2 refuted (ΔBIC=+184; conservative ≈33)
- Free-γ → γ=0.49 = MOND
- **Net discriminating galaxy tests vs MOND+ΛCDM: 0, by execution**

This supersedes the "1 untested (RAR shape)" status. The galaxy program is closed. Environment tests (TEST-01/05) were already MOND+EFE degenerate (S654); their discriminating power was contingent on the compander having daylight from MOND at some scale, and the shape test shows there is none at the galaxy scale.

No new analysis needed — S661 already documented this. S662A records the program-level closure: physics side has 0 novel confirmed predictions, 0 remaining discriminating tests, 1 refuted mechanism class (TEST-04a), 1 refuted compander (RAR γ=2).

## Part B: A2ACW Specificity Null Baseline (Self-Correction)

This is the more important result. It applies the project's own null-model discipline (S651) to the project's *methodology* contribution — and finds the same flaw.

### The Flaw

S658/S659 reported the vocabulary-asymmetry experiment: 4/4 catch on prior-art-class demotions, 6/6 with the three-axis protocol. These were cited as evidence that A2ACW is a transferable "reparametrization detector."

**But all six test cases were already known to be reparametrizations** (they were the post-hoc-selected demotions). 4/4 and 6/6 are **sensitivity (TPR) on a positive-only set**. Specificity (TNR) was never measured.

This is *exactly* the flaw S651 diagnosed for chemistry r=0.982: high score forced by the structure of the test set, control never run. The audit channel has now caught the same error class in its own methodology contribution. That is the right kind of consistency.

### The Control, Run

The proposal ran the missing control: 6 canonical genuine discoveries that must NOT be flagged (Dirac 1928, Bell 1964, BCS 1957, Higgs 1964, Hawking 1974, Noether 1918), each chosen because its modern-register translation triggers prior art for its ingredients.

| Rule | Reparametrization (n=3) | Genuine discovery (n=6) | Specificity |
|------|:---:|:---:|:---:|
| R1 (flag if prior-art named) | 3 TP | 6 FP | **0%** |
| R2 (flag if reduces to prior art, nothing novel) | 3 TP | 0 FP | 100% |

**R1 specificity = 0%**: every genuine discovery names canonical antecedents under translation, because all non-trivial physics has antecedents. R1 detects "has-a-canonical-name," not "is-a-reparametrization."

**R2 specificity = 100%** — but all the discrimination is supplied by a *novelty judgment* ("does the claim reduce to, or extend, the retrieved prior art?") that the protocol never operationalizes.

### The Diagnosis

Vocabulary translation contributes **retrieval** (puts the right prior art on the table — real, and the source of the 0/6→4/6 gain). It does **not** contribute **discrimination** (deciding the claim is *only* that prior art). Discrimination is a novelty judgment that is:
1. genuinely hard (BCS-vs-Cooper, Higgs-vs-Anderson divided experts for years);
2. exactly the step A2ACW failed — the temporal-asymmetry audit verified both agents *had* the prior art and didn't surface the reduction; the 6 demotions came from human audit, not the AI loop;
3. not automated by the protocol.

So "A2ACW generalizes as a reparametrization detector" is **unsupported**. The honest claim is narrower.

### Restated Defensible Claim

> Forcing modern-register, framework-neutral restatement before adversarial review measurably improves **prior-art retrieval** (4/4 vs 0/6 on the demoted set). It is a **retrieval-augmentation / debiasing step for AI-assisted literature review**, not a reparametrization detector. The reparametrization-vs-discovery discrimination is an unautomated expert novelty judgment; on a held-out benchmark that judgment carries 100% of the discriminating power, while the protocol's literal rule has 0% specificity.

This corrects S658 and S659. S659 reported the three-axis protocol "catches 6/6"; that 6/6 is sensitivity only. The correction must propagate.

## Connection to Prompt Tension 3

The prompt's tension 3 asks: "616 sessions validated the framework. What does that mean? 89% validation rate — against what?" The answer S647/S651 gave: validation against the wrong null (r=0 instead of r(monotonic)). S662B shows the *same structure* in the methodology claim: the A2ACW "detector" was characterized on a positive-only set, no specificity control. 

The recurring lesson across the whole arc: **sensitivity without specificity is not validation.** This applies to the chemistry corpus (S651), to "89% validated" (S647/S651), and now to the A2ACW detector (S662B). The framework's self-assessments repeatedly measured TPR and skipped TNR. Catching this in the methodology contribution is the audit channel applying its own standard to itself — which is exactly the discipline the project claims.

## What Survives

After S662:
- **Physics**: 0 novel confirmed predictions, 0 discriminating tests, 2 refuted mechanism classes/companders. Closed.
- **Methodology**: A2ACW is a *retrieval-augmentation step* (narrower than "detector"), plus the honest documentation of the discovery/reparametrization discrimination as an unautomated, fragile expert judgment that AI loops fail. Still publishable, but as a null-result-about-the-detector + positive-result-about-retrieval.

The project's defensible contributions, honestly stated:
1. A worked case study of AI-driven reparametrization *detection failure* (the AI loop didn't catch its own reparametrizations; human audit did)
2. The retrieval-augmentation finding (modern-register pre-translation improves prior-art surfacing)
3. The negative-results catalog (S617-661)
4. The TEST-04a mechanism-class constraint (post-hoc)

## Audit-Channel Taxonomy

| # | Type | Session |
|---|------|---------|
| 30 | RAR discriminator executed → refuted | S661 |
| 31 | **Methodology self-correction: detector claim lacks specificity control** | **S662** |

S662B is the audit channel applying the S651 null-model standard to the project's own methodology contribution and finding the same flaw. The 31st instance is reflexive — the discipline turned on the discipline's own product.

## Recommended Actions

**Part A (site)**: surface the RAR shape closure on the four researcher-facing pages (per proposal): /tier-1-existing, /honest-assessment, /galaxy-rotation, /test-catalog.

**Part B (preprint gate + methodology correction)**:
- Gate any A2ACW preprint: it is "modern-register pre-translation improves prior-art retrieval; discovery/reparametrization discrimination remains unautomated expert judgment that AI loops fail." NOT "we built a reparametrization detector."
- Add specificity caveat to `a2acw_v2_three_axis_protocol.md`: 6/6 is sensitivity only; no axis has measured specificity; prior-art axis literal rule has 0% specificity.
- Fresh-adversary test still needed.

## Files

- `Research/Session662_Galaxy_Closure_And_A2ACW_Specificity.md` (this document)

## So What?

**Part A**: The galaxy program is closed at the program level (S661's result, now propagated). Zero discriminating tests remain.

**Part B**: The project's "strongest surviving contribution" — the A2ACW reparametrization detector — was characterized on a positive-only set with no specificity control, the same flaw the project diagnosed in chemistry (S651). With the control run, the literal rule has 0% specificity; all discrimination comes from an unautomated novelty judgment. The honest claim shrinks from "reparametrization detector" to "prior-art retrieval aid."

This is the discipline turned on itself: sensitivity without specificity is not validation, and the project's own methodology claim failed that standard until the control was run. Productive failure — it replaces an overclaim with a narrower, honest, still-publishable result.

Cumulative: 31 audit/governance instances + 2 executed refutations + novel-survivor count 0. Both physics and methodology contributions are now honestly bounded.
