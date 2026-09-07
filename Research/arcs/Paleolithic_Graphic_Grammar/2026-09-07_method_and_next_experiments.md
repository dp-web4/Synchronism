# Paleolithic Graphic Grammar - Method and Next Experiments

**Date:** 2026-09-07  
**Arc:** [Paleolithic Graphic Grammar](README.md)  
**Purpose:** Turn an intriguing structural pattern into falsifiable computational tests

## Research objective

The arc should answer a narrow scientific question:

> **What representational and grammatical model best predicts unseen Paleolithic graphic compositions without assuming their meaning in advance?**

The primary comparison is not "grammar versus no grammar." Several published analyses already show non-random structure. The sharper question is **which grammar**, at what scale, and with how much transfer across corpora.

The most important methodological rule is to prevent the model from rediscovering labels we supplied ourselves.

---

## Competing model families

At minimum, every serious experiment should compare these families.

### M0 - Frequency baseline

Predict motif/sign occurrence from marginal frequencies only.

This asks whether apparent structure disappears once common motifs are accounted for.

### M1 - Pairwise co-occurrence baseline

Predict from pairwise association statistics such as:

- Jaccard similarity;
- pointwise mutual information;
- conditional frequency;
- network centrality / neighborhood.

This captures ordinary association without grammar.

### M2 - Surface-shape categories

Represent signs by conventional morphology only:

- dots;
- lines;
- chevrons;
- quadrangles;
- claviforms;
- tectiforms;
- circles/ovals;
- branching/Y forms;
- etc.

This tests how much structure lives in visual category alone.

### M3 - Generative geometry / MDL grammar

Represent a motif as a program constructed from primitives and operators:

```text
PRIMITIVE := point | line | curve | closed_curve | junction

OP :=
    repeat
  | repeat_with_variation
  | concatenate
  | superpose
  | embed
  | reflect
  | rotate
  | scale
  | branch
```

Optimize for minimum description length or a related complexity criterion.

### M4 - Typed scene grammar

Model a panel as a typed relational graph:

```text
nodes:
  figurative themes
  signs
  repeated groups
  support/context nodes

edges:
  adjacent
  superposed
  embedded
  attached
  contained
  oriented_toward
  repeated_with
  same_panel
```

Learn latent node/edge roles from distribution rather than names.

### M5 - STAR-constrained grammar

Force each relevant element or latent role into one of:

- Thing
- Relationship
- State
- Action

STAR should be scored only on predictive/compressive performance and stability. It gets no bonus for conceptual neatness.

### M6 - Unconstrained latent-role model

Allow the number and nature of latent roles to be selected from the data.

Candidate approaches:

- stochastic block models;
- nonnegative matrix factorization;
- latent class models;
- Bayesian mixture models;
- graph embeddings followed by clustering;
- minimum-description-length community models;
- role discovery based on structural equivalence.

This is the key comparator for STAR. If M6 repeatedly converges on four stable roles whose behavior maps cleanly onto STAR, that is evidence. If it does not, STAR loses.

---

## Experiment 1 - Reproduce the published open datasets locally

### Goal

Establish a clean, versioned baseline before adding interpretation.

### Data

1. Intxaurbe, Garate & Arriolabengoa (2024) open 500-GU dataset.
2. Intxaurbe (2026) co-occurrence workflow and Gephi exports.
3. SignBase metadata where accessible under its published/open terms.
4. Bacon et al. aggregate and supplementary tables, retained as a **contested corpus** with explicit provenance.

### Deliverables

- immutable raw-data hashes;
- schema documentation;
- corpus provenance table;
- exact reproduction of published summary statistics where possible;
- a `KNOWN_ISSUES.md` containing every ambiguity that could alter sign/animal/panel association.

### Falsifier

If the published open analyses cannot be reproduced within reasonable numerical tolerance, stop all higher-level inference until the discrepancy is understood.

---

## Experiment 2 - Reproduce the first aggregate Y-sign probe

### Input

Bacon et al. Table 1:

| Group | No Y | Y |
|---|---:|---:|
| Aurochs | 30 | 15 |
| Bison | 89 | 41 |
| Caprid | 41 | 17 |
| Cervid | 102 | 50 |
| Fish | 132 | 8 |
| Horse | 199 | 104 |
| Mammoth | 13 | 0 |
| Bird | 0 | 21 |

### Tests

1. Chi-square Y-presence heterogeneity across groups.
2. Terrestrial-herbivore-only heterogeneity.
3. Exhaustive enumeration of all 4,140 set partitions of the eight groups.
4. Binomial maximum likelihood per latent class.
5. BIC, AIC, and leave-one-group-out predictive likelihood.
6. Bootstrap stability of the selected partition.

### Registered first-pass expectations

The initial calculation found:

```text
all groups:
  chi-square = 98.123
  df = 7
  p = 2.63e-18
  Cramer's V = 0.337

terrestrial herbivores only:
  chi-square = 0.726
  df = 4
  p = 0.948

best BIC partition:
  [aurochs, bison, caprid, cervid, horse] -> 227/688 = 33.0%
  [fish, mammoth]                        ->   8/153 = 5.2%
  [bird]                                 ->  21/21  = 100%

BIC:
  best 3-class = 955.64
  8 independent rates = 987.24
  1 universal rate = 1055.45
```

### Critical caveat

This tests the **constructed Bacon corpus**, not Paleolithic reality directly. Published criticism of tracing and association choices must be encoded as uncertainty or alternative corpus variants before semantic weight is placed on the result.

### Falsifier

If modest plausible corrections to disputed cases destroy the class structure, treat the aggregate pattern as corpus-sensitive rather than archaeological evidence.

---

## Experiment 3 - Geometry as program: MDL encoding of recurring signs

### Goal

Test Dehaene/Sauvet-style generative construction directly.

### Step A - Define a deliberately small primitive language

Start with:

```text
P = {
  point,
  line,
  curve,
  closed_curve,
  junction
}
```

and operators:

```text
O = {
  repeat(n),
  translate(dx,dy),
  rotate(theta),
  reflect(axis),
  scale(s),
  concatenate,
  superpose,
  embed,
  branch(k)
}
```

Do not add a primitive merely to make one difficult glyph cheap. Every new primitive pays a global model-complexity cost.

### Step B - Encode the recurring sign inventory

For each sign family, store:

- shortest known program;
- alternate near-shortest programs;
- program length;
- primitive count;
- operator count;
- recursion depth;
- symmetry;
- repetition structure;
- branching structure.

### Step C - Compare against surface features

Ask whether program features predict:

- site/region;
- chronology;
- co-occurring themes;
- panel role;
- cave depth/context;
- sign-sign compatibility;

better than raster/shape descriptors or conventional sign names.

### Strong pass

A compact generative description yields significant held-out predictive gain and transfers across sign families or sites.

### Kill

Program complexity is merely another way to describe visual similarity and adds no held-out predictive information.

---

## Experiment 4 - Panel completion as a grammar test

### Core task

Given an incomplete panel representation, predict the withheld element.

Examples:

```text
known: horse + ibex + [MASK]
predict: likely theme/sign class
```

or:

```text
known graph:
  horse --adjacent--> sign_A
  horse --orientation--> left
  sign_A --[MASK]--> sign_B

predict missing edge / node type
```

### Why this matters

Held-out completion converts "interesting association" into a measurable grammar problem.

### Split discipline

Use increasingly difficult splits:

1. random panel holdout;
2. cave holdout;
3. regional holdout;
4. chronological holdout where the data support it.

A grammar that only works under random holdout may be memorizing site conventions.

### Metrics

- negative log likelihood;
- top-k accuracy;
- mean reciprocal rank;
- Brier score / calibration;
- description length;
- transfer degradation by held-out cave/region.

---

## Experiment 5 - Latent role discovery

### Goal

Ask whether the corpus contains stable "parts of speech" without supplying their names.

### Feature families

For every motif/sign/theme:

- marginal frequency;
- panel degree;
- betweenness / hubness;
- entropy of co-occurring neighbors;
- substitution profile;
- cave/site dispersion;
- association with dominant figurative anchors;
- spatial relation distribution;
- orientation/inclination behavior;
- repeatability;
- visual/generative complexity.

### Analysis

Fit role models for a range of latent dimensions/classes.

Select complexity using:

- BIC / integrated completed likelihood;
- MDL;
- held-out likelihood;
- stability under bootstrap resampling;
- stability under cave/region holdout.

### STAR test

Only after roles are learned, ask independent annotators/models to characterize them behaviorally.

Then compare a four-role STAR interpretation against alternatives.

A serious STAR result would require all of the following:

1. approximately four roles emerge without being requested;
2. the roles are stable under resampling and site holdout;
3. mapping to Thing/Relationship/State/Action is low-ambiguity;
4. the mapping improves prediction or compression;
5. the same role structure appears in more than one corpus.

Anything less is suggestive at best.

---

## Experiment 6 - Is the grammar linear, graphical, or spatial?

Compare three representations of the same evidence:

### Sequence model

```text
A B C D
```

appropriate where marks form clear linear sequences.

### Graph model

```text
A --relation1--> B
A --relation2--> C
```

appropriate for panels with multiple spatial relations.

### Scene-plus-context model

```text
panel_graph + cave_sector + visibility + accessibility + support geometry
```

### Question

Which representation best predicts held-out content?

### Expected informative outcome

It is entirely possible that different Paleolithic systems use different grammars:

- mobile notational sequences may be more linear;
- cave panels may be graph-like;
- deep-cave compositions may require physical context.

The arc should not force one universal representation across all media.

---

## Experiment 7 - Corpus chronology and semantic drift

### Problem

The Upper Paleolithic spans tens of millennia. A recurring shape does not guarantee a recurring meaning.

### Test

For each sign family or latent role, estimate distributional embedding by period/region.

Measure drift:

```text
role_embedding(sign, period_1)
vs
role_embedding(sign, period_2)
```

### Pass

Some signs retain stable distributional roles over long periods while others drift. This would permit a nuanced notion of inherited grammar with changing lexicon/semantics.

### Kill for "single language" hypotheses

If distributional roles fragment strongly by period and region, reject pan-European continuity claims even if the shapes recur.

---

## Experiment 8 - Negative controls

A grammar pipeline can manufacture patterns from almost anything. Controls are mandatory.

### Control A - degree-preserving panel randomization

Shuffle graph edges while preserving motif frequencies/degrees.

### Control B - within-cave shuffle

Preserve site-level abundance while destroying panel composition.

### Control C - synthetic grammars

Generate corpora from known grammars and verify that the pipeline recovers the correct class count and operators.

### Control D - synthetic non-grammar

Generate frequency-matched independent motifs and verify that the pipeline does **not** hallucinate stable roles.

### Control E - analyst-blind labels

Run the role-discovery stage on anonymous motif IDs. No labels like `horse`, `Y`, `tectiform`, or `action-like` should be visible to the classifier.

---

## Experiment 9 - Bayesian treatment of disputed archaeological associations

For contested corpora, a binary "include/exclude" decision is too brittle.

Represent uncertain associations as probabilities:

```text
P(sign belongs to animal/panel) = q
```

Run posterior inference or Monte Carlo corpus realizations.

Report:

- findings stable across plausible annotation uncertainty;
- findings that depend on a handful of disputed cases.

This is especially important for Bacon-style sign/animal associations.

---

## Experiment 10 - Cross-corpus transfer

The most ambitious test asks whether a grammar learned in one domain helps another.

Candidate transfer directions:

```text
Magdalenian cave panels -> other Magdalenian regions
Aurignacian mobile signs -> later mobile signs
sign-internal geometry -> sign co-occurrence role
figurative panel grammar -> mixed figurative/nonfigurative panels
```

A grammar that transfers is much more interesting than one that merely compresses one curated dataset.

No transfer should be expected automatically. Failure to transfer may itself reveal multiple distinct graphic traditions.

---

## Suggested implementation order

### Stage 0 - Reproducibility foundation

- mirror/open-reference the public datasets;
- checksum versions;
- reproduce Intxaurbe network statistics;
- reproduce published HCPC/FAMD summaries where practical;
- reproduce the Bacon Table-1 probe in committed code.

### Stage 1 - Cheap latent structure

- co-occurrence baselines;
- latent role discovery;
- bootstrap stability;
- cave-holdout prediction.

### Stage 2 - Generative geometry

- encode a manageable subset of recurring geometric signs;
- build MDL grammar;
- test whether program features predict use.

### Stage 3 - Typed scene grammar

- panel graphs;
- masked-node/edge prediction;
- compare linear versus graph versus context-aware representations.

### Stage 4 - STAR adjudication

- map learned roles after the fact;
- compare constrained four-role model against unconstrained alternatives;
- publish a positive or negative verdict without moving the goalposts.

---

## What would count as a significant result?

### Result A - Grammar confirmed, semantics unknown

A compact compositional model robustly predicts held-out panel/sign structure across sites.

This alone would be important.

### Result B - Stable latent role system

A small set of distributional roles recurs across sites and periods.

This would justify talking about visual "parts of speech" in a careful structural sense.

### Result C - STAR survives

Four roles emerge independently, transfer, and map cleanly onto State/Thing/Action/Relationship.

This would be genuinely surprising and worth separate publication.

### Result D - STAR fails, another grammar wins

This is at least as interesting. The goal is discovery, not defense of the trigger hypothesis.

### Result E - No transferable grammar

Local panel rules exist, but there is no stable cross-site or cross-period grammar.

That would constrain stronger claims and still improve archaeological methodology.

---

## Current best hypothesis to try to falsify

> **Paleolithic graphic communication is heterogeneous but contains nested compositional grammars. Complex signs are generated from a small geometric operation set; panels are typed relational compositions organized around anchor themes; physical cave context constrains valid composition and use. Semantic translation is not required for this structure to be measurable.**

### Kill conditions

Kill or sharply narrow that hypothesis if:

1. held-out prediction gains disappear under degree/frequency controls;
2. learned roles are unstable across resampling;
3. grammar does not transfer even within tightly matched regional/chronological subsets;
4. MDL geometry adds no information beyond visual similarity;
5. scene-graph structure performs no better than pairwise frequency;
6. apparent structure is dominated by analyst classification choices.

---

## Publication discipline

Any future write-up should maintain three separate claim levels:

### Level 1 - Observed

Examples:

- motif A co-occurs with motif B above chance;
- a model predicts withheld themes with X accuracy;
- a latent class is stable under bootstrap.

### Level 2 - Structural interpretation

Examples:

- motif A behaves like an operator/modifier;
- the corpus is well-described by a recursive compositional grammar.

### Level 3 - Semantic interpretation

Examples:

- motif A means birth;
- motif B denotes a state;
- a panel communicates a specific proposition.

The evidentiary bar rises sharply at each level. Most of the arc should remain at Levels 1 and 2 until independent semantic anchors exist.

---

## Closing principle

The right question is not:

> "Can we make these signs fit our ontology?"

It is:

> **"If the signs had a grammar and we did not know its vocabulary, what measurable regularities would let the grammar reveal itself?"**

That is the experiment.
