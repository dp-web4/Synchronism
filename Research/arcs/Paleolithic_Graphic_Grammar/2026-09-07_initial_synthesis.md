# Paleolithic Graphic Grammar - Initial Synthesis and Computational Probe

**Date:** 2026-09-07  
**Arc:** [Paleolithic Graphic Grammar](README.md)  
**Status:** First-pass synthesis, active  
**Posture:** Structure first, semantics later

## Executive summary

The initial question was whether recurring Upper Paleolithic cave symbols might contain something legitimately describable as **grammar**, and whether Jean-Jacques D.'s STAR categories - State, Thing, Action, Relationship - might capture it.

The answer after a first literature and quantitative pass is asymmetric:

- The case for **non-random, rule-governed compositional structure** is substantially stronger than expected.
- There is a decades-long literature explicitly treating Paleolithic graphic systems in syntactic and formal-grammar terms.
- Independent modern work using open network, multivariate, information-theoretic, and cognitive-science methods converges on structured composition without requiring decipherment.
- The most robust primitive operations appear to be variants of **repetition, concatenation/juxtaposition, superposition, integration/fusion, and recursive embedding**.
- The data do **not** currently support a privileged four-way STAR decomposition. Four classes can be imposed, and four clusters appear in some analyses for unrelated reasons, but the number four is not emerging as a stable semantic invariant.
- A more promising model is a **typed scene/composition grammar** in which figurative themes act as anchors, geometric marks can act as modifiers/operators/relations/quantities, and spatial/compositional relationships participate in the message.

The strongest working hypothesis from this pass is:

> **Some Paleolithic graphic systems are best approached as structured programs or typed visual compositions, not as strings of proto-words waiting for one-to-one translation.**

This is not a claim that we have decoded them.

---

## 1. The triggering hypothesis: STAR

The exploration began from a public post describing a STAR quadrant:

- **State**
- **Thing**
- **Action**
- **Relationship**

and an attempt to classify a recurring inventory of European Paleolithic cave signs into those four roles.

The immediate attraction was obvious: if symbols that differ in shape nevertheless occupy stable functional roles across caves and generations, a grammatical interpretation becomes plausible.

The immediate methodological problem was equally obvious:

> **Shape is not role.**

A zigzag is not evidence of an Action role merely because it looks dynamic. A circle is not evidence of a State role because it looks closed or static. If STAR is real, its categories should be recoverable from **behavior in context**:

- co-occurrence;
- substitution;
- compatibility/incompatibility;
- relative placement;
- repetition patterns;
- association with figurative themes;
- distribution across site, period, and context;
- held-out predictive value.

Therefore this arc treats STAR as a **candidate model to test**, not as the labeling scheme used to construct the evidence.

---

## 2. The older literature already asked the grammar question

### 2.1 Sauvet, Sauvet & Wlodarczyk (1977)

The most important historical find in this exploration is that the syntactic question is not new.

Georges Sauvet, Suzanne Sauvet, and André Wlodarczyk published *Essai de sémiologie préhistorique (Pour une théorie des premiers signes graphiques de l'homme)* in 1977. Rather than beginning with ethnographic translation, they treated the signs as a semiological system and analyzed their morphological decomposition and combination.

Their work identified recurrent construction relations including forms corresponding to:

- **juxtaposition**;
- **superposition**;
- **integration/fusion**.

The important point is methodological: they attempted to characterize **how signs combine before claiming what signs mean**.

That is exactly the epistemic order this arc adopts.

### 2.2 Sauvet & Wlodarczyk (2008)

The 2008 paper *Towards a Formal Grammar of the European Palaeolithic Cave Art* is even stronger.

Their corpus contained:

- **416 polythematic panels**;
- combinations of **2 to 6 themes**;
- **14 principal figurative motifs**;
- a statistical reduction to **five thematic classes**.

Only a small fraction of mathematically possible thematic combinations actually occurred. A simple system of rewriting rules accounted for **98% of the observed figurative compositions**.

That is formal grammar in a literal technical sense: a compact set of production constraints generates nearly the whole observed composition space.

It does **not** prove spoken-language-like syntax or reveal semantic translations. It does show that "these are just independent pictures on walls" is a poor null model.

---

## 3. Cognitive-science convergence: geometry may itself have a grammar

In a 2024 Collège de France lecture, Stanislas Dehaene explicitly revisited Paleolithic geometric signs and linked the archaeological literature to his laboratory's work on a human **language of geometry**.

His proposed elementary vocabulary is built from primitives such as:

- point;
- line;
- curve;

combined by operations such as:

- **repetition, with or without variation**;
- **concatenation**;
- **recursive embedding**.

The convergence with the Sauvet construction operations is striking because it comes from a very different direction: modern cognitive experiments and minimum-description-length representations rather than archaeological semiotics.

This suggests an important distinction:

> The grammar may begin **inside the glyph itself**, before any glyph-to-glyph "sentence" grammar is considered.

A complex Paleolithic sign may be a small generative program rather than an atomic token.

For example, instead of representing a motif as an indivisible label `SIGN_17`, a generative description might look like:

```text
primitive(line)
repeat(line, 4, parallel)
embed(curve, repeated_lines)
rotate(result, theta)
```

The relevant research question then becomes whether this program-like description predicts where and how the sign is used better than a purely visual category label.

---

## 4. Modern spatial and network evidence

### 4.1 Intxaurbe, Garate & Arriolabengoa (2024)

The 2024 open study *Drawing in the depths: spatial organization patterns related to Magdalenian cave art* analyzed **500 graphic units in nine caves** using GIS, iconographic variables, spatial variables, Factor Analysis for Mixed Data (FAMD), and Hierarchical Clustering on Principal Components (HCPC).

The initial analysis produced **four clusters**, but this is not evidence for STAR. One cluster consisted of nonfigurative representations. When those elements were removed and the figurative subset was reanalyzed, the system produced **three clusters**.

The substantive result was instead that cave art was strongly structured by context:

- some figures were placed and rendered to facilitate visibility and public comprehension;
- others were difficult to access or see;
- technique, completion, accessibility, cave depth, and spatial placement covaried.

This matters because it broadens the potential grammar. The physical support and location may function as part of the communication system.

A panel is not necessarily a sentence written on a neutral page. The cave itself may participate in the syntax.

### 4.2 Intxaurbe (2026)

The 2026 open paper *Mapping the Symbolic Structure of Palaeolithic Rock Art Using Co-occurrence Network Analysis* uses the same broad Magdalenian corpus to construct:

- frequency-weighted theme co-occurrence networks;
- Jaccard-normalized networks;
- filtered networks;
- minimum spanning trees;
- panel-theme bipartite models;
- orientation and inclination tests.

The key finding is stable **non-random, hierarchical, modular organization**. Large herbivores, particularly **bison, horse, and ibex**, occupy recurrent central positions. Formal variables such as orientation and inclination also show structured distributions.

Crucially, the paper is explicit that these methods operate at a **syntactic/structural level** and do not directly recover social meaning or symbolic intention.

That is exactly the distinction this arc needs:

```text
syntax / structure: measurable
semantics / translation: underdetermined
```

The associated code and selected data are open, making this the best current substrate for a reproducible follow-up experiment.

---

## 5. Information-theoretic evidence from much earlier mobile signs

Bentz & Dutkiewicz (2026), in *Humans 40,000 y ago developed a system of conventional signs*, analyzed several thousand engraved signs on more than 200 Aurignacian mobile objects dating roughly **43,000 to 34,000 years ago**.

Their results are useful precisely because they are conservative:

- the sign sequences are statistically distinguishable from modern writing;
- some statistical properties are comparable to early protocuneiform;
- information density varies systematically by object type;
- figurines carry higher information density than tools;
- the signs were used deliberately, systematically, and conventionally.

The authors explicitly stop short of proving numero-ideographic semantics or writing. They also note that the Aurignacian sequences lack some hallmarks of full writing, including strong productive combinatoriality and the rebus principle.

This is an independent evidence stream. It should **not** be merged naively with Magdalenian cave panels, but it establishes that conventional geometric sign systems in Europe are very old and can be studied quantitatively without translation.

The underlying SignBase project, published in 2020, provides a broader open catalog of geometric signs on mobile Paleolithic objects.

---

## 6. The Bacon et al. Y-sign hypothesis: useful even if the decipherment fails

Bacon et al. (2023) proposed a specific semantic reading for sequences of dots/lines and the Y sign associated with animal figures.

Their corpus contained:

- **606 sequences without Y**;
- **256 sequences with Y**;
- **862 total sequences**.

Published Table 1 gives:

| Group | Without Y | With Y | Total | Y rate |
|---|---:|---:|---:|---:|
| Aurochs | 30 | 15 | 45 | 33.3% |
| Bison | 89 | 41 | 130 | 31.5% |
| Caprid | 41 | 17 | 58 | 29.3% |
| Cervid | 102 | 50 | 152 | 32.9% |
| Fish | 132 | 8 | 140 | 5.7% |
| Horse | 199 | 104 | 303 | 34.3% |
| Mammoth | 13 | 0 | 13 | 0% |
| Bird | 0 | 21 | 21 | 100% |
| **Total** | **606** | **256** | **862** | **29.7%** |

The paper argues that dot/line counts encode lunar-month information and that Y marks parturition/birth timing.

That semantic interpretation is controversial. García-Bustos, Rivero, Sauvet & García Bustos (2023) identify methodological problems including questionable tracings and disputed associations between signs and animal figures. Other critiques question the statistical and cultural-calendar assumptions.

For this arc, the useful move is to **remove the proposed translation** and ask a more primitive question:

> Does Y behave like a reusable context-sensitive element rather than arbitrary decoration?

### 6.1 Aggregate heterogeneity test

Using only Bacon et al.'s published Table 1 counts, a chi-square test of Y presence against the eight analytical groups gives:

```text
chi-square = 98.123
_df = 7
p = 2.63e-18
Cramer's V = 0.337
N = 862
```

So Y occurrence is strongly non-random across the published animal groups.

This by itself says nothing about "birth". It says the sign is distributed contextually in the constructed corpus.

### 6.2 A surprising invariance among major terrestrial herbivores

Restrict the test to:

- aurochs;
- bison;
- caprid;
- cervid;
- horse.

Their Y-use rates lie in a narrow interval from **29.3% to 34.3%**.

A heterogeneity test gives:

```text
chi-square = 0.726
_df = 4
p = 0.948
```

Within this subset, there is essentially no detectable taxon-specific difference in whether a sequence contains Y.

That pattern is compatible with a reusable grammatical/operator-like element: the referent changes while the probability of using the element remains stable across a broad class of terrestrial herbivore anchors.

It is not proof of grammatical function. Corpus selection could also produce it.

### 6.3 Exhaustive latent-class partition probe

A deliberately simple model-selection experiment was then performed on the eight published groups.

Treat every sequence as Bernoulli `Y present / absent`. Allow taxa to share a latent Y-use probability. Enumerate **all 4,140 set partitions** of the eight taxa. For every partition:

1. estimate one maximum-likelihood Y probability per latent class;
2. compute binomial log likelihood;
3. score with BIC using the number of class probabilities as model complexity.

The best BIC solution has **three classes**:

```text
Class A:
  aurochs, bison, caprid, cervid, horse
  Y = 227 / 688 = 33.0%

Class B:
  fish, mammoth
  Y = 8 / 153 = 5.2%

Class C:
  bird
  Y = 21 / 21 = 100%
```

BIC comparison:

```text
3-class best partition: 955.64
8 independent taxon rates: 987.24
1 universal rate: 1055.45
```

The point is not that "three is the true grammar." Bird and mammoth sample sizes are small, and the corpus itself is disputed. The point is narrower:

> In this aggregate dataset, the information prefers a small number of behavioral classes over both one universal class and eight unrelated species classes - and it does not naturally select four.

This is evidence against treating the number four as privileged before analysis.

---

## 7. First STAR challenge

The Y example exposes the main problem with shape-first STAR assignment.

In the triggering STAR diagram, Y was placed under **Relationship**.

But under Bacon et al.'s proposed interpretation, Y would be closer to an **Action/event predicate**: birth / giving birth.

If instead Y encoded pregnancy, season, readiness, category membership, or something else, it could look more like **State** or **Relationship**.

The same physical glyph is compatible with several semantic roles. Its shape does not decide among them.

Therefore the correct STAR test is:

```text
1. infer behavioral/latent classes without STAR labels
2. characterize their distributional roles
3. only then ask whether State/Thing/Action/Relationship
   is a low-loss interpretation of those roles
```

Not:

```text
1. label signs as STAR from appearance
2. observe differences between the labeled groups
3. claim that STAR was discovered
```

The latter is circular.

---

## 8. First formal-feature clustering probe

A second exploratory probe used published/open orientation and inclination behavior from the Intxaurbe work as a crude behavioral representation of themes.

Each theme was represented by formal tendencies such as:

```text
left/right orientation bias
horizontal inclination tendency
inclined tendency
vertical tendency
inverted tendency
```

A preliminary standardized clustering sweep from `k=2..8` did **not** show a special preference for `k=4`. In the first reconstruction, K-means silhouette peaked near **k=5**, while Ward-style grouping was similarly competitive around **k=5..6**. Four-way grouping was possible but not privileged.

This result is intentionally tagged **PROVISIONAL / REPRODUCTION REQUIRED** because the first pass reconstructed from published/open summary variables rather than preserving a committed raw-data analysis script in this repository.

It is directionally consistent with the stronger published 2024 result:

- four HCPC clusters with nonfigurative units included;
- three clusters after removing the nonfigurative cluster.

The important conclusion is only that **"four" has not independently emerged as a stable number of semantic roles**.

A reproducible raw-data replication is in the next-stage plan.

---

## 9. What grammar currently seems to fit the evidence

The evidence increasingly favors a grammar operating at several nested levels.

### 9.1 Glyph-internal generative grammar

```text
SHAPE :=
    primitive
  | repeat(SHAPE, n, variation?)
  | concatenate(SHAPE, SHAPE)
  | superpose(SHAPE, SHAPE)
  | embed(SHAPE, SHAPE)
  | transform(SHAPE, orientation/scale)
```

This is close to the convergence between Sauvet's sign-combination work and Dehaene's geometry language.

### 9.2 Panel-level scene grammar

```text
COMPOSITION :=
    ANCHOR
  + ARTICULATION*
  + SPATIAL_RELATIONS
  + CONTEXT

ANCHOR :=
    dominant figurative theme
  | complex sign
  | repeated thematic center

ARTICULATION :=
    sign
  | quantity/repetition
  | secondary theme
  | modifier/operator

SPATIAL_RELATION :=
    juxtapose
  | superpose
  | integrate/embed
  | orient
  | contain
  | attach
```

### 9.3 Physical-context grammar

```text
CONTEXT :=
    cave
  + sector
  + panel
  + visibility
  + accessibility
  + support geometry
  + technique
  + chronology / regional convention
```

The physical context may constrain which compositions are valid, visible, restricted, public, or meaningful.

The result is less like:

```text
WORD WORD WORD WORD
```

and more like:

```text
TYPED GRAPH + GENERATIVE CONSTRUCTION + PHYSICAL CONTEXT
```

---

## 10. Why "program" is a useful metaphor

A Paleolithic graphic unit might encode meaning partly through **how it is constructed**, not just which named glyph it resembles.

For example:

```text
repeat(line, 5)
```

can carry a different informational affordance from:

```text
embed(curve, repeat(line, 5))
```

without requiring either construction to map one-to-one onto a spoken word.

At the panel level:

```text
anchor(horse)
attach(repeat(dot, 3), horse)
orient(horse, left)
place(panel, deep_restricted_sector)
```

is a structured object even before any English semantics are assigned.

This framing has several advantages:

1. it explains why spatial arrangement matters;
2. it allows recursion and composition;
3. it accommodates quantities and modifiers naturally;
4. it permits stable syntax with unknown semantics;
5. it offers a direct minimum-description-length test;
6. it makes held-out prediction possible.

---

## 11. What the evidence does NOT currently justify

### 11.1 Not a decipherment

No result here establishes that a given sign means "birth," "hunt," "water," "female," "movement," or any other English concept.

### 11.2 Not proof of writing

Bentz & Dutkiewicz explicitly distinguish their Aurignacian sign sequences from modern writing. The Paleolithic systems may be conventional external memory, notation, visual grammar, ritual composition, mnemonic systems, or several different things across time and region.

### 11.3 Not one pan-European language

The evidence spans enormous time and geography. Similar signs may be inherited, reinvented, repurposed, or semantically shifted. Regional and chronological structure must remain explicit.

### 11.4 Not evidence for Synchronism

This arc currently has no evidentiary role in the Synchronism physics framework.

### 11.5 Not evidence for STAR yet

The evidence supports structured composition much more strongly than it supports STAR's exact semantic partition.

---

## 12. Current evidence hierarchy

### Strong

- repeated, conventional sign inventories exist;
- motif/theme placement is non-random;
- panel composition obeys restricted combination rules;
- a small formal grammar can account for a very large fraction of figurative panel combinations in Sauvet & Wlodarczyk's corpus;
- spatial context and production choices are structured;
- open modern network analyses recover hierarchical/modular organization;
- early mobile sign systems show systematic information-density differences by artifact type.

### Moderate

- some geometric marks function as reusable modifiers/operators rather than independent referents;
- construction operations such as repetition/concatenation/embedding are cognitively fundamental and archaeologically relevant;
- typed scene grammar is a better modeling family than linear token sequence alone.

### Weak / unresolved

- exact semantic roles of specific signs;
- whether STAR's State/Thing/Action/Relationship categories correspond to stable latent roles;
- whether apparently similar sign systems across tens of millennia preserve semantics;
- whether Bacon et al.'s Y translation survives stronger corpus auditing.

---

## 13. Updated confidence after the first pass

These are working research calibrations, not literature-derived probabilities:

- **Real formal/compositional syntax in at least some Paleolithic graphic traditions:** >90%.
- **Some signs serving operator/modifier-like roles rather than only standalone referents:** roughly 50-70%.
- **Typed scene/composition grammar as a productive modeling frame:** roughly 70%.
- **STAR as the major four latent semantic roles:** roughly 10-15% at present.
- **Current English-language translations of individual signs:** very low confidence unless independently grounded.

The important movement from the start of the exploration is that confidence in **grammar** went up while confidence in **STAR specifically** went down.

That is a healthy result.

---

## 14. The next experiment that matters

The highest-value next step is not another interpretive essay. It is a predictive computational test:

> **Can a compact generative grammar of sign construction and scene composition compress and predict held-out Paleolithic graphic data better than surface-shape categories, raw co-occurrence, and STAR?**

The detailed protocol is in [`2026-09-07_method_and_next_experiments.md`](2026-09-07_method_and_next_experiments.md).

If the learned grammar wants four STAR-like roles, that will be interesting because we did not put them there.

If it wants three, six, continuous latent factors, or no stable role system at all, that is equally useful.

---

## Bottom line

The first pass finds a much more substantial pin than expected:

> **Upper Paleolithic graphic systems show multiple independent signatures of constrained, conventional, hierarchical composition. The defensible frontier is no longer simply "were the marks meaningful?" but "what formal grammar best describes their construction and use?"**

The signs may not be proto-words. They may be compositional programs embedded in scenes.

That possibility is worth a real arc.
