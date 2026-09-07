# Paleolithic Graphic Grammar - Source and Data Ledger

**Arc:** [Paleolithic Graphic Grammar](README.md)  
**Last updated:** 2026-09-07

This ledger separates evidence streams that are easy to conflate: recurring sign inventories, sign-construction syntax, panel-level thematic grammar, cave spatial organization, proposed semantic decipherments, and early mobile sign systems.

The arc should prefer primary papers and open data/code wherever available.

---

## A. Structural/semiotic analysis of Paleolithic signs

### Sauvet, Georges; Sauvet, Suzanne; Wlodarczyk, André (1977)

**Title:** *Essai de sémiologie préhistorique (Pour une théorie des premiers signes graphiques de l'homme)*  
**Journal:** Bulletin de la Société préhistorique française, 74(2), 545-558  
**DOI:** https://doi.org/10.3406/bspf.1977.8467  
**Open record/full text:** https://www.persee.fr/doc/bspf_0249-7638_1977_hos_74_2_8467

**Why it matters:** Early explicit attempt to analyze Paleolithic signs as a semiological system without relying primarily on object-resemblance interpretations. Important for sign decomposition, compatibility, and compositional operations such as juxtaposition, superposition, and integration/fusion.

**Arc use:** Historical and methodological foundation for structure-before-semantics.

**Caveat:** Classification choices are from an earlier structuralist tradition and should not be treated as modern ground truth.

---

## B. Formal grammar of figurative panel composition

### Sauvet, Georges; Wlodarczyk, André (2008)

**Title:** *Towards a Formal Grammar of the European Palaeolithic Cave Art*  
**Journal:** Rock Art Research, 25(2), 165-172  
**DOI:** https://doi.org/10.69978/rar.v25i2.30  
**Article:** https://rockartresearch.com/index.php/rock/article/view/30

**Key published facts:**

- 416 polythematic panels;
- 14 main figurative motifs;
- factor analysis + hierarchical classification yielded five classes;
- only a small subset of possible combinations occurs;
- a few rewriting rules account for 98% of the observed figurative compositions.

**Why it matters:** Strongest direct precedent for using the term **formal grammar** on European Paleolithic cave-art composition.

**Arc use:** Benchmark for panel-level grammar and rewriting-rule approaches.

**Caveat:** The 98% figure is corpus/model specific and does not establish linguistic semantics.

---

## C. Cognitive-science convergence on geometric composition

### Dehaene, Stanislas (2024)

**Lecture:** *The origin of geometric symbols since prehistoric times: a language of thought?*  
**Venue:** Collège de France, 26 January 2024  
**Page/audio/support:** https://www.college-de-france.fr/en/agenda/lecture/the-perception-of-elementary-mathematical-objects-geometric-shapes-patterns-and-graphs/the-origin-of-geometric-symbols-since-prehistoric-times-language-of-thought

**Relevant proposal:** Human geometric cognition can be represented using elementary primitives such as points, lines, and curves, recombined through operations including repetition with variation, concatenation, and recursive embedding.

**Why it matters:** Independent convergence between cognitive-science models of geometric thought and archaeological sign-composition analyses.

**Arc use:** Basis for the proposed MDL/generative-geometry experiment.

**Caveat:** This is a lecture-level synthesis and cognitive hypothesis, not evidence that a specific Paleolithic motif carried a specific semantic meaning.

---

## D. Spatial organization of Magdalenian cave art

### Intxaurbe, Iñaki; Garate, Diego; Arriolabengoa, Martin (2024)

**Title:** *Drawing in the depths: spatial organization patterns related to Magdalenian cave art*  
**Journal:** Archaeological and Anthropological Sciences, 16, 104  
**DOI:** https://doi.org/10.1007/s12520-024-02007-3  
**Article:** https://link.springer.com/article/10.1007/s12520-024-02007-3  
**Data/code repository:** https://github.com/inakiintxaurbe/spatial-organization-patterns-related-to-magdalenian-cave-art

**Key published facts:**

- 500 graphic units;
- nine caves in the Cantabrian/Pyrenean region;
- GIS + iconographic/spatial variables + FAMD + HCPC;
- first analysis yielded four clusters;
- one cluster was nonfigurative;
- reanalysis excluding nonfigurative units yielded three clusters;
- visibility, accessibility, cave depth, technique, and completion show structured relationships.

**Why it matters:** Demonstrates that the physical cave setting and production strategy are structured variables, suggesting that a graphic grammar may extend beyond marks on a panel.

**Arc use:** Context-aware grammar and cluster-number sanity check.

**Caveat:** The clusters are functional/icono-topographic, not grammatical parts of speech.

---

## E. Reproducible co-occurrence/network analysis

### Intxaurbe, Iñaki (2026)

**Title:** *Mapping the Symbolic Structure of Palaeolithic Rock Art Using Co-occurrence Network Analysis*  
**Journal:** Journal of Archaeological Method and Theory, 33, article 59  
**Published:** 28 May 2026  
**DOI:** https://doi.org/10.1007/s10816-026-09796-y  
**Article:** https://link.springer.com/article/10.1007/s10816-026-09796-y  
**Code/graph repository:** https://github.com/inakiintxaurbe/Rock-Art-Theme-Co-occurrence-Network

**Methods:**

- panel-scale theme co-occurrence;
- frequency weights;
- Jaccard normalization;
- filtered networks;
- minimum spanning trees;
- panel-theme bipartite models;
- theme/orientation and inclination statistics.

**Key result:** Stable non-random hierarchical and modular organization, with large herbivores such as bison, horse, and ibex repeatedly occupying central positions.

**Why it matters:** Best current reproducible substrate for testing syntax without semantic assumptions.

**Arc use:** Primary dataset/code target for the next computational stage.

**Caveat:** Network centrality is not semantics. The author explicitly does not claim direct access to symbolic intention.

---

## F. Specific semantic/calendar hypothesis

### Bacon, Bennett; Khatiri, Azadeh; Palmer, James; Freeth, Tony; Pettitt, Paul; Kentridge, Robert (2023)

**Title:** *An Upper Palaeolithic Proto-writing System and Phenological Calendar*  
**Journal:** Cambridge Archaeological Journal, 33(3), 371-389  
**Published online:** 5 January 2023  
**DOI:** https://doi.org/10.1017/S0959774322000415  
**Article:** https://www.cambridge.org/core/journals/cambridge-archaeological-journal/article/an-upper-palaeolithic-protowriting-system-and-phenological-calendar/6F2AD8A705888F2226FE857840B4FE19

**Published corpus summary:**

- 606 dot/line sequences without Y;
- 256 sequences with Y;
- 862 total;
- largely France/Spain, spanning a long Upper Paleolithic interval.

**Authors' semantic proposal:** line/dot counts encode lunar-month information relative to a seasonal anchor; Y indicates parturition/birth timing.

**Arc use:** Contested but useful test case for whether a sign behaves contextually after stripping away the proposed translation. Published Table 1 supports the first aggregate Y-presence analysis in this arc.

**Caveat:** Treat as a **contested corpus and semantic hypothesis**, not established decipherment.

---

## G. Published critique of the Bacon corpus/interpretation

### García-Bustos, Miguel; Rivero, Olivia; Sauvet, Georges; García Bustos, Paula (2023)

**Title:** *Discussion: “An Upper Palaeolithic Proto-writing System and Phenological Calendar” by Bennett Bacon et al. (2023)*  
**Journal:** Journal of Paleolithic Archaeology, 6, article 32  
**Published:** 21 October 2023  
**Article:** https://link.springer.com/article/10.1007/s41982-023-00158-8

**Why it matters:** Documents methodological concerns including problematic tracings and uncertain associations between motifs/signs and animal depictions.

**Arc use:** Required uncertainty source for any Bacon-derived analysis.

**Rule:** No semantic conclusion may rely on Bacon-derived results without testing sensitivity to disputed examples and corpus construction choices.

---

## H. Early conventional sign systems and information theory

### Bentz, Christian; Dutkiewicz, Ewa (2026)

**Title:** *Humans 40,000 y ago developed a system of conventional signs*  
**Journal:** Proceedings of the National Academy of Sciences, 123(9), e2520385123  
**Published:** 23 February 2026  
**DOI:** https://doi.org/10.1073/pnas.2520385123  
**Open full text:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12956821/

**Corpus:** More than 200 mobile Aurignacian objects, roughly 43,000-34,000 years old, bearing several thousand geometric signs.

**Key result:** Sign sequences are clearly different from modern writing but display deliberate, systematic, conventional use; some information-theoretic properties are comparable to protocuneiform; information density differs by artifact class.

**Why it matters:** Strong independent evidence that very early geometric sign systems can be quantitatively structured without being full writing.

**Arc use:** Separate cross-corpus evidence stream and future information-theory comparison.

**Caveat:** Mobile Aurignacian signs are not interchangeable with later cave-art corpora.

---

## I. SignBase

### Dutkiewicz, Ewa; Russo, Gabriele; Lee, Saetbyul; Bentz, Christian et al. (2020)

**Title:** *SignBase, a collection of geometric signs on mobile objects in the Paleolithic*  
**Journal:** Scientific Data, 7, 364  
**Published:** 23 October 2020  
**DOI:** https://doi.org/10.1038/s41597-020-00704-x  
**Article:** https://www.nature.com/articles/s41597-020-00704-x  
**Open mirror:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7585433/

**Why it matters:** Open structured data source for geometric signs on mobile Paleolithic artifacts across multiple regions/periods.

**Arc use:** Candidate for cross-corpus geometric-complexity and role-transfer analysis.

---

## J. von Petzinger recurring-sign inventory

### Genevieve von Petzinger

**Relevant work:** Pan-European cataloging of recurring geometric signs in Upper Paleolithic cave art, commonly summarized as a repertoire of roughly 32 recurring sign types across a large geographic and chronological range.

**Arc use:** Motivating inventory and candidate target for generative-geometry encoding.

**Caveat:** The full research database is not presently treated as an openly downloadable analysis corpus in this arc. Public summaries should not be mistaken for complete occurrence-level data.

**Important distinction:** A recurring sign-type inventory demonstrates shape recurrence, not stable semantics or grammar by itself.

---

## K. Triggering STAR proposal

### Jean-Jacques D. (public post, observed 2026-09-07)

**Proposal:** STAR = State, Thing, Action, Relationship, applied as a four-quadrant proto-grammar and tentatively mapped onto recurring Paleolithic signs.

**Arc use:** Hypothesis generator only.

**Epistemic status:** No evidence from this source is used to validate STAR. The arc deliberately tests whether STAR emerges after behavior-based analysis.

---

# Evidence-stream separation

| Evidence stream | Main source(s) | What it can support | What it cannot support alone |
|---|---|---|---|
| Recurring sign shapes | von Petzinger / SignBase | Stable/repeated graphic inventories | Shared meaning |
| Sign construction | Sauvet 1977; Dehaene 2024 | Composition operators / geometry grammar | Specific semantics |
| Figurative panel combinations | Sauvet & Wlodarczyk 2008 | Formal rewriting grammar | Spoken-language syntax |
| Cave spatial context | Intxaurbe et al. 2024 | Structured placement/function | Glyph meaning |
| Theme networks | Intxaurbe 2026 | Hierarchy, modularity, co-occurrence syntax | Social intention |
| Dot/line/Y sequences | Bacon et al. 2023 | Testable sign/context association | Accepted calendar decipherment |
| Bacon critique | García-Bustos et al. 2023 | Corpus uncertainty / methodological limits | Absence of all structure |
| Aurignacian mobile sequences | Bentz & Dutkiewicz 2026 | Conventional/systematic external signs | Direct continuity with Magdalenian cave grammar |

---

# Data priority for next stage

1. **Intxaurbe 2026 repository** - highest priority because the workflow is explicitly reproducible and structurally focused.
2. **Intxaurbe et al. 2024 repository** - raw graphic-unit/context variables for spatial grammar and cluster replication.
3. **SignBase** - independent mobile-sign corpus for generative geometry and information-theory transfer.
4. **Bacon supplementary corpus** - only with a parallel uncertainty/critique layer.
5. **von Petzinger occurrence data** - pursue only if a legitimate occurrence-level dataset becomes publicly or directly available.

---

# Source-handling rule

Every future arc result should identify whether it is based on:

- **primary occurrence data**;
- **published aggregate table**;
- **author-coded category**;
- **our derived feature**;
- **our latent model output**;
- **semantic interpretation**.

These levels must never be collapsed into one another.
