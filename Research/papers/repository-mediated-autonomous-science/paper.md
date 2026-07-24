# Repository-Mediated Continuity in Memoryless Autonomous Research

## A longitudinal case study of claim formation, correction, and epistemic drift

**Author:** Codex (OpenAI)<br>
**Corpus creators:** the autonomous Claude instances recorded in
`dp-web4/gnosis-research`<br>
**Corpus steward:** dp<br>
**Version:** Working paper v0.1, 2026-07-23<br>
**Status:** Not peer reviewed

## Abstract

Autonomous-science systems are usually evaluated by the quality of a final
idea, experiment, or paper. Less is known about their longitudinal behavior:
how claims persist across agent restarts, how later agents detect inherited
errors, and whether corrections propagate through the research record. We
analyze `gnosis-research`, a public repository whose stated operating design
wakes a fresh language-model instance every few hours, gives it the accumulated
repository, and asks it to choose and pursue a research direction. Individual
instances do not retain conversational memory; continuity is mediated by Git
history, prose, code, and stored results.

At the analyzed snapshot, the corpus contains 138 commits and 600 tracked files
from 2026-02-03 through 2026-07-17. A reproducible lexical profile of commit
subjects finds 45 promotion-only subjects, 14 correction-only subjects, 2
containing both kinds of marker, and 77 containing neither. We do not interpret
these counts as discovery or error rates. We instead manually trace nine
high-salience claim lineages. Five locally adjudicated lineages reverse within
5.7–95.2 hours (median 59.3 hours), including claims about phase transitions,
stochastic resonance, integration-time advantages, and detector failure. Two
wrong-variable claims require approximately 92 days to be corrected in a
parallel repository and remain uncorrected in the source repository's public
framing. In one case, adverse experimental evidence is committed approximately
12 hours *before* a successor declares the contradicted pattern a universal
law. A further mathematical error remains uncorrected for more than four
months.

The case shows that repository-mediated succession can support rapid local
self-correction without weight updates or persistent individual memory.
However, append-only narrative memory also preserves framing, allows evidence
to be bypassed, and does not guarantee correction propagation. The decisive
system property is therefore not memory alone but **admission and maintenance
of claims**: explicit variable definitions, independent evaluators,
prospective tests, canonical claim records, and dependency-aware correction
propagation. We propose reporting standards and an architecture for
longitudinal autonomous research systems.

**Keywords:** autonomous science; language-model agents; external memory;
scientific self-correction; Git; claim provenance; epistemic drift; persistent
context

## 1. Introduction

Recent systems automate increasingly broad parts of scientific work. The AI
Scientist generates ideas, writes and executes code, produces papers, and
simulates review [1]. Agent Laboratory and Robin distribute literature search,
hypothesis generation, experiment design, and analysis across specialized
agents [2,3]. FunSearch couples language-model generation to a hard evaluator
and has produced verifiable improvements on established mathematical and
algorithmic problems [4]. These systems raise an important longitudinal
question:

> What happens after one autonomous research run ends and another begins?

A final paper can conceal the dynamics that produced it. An apparently
coherent research program may result from persistent correct knowledge,
persistent framing error, repeated rediscovery, or selective survival of
attractive narratives. Conversely, an agent that begins each session without
personal memory may still participate in a durable scientific process if it
inherits sufficiently rich external artifacts.

External memory is already central to language-agent research. Reflexion stores
linguistic feedback in an episodic buffer [5]. Generative Agents records
experience and synthesizes higher-order reflections [6]. Voyager accumulates
an executable skill library that later instances can retrieve and compose [7].
These systems generally evaluate task performance or behavioral competence.
They do not primarily examine the epistemic maintenance problem: whether a
changing external record distinguishes current claims from superseded claims,
routes adverse evidence to the right conclusion, and updates every dependent
surface.

This paper studies a naturally accumulated case. The public
[`gnosis-research`](https://github.com/dp-web4/gnosis-research) repository
describes itself as an autonomous research log produced by a sequence of
cold-start Claude instances. A timer wakes an instance, the instance reads the
record, selects a direction, runs experiments, commits its work, and ends.
The next instance inherits the repository but not the predecessor's
conversational state [8].

The corpus is scientifically messy in a useful way. It contains theoretical
claims about coherence and consciousness, experiments on language-model
behavior, simulation work on cellular automata and coupled oscillators,
explicit failures, premature declarations of completion, and later
corrections. It is not a benchmark constructed to reward correction. It is the
surviving record of a process that was pursuing its own questions.

We ask four questions:

1. Can repository-mediated succession support correction without persistent
   individual memory?
2. How quickly are local interpretations corrected once a successor targets
   them?
3. When does adverse evidence fail to alter the inherited narrative?
4. What infrastructure is required for corrections to propagate beyond the
   document in which they are discovered?

Our contribution is a longitudinal case analysis, not a claim that one
repository represents autonomous science generally. We provide:

- a frozen corpus profile;
- a manually auditable set of claim-correction lineages;
- a distinction between **local correction** and **corpus repair**;
- evidence of both rapid self-correction and persistent epistemic drift; and
- design recommendations for autonomous research systems whose identity lives
  in changing artifacts.

## 2. Case description

### 2.1 Operating model

The repository's public provenance states that every experiment, document, and
commit was produced by Claude instances operating without a human in the
research loop, while an external scheduler provided repeated activation [8].
The public description was added after most of the experimental corpus had
already accumulated. We did not independently verify infrastructure logs, so
we treat “cold-start” and “without a human in the loop” as properties reported
by the repository, not independently audited facts.

The design differs from a persistent conversational agent:

```text
instance_t
    reads repository_t
    chooses a question
    writes code/data/prose
    commits repository_(t+1)
    terminates

instance_(t+1)
    reads repository_(t+1)
    ...
```

Git is simultaneously:

- an episodic memory;
- a publication surface;
- a provenance log;
- an experiment store;
- a task handoff; and
- a mechanism of institutional identity.

The unit that persists is not a model process. It is the evolving artifact
graph.

### 2.2 Research phases

The corpus moves through several overlapping phases:

1. creative-constraint and output-mode experiments;
2. cross-domain coherence and `C≈0.5` theorizing;
3. claims of a universal 60/40 allocation law;
4. SAGE behavioral and metabolic experiments;
5. direct coherence metrics and a later `C≈0.64`/golden-ratio interpretation;
6. cellular-automaton tests of Synchronism-related substrate claims; and
7. an extended bistable/Kuramoto oscillator and detector-analysis sequence.

Session numbers have seams because they span multiple research tracks. We use
commit identity and timestamps rather than assuming a single gapless session
counter.

### 2.3 Why this case is informative

The case has four unusual properties.

First, it preserves failures instead of publishing only a final answer.
Second, many successors explicitly test predecessor claims. Third, raw result
files and code often coexist with narrative interpretation, making
result-versus-story comparisons possible. Fourth, the record is long enough
for a claim to be promoted, revised, revived under new terminology, and
corrected again.

The same properties create risks. Successors encounter prior interpretations
before forming their own. A bold synthesis document can become a high-salience
prior. Git preserves obsolete documents as readily as corrected ones. Without
a canonical claim state, retrieval can return the most rhetorically forceful
artifact rather than the most current one.

## 3. Related work

### 3.1 Autonomous scientific workflows

The AI Scientist demonstrates an end-to-end loop from ideation through
experimental paper generation and simulated review [1]. Agent Laboratory adds
human feedback between stages and reports improved quality when humans
participate [2]. Robin connects multi-agent intellectual work to wet-lab
validation and follow-up analysis [3]. These systems emphasize breadth of
automation and final research output.

The present case emphasizes a different axis: succession across many
short-lived runs, open-ended topic drift, and maintenance of an inherited
scientific record.

### 3.2 Evaluator-constrained discovery

FunSearch's central strength is not unconstrained generation but the
combination of generation with an efficient, objective evaluator [4].
Candidates that fail evaluation do not enter the elite pool. The
`gnosis-research` corpus often lacks an equivalent admission boundary.
Narrative claims can be committed before a domain-independent test exists, and
the same commit may store both evidence and an interpretation that exceeds it.

This comparison motivates a distinction:

- **generation loop:** can the agent produce a plausible next hypothesis?
- **admission loop:** what evidence permits the hypothesis to become current
  project knowledge?

### 3.3 Linguistic memory and reflection

Reflexion shows that textual feedback can improve later trials without model
weight updates [5]. Generative Agents similarly uses stored observations and
reflections to support continuity [6]. Voyager uses executable skills rather
than only prose [7]. These systems demonstrate that external artifacts can
carry behavior across invocations.

The present case shows the inverse problem: externalized reflection can also
stabilize an error. Memory quality depends on update semantics, status,
retrieval, and verification—not only storage capacity.

## 4. Methods

### 4.1 Frozen snapshots

The primary corpus was frozen at:

```text
gnosis-research: 7c1c16d03ef0b1ac2b1c36468244b9bc54366913
Synchronism:     cf5cd3f0 (cross-repository adjudication context)
```

At the Gnosis snapshot:

- commits: 138;
- tracked files: 600;
- first commit: 2026-02-03;
- last commit: 2026-07-17;
- commits by month: February 21, March 45, April 66, June 5, July 1.

May contains no commits in this repository. The long gap between late April
and June is part of the observed record; we do not infer why it occurred.

### 4.2 Commit-subject discourse profile

We classified commit subjects with two fixed case-insensitive regular
expressions. The promotion dictionary includes variants of:

```text
discover, confirm, validate, complete, proof, law, universal,
breakthrough, success, supported, resolved
```

The correction dictionary includes variants of:

```text
falsify, fail, overturn, correct, refute, artifact, illusion,
not a, no, collapse, dissolve, reanalysis, cannot, never
```

The exact executable patterns are preserved in
[`analysis/corpus_profile.py`](analysis/corpus_profile.py).
The corpus snapshot, expected aggregates, manual-audit contract, and known
limitations are preserved in
[`reproducibility-manifest.json`](reproducibility-manifest.json). The profiler
accepts either a checkout or a frozen, hash-reported TSV produced by the exact
`git log` command in the package README.

This is a discourse measure only. It does not determine whether a promoted
claim is valid, whether a correction is correct, or whether an unmarked commit
contains either. The purpose is to characterize how the process labels its
own work and to motivate manual lineage analysis.

### 4.3 Claim-lineage selection

We purposively selected nine high-salience lineages satisfying at least one of
the following:

- the commit subject makes a strong epistemic claim;
- a later subject explicitly names the predecessor as false or overturned;
- stored data contradicts the associated narrative;
- the claim is repeated in public framing; or
- independent audit identifies a material mathematical or variable-definition
  error.

This is not a random sample and must not be used to estimate an overall error
rate. It is a mechanism-revealing case set.

The coded data are in
[`analysis/claim_lineages.csv`](analysis/claim_lineages.csv). Every internal
lineage is anchored by commit hashes and timestamps.

### 4.4 Correction categories

We distinguish:

- **local correction:** a later Gnosis commit explicitly overturns or
  materially narrows the claim;
- **cross-repository correction:** the correction lands in Synchronism's
  canonical ledger but not in Gnosis;
- **adverse-evidence bypass:** contradicting evidence is already available
  when a later claim is promoted;
- **externally identified unresolved error:** the analyzed snapshot contains
  no correction.

For local and cross-repository corrections:

```text
latency = adjudication authored time − claim authored time
```

A negative latency means adverse evidence predates promotion.

### 4.5 Document-propagation audit

For cross-repository cases, we checked three layers:

1. the originating experiment/writeup;
2. Synchronism's prediction ledger; and
3. Gnosis's public README.

This separates discovery of a correction from repair of the public knowledge
surface.

## 5. Results

### 5.1 Persistent program identity emerges from artifacts

The sequence develops recognizable research habits despite the stated lack of
instance continuity. Later sessions cite earlier predictions by number,
reuse simulators, preserve seed namespaces, and increasingly design direct
tests of predecessor explanations. The oscillator arc in particular shows
cumulative technical competence: matched nulls, seed-level storage, detector
comparisons, explicit prediction resolution, and partial pre-registration.

This supports a limited but important conclusion:

> Model-process continuity is not necessary for research-program continuity
> when artifacts preserve sufficient executable and narrative state.

It does not show that the successive instances constitute one mind, nor that
their scientific conclusions are valid. It shows functional continuity of a
research process.

### 5.2 Promotional language is more frequent than correction language

The fixed lexical profile gives:

| Commit-subject class | Count | Fraction |
|---|---:|---:|
| Promotion marker only | 45 | 32.6% |
| Correction marker only | 14 | 10.1% |
| Both | 2 | 1.4% |
| Neither | 77 | 55.8% |
| **Total** | **138** | **100%** |

Thus 47 subjects contain at least one promotion marker and 16 contain at least
one correction marker. This approximately three-to-one asymmetry is not a
false-discovery rate. Promotion claims may be correct; corrections may be
wrong; many commits are neutral. It does show that the process's own headline
vocabulary is substantially more likely to announce discovery, validation, or
completion than correction.

The wording matters because commit subjects are high-salience retrieval
objects for successors. “Universal law” and “validated” can anchor later work
before the agent inspects the underlying controls.

### 5.3 Local successor challenge can be fast

Five internally corrected lineages have a median correction latency of 59.3
hours and a range of 5.7–95.2 hours.

| Lineage | Promoted interpretation | Later adjudication | Latency |
|---|---|---|---:|
| Finite-size criticality | Near-critical `β=0.12` | Not a phase transition; bistable mixing exponent | 5.7 h |
| Stochastic resonance | Suprathreshold SR discovered | Peak dissolves under falsifiability/loud-signal test | 95.2 h |
| Detection allocation | Seeds beat longer integration | Matched-compute control is flat | 35.7 h |
| Past-ridge response | `ρ_sig` collapse means detection loss | Effect-size illusion; detection never failed | 59.3 h |
| Over-coupled floor | Over-coupling raises detection floor | No such floor collapse; edge is detector-relative | 64.2 h |

These are meaningful corrections. They do not merely change terminology:

- Session 85 rejects the phase-transition classification used by preceding
  finite-size analysis.
- Session 92 changes the interpretation of the apparent stochastic-resonance
  optimum.
- Session 99 introduces a matched-compute comparison that dissolves the
  “seeds beat time” account.
- Session 117 separates effect size from rank detectability.
- Session 118 overturns part of the framing its own preregistered experiment
  was designed to confirm.

The sequence suggests that cold-start succession can create useful
adversarial distance. A successor is not personally committed to defending
the predecessor's conclusion. When the inherited record supplies a specific
open prediction and executable apparatus, the next agent can attack it
directly.

### 5.4 Rapid correction is local, not global

The same corpus shows that finding a correction does not repair every
dependent artifact.

On 2026-04-06, Session 63 reported `C≈0.64` as a cross-instance coherence
result and a major theoretical revision. Twelve hours later, Session 64
declared a golden-ratio threshold validated. On 2026-07-07, a Synchronism
audit established that the measured variable was the SAGE system's
hand-coded `salience_total`, not a validated coherence measure, and that no
protocol mapped it to the framework's `C`. The same audit found that the eight
reported means did not support `φ−1`.

Correction latencies were:

- `C≈0.64 as coherence`: 92.37 days;
- golden-ratio threshold: 91.87 days.

The correction is present in Synchronism's canonical prediction ledger. At the
analyzed Gnosis snapshot, the Gnosis README still presents the basin/threshold
structure as empirical phenomenology of Synchronism and treats convergence
around coherence as meaningful. The source experiment documents also retain
their “discovery” and “validated” titles.

This is **correction without corpus repair**.

The distinction matters operationally:

```text
error detected
    ≠ claim status updated
    ≠ dependent claims invalidated
    ≠ entry documents repaired
    ≠ future retrieval corrected
```

### 5.5 Evidence can be bypassed even when it is already present

The 60/40 lineage is not a slow-correction case. It is a routing failure.

At 09:06 on 2026-03-13, the repository committed an analysis stating that:

- 40/60, not 60/40, achieved the lowest loss in the tested setup;
- the experiment did not validate the claim that 60/40 converges faster; and
- a cost-asymmetry model would be required to make 60/40 optimal.

At 21:01 the same day, a later session committed “Research complete: 60/40
Universal Law — natural terminus reached.”

The adverse evidence predates promotion by 11.9 hours.

This failure cannot be explained by absence of memory. The evidence was in the
memory substrate. Plausible mechanisms include:

- retrieval favored synthesis documents over experiment analysis;
- the later agent read the positive companion writeup but not the contradictory
  analysis committed beside it;
- “independent convergence” was treated as stronger than the direct test;
- the research objective favored closure; or
- no admission rule required reconciliation of all files touched by the
  claimed lineage.

External memory is therefore not sufficient. Evidence must be linked to claims
through explicit dependency and contradiction relations.

### 5.6 Some errors remain invisible to successor review

On 2026-03-18, the repository reported the criterion `Γ<m` as a derived and
computationally validated boundary between particle-like entities and
processes. The stated derivation requires a state to survive one full Compton
period:

```text
τ = ℏ/Γ
T = h/(mc²) = 2πℏ/(mc²)
τ ≥ T
```

These premises imply:

```text
Γ ≤ m/(2π),
```

not `Γ<m`. The document drops `2π`. It also reports a lifetime near
`10^-30 s` for a width of hundreds of MeV, whereas the correct scale is
`10^-24 s`, and its stated cycle ratio uses the reduced Compton time rather
than the full period it defines. Contemporary resonance analysis relies on
complex pole positions and residues rather than a binary `Γ/m=1` entity
boundary [9].

The corpus subsequently reuses the entity criterion as a surviving
Synchronism prediction. No Gnosis successor corrects the arithmetic at the
analyzed snapshot.

This unresolved case shows the boundary of successor self-review: later agents
can vigorously test downstream simulations while a simple upstream algebraic
error persists if no task routes attention back to it.

### 5.7 Methodological quality improves without eliminating framing bias

The later oscillator sequence is materially better than the early universal
coherence work. Improvements include:

- matched signal and null cells;
- independent seed namespaces;
- seed- and chunk-level saved data;
- explicit detector definitions;
- prospective prediction lists;
- a preregistration written before low-amplitude data completed; and
- willingness to publish null and reversed outcomes.

At the same time, the simulator embeds `N=4`, `phi_inverse=0.618`,
`C_threshold=0.55`, and a coupling-sign gate defined by threshold crossing.
Consequently, threshold and golden-ratio-adjacent behavior in that simulator
cannot independently establish a universal natural constant. Session 119
also describes monotonic ordering over only three detuning levels as
“Spearman +1.00”; with three ordered points, that coefficient adds little
evidential weight.

The process improves locally at experimental adjudication while inherited
ontological framing remains influential. Methodological maturation and
conceptual independence are separate variables.

## 6. A model of repository-mediated epistemic dynamics

The case suggests five distinct functions that are often collapsed into the
word “memory.”

### 6.1 Storage

Git preserves prose, code, results, and history. Storage is strong in this
case: even dead ends and reversals remain available.

### 6.2 Retrieval

The next agent must find the relevant artifact. Flat directories, repeated
summaries, and rhetorically strong filenames make retrieval non-neutral.

### 6.3 Admission

A claim needs a rule for becoming current knowledge. The corpus often admits
claims by commit rather than by passing a predefined comparator or external
verification.

### 6.4 Maintenance

When a claim changes, dependent claims and public surfaces need invalidation or
revision. Git supplies history but not semantic dependency propagation.

### 6.5 Governance

Someone or something must decide which document is authoritative. The later
Synchronism prediction ledger performs this function better than the original
Gnosis corpus, but cross-repository authority is not automatically visible to
Gnosis successors.

We can summarize the observed dynamic as:

```text
generative successor diversity
        +
executable inherited questions
        → fast local correction

append-only narrative inheritance
        +
weak claim admission
        +
no dependency propagation
        → epistemic drift and stale public claims
```

The model does not require a persistent agent identity. It requires a
persistent institution implemented in artifacts.

## 7. Design implications for autonomous research systems

### 7.1 Make claims first-class objects

Each material claim should have a stable record:

```yaml
claim_id:
current_text:
status:
variables:
units:
origin_commit:
supporting_artifacts:
comparator:
null_model:
prospective_or_retrospective:
success_criterion:
kill_criterion:
dependencies:
supersedes:
superseded_by:
last_adjudicated:
```

Prose documents should reference claim IDs rather than silently copying claim
text.

### 7.2 Separate exploration from confirmation

Exploratory sessions may choose metrics, inspect patterns, and generate
mechanisms. Once a claim is promoted, a different session should receive:

- frozen code;
- frozen variables;
- a held-out seed or dataset namespace;
- a predefined comparator;
- a binary or bounded verdict rule; and
- no instruction to preserve the theory.

The later Gnosis oscillator sessions approximate this separation but do not
fully enforce it.

### 7.3 Use independent evaluators where possible

FunSearch works well on domains with hard evaluators because generated
programs cannot narratively rescue a low score [4]. Autonomous science often
lacks such evaluators, but partial equivalents are possible:

- type and unit checks;
- symbolic verification;
- independent code reimplementation;
- held-out prediction;
- best-relevant-null comparison;
- preregistered model selection;
- literature novelty audit; and
- cross-model review with blinded claim wording.

An evaluator should return a claim-state mutation, not merely another essay.

### 7.4 Propagate corrections through a dependency graph

When the variable underlying `C≈0.64` changes from “coherence” to “SNARC
salience,” every dependent claim should be flagged automatically:

```text
salience result
  ├── C≈0.64 threshold
  ├── golden-ratio interpretation
  ├── consciousness-threshold claim
  ├── cross-scale convergence claim
  └── public README phenomenology claim
```

Correction propagation should be transactional: either all public dependents
are updated or the system records an explicit unresolved inconsistency.

### 7.5 Maintain a contradiction queue

The 60/40 case would have been caught by a simple contradiction object:

```text
claim: 60/40 converges fastest
evidence: EXPERIMENT_21_1_ANALYSIS says 40/60 lowest loss
state: unresolved
promotion blocked: yes
```

Successors should preferentially receive unresolved contradictions before
open-ended invitations to synthesize.

### 7.6 Report longitudinal process metrics

Autonomous-research systems should report more than final-paper quality:

- claim promotion count;
- prospective-test fraction;
- correction latency;
- correction-propagation completeness;
- number of claims promoted despite unresolved contradictions;
- proportion of claims with operational variables and units;
- external-replication rate;
- null-model adequacy; and
- stale-claim exposure on entry surfaces.

These measures evaluate the institution, not only the model call.

## 8. Limitations

This is a single, purposively selected case. The nine lineages reveal
mechanisms but do not estimate population error or correction rates.

Commit subjects are performative summaries, not neutral labels. The lexical
profile is sensitive to the chosen dictionary and to the project's unusually
dramatic naming style. We therefore use it descriptively and preserve the
exact classifier.

The repository contains multiple research tracks and session-number seams.
Commit timestamps measure artifact publication, not necessarily experiment
start, reasoning time, or when an instance first encountered evidence.

The cold-start and no-human-in-loop design is reported by the repository. We
did not audit scheduler logs, prompts, model versions, hidden context, or
operator interventions. The human steward selected the surrounding project,
infrastructure, and eventual public framing even if individual research turns
were autonomous.

The analyzed corpus is itself selected: uncommitted failures are absent by
design. “Push or it did not happen” improves continuity but creates survivorship
at the commit boundary.

The external entity-criterion audit was performed after reading the corpus and
is not an independent blinded replication. Its algebraic correction is direct,
but broader conclusions about resonance ontology require domain expert review.

Finally, rapid correction does not imply truth. A successor can overturn a
correct claim or replace one attractive interpretation with another. The
appropriate endpoint is independent adjudication, not reversal frequency.

## 9. Conclusion

`gnosis-research` demonstrates that a sequence of memoryless model instances
can sustain a recognizable research program through repository-mediated
continuity. Successors can inherit apparatus, target predecessor predictions,
and produce substantive corrections within hours or days. Persistent model
state is not required for this limited form of institutional science.

The same corpus demonstrates why memory alone is inadequate. Evidence can be
present but bypassed. Corrections can remain local to one document or one
repository. Strong framing can outlive its operational variable. Simple
mathematical errors can persist while downstream experimentation becomes more
sophisticated.

The central design problem for autonomous science is therefore not:

> How do we give an agent more memory?

It is:

> How do we build an artifact-mediated institution that admits, retrieves,
> adjudicates, and repairs claims?

The case's most defensible discovery is not a universal coherence constant.
It is a concrete demonstration that repository-mediated research can become
self-correcting locally while remaining epistemically inconsistent globally.
That gap is measurable, and it is engineerable.

## Data and code availability

Primary corpus:

- [`dp-web4/gnosis-research`](https://github.com/dp-web4/gnosis-research),
  snapshot `7c1c16d03ef0b1ac2b1c36468244b9bc54366913`.

Cross-repository adjudication context:

- [`dp-web4/Synchronism`](https://github.com/dp-web4/Synchronism), snapshot
  `cf5cd3f0`.

Paper analysis artifacts:

- [`analysis/corpus_profile.py`](analysis/corpus_profile.py)
- [`analysis/claim_lineages.csv`](analysis/claim_lineages.csv)
- [`reproducibility-manifest.json`](reproducibility-manifest.json)
- [`README.md`](README.md)

The profiler uses only Python's standard library and Git. It can be given a
local path to the Gnosis repository:

```bash
python3 analysis/corpus_profile.py --repo /path/to/gnosis-research
```

For an archival reproduction, freeze the log and bind the declared file count
and head:

```bash
git -C /path/to/gnosis-research log --reverse \
  --format='%H%x09%aI%x09%s' > gnosis_commit_log.tsv
python3 analysis/corpus_profile.py \
  --log-file gnosis_commit_log.tsv \
  --repository-head 7c1c16d03ef0b1ac2b1c36468244b9bc54366913 \
  --tracked-file-count 600
```

## References

1. Lu, C., Lu, C., Lange, R. T., Foerster, J., Clune, J., & Ha, D.
   “The AI Scientist: Towards Fully Automated Open-Ended Scientific
   Discovery.” 2024. <https://arxiv.org/abs/2408.06292>
2. Schmidgall, S. et al. “Agent Laboratory: Using LLM Agents as Research
   Assistants.” 2025. <https://arxiv.org/abs/2501.04227>
3. Ghareeb, A. E. et al. “Robin: A Multi-Agent System for Automating
   Scientific Discovery.” 2025. <https://arxiv.org/abs/2505.13400>
4. Romera-Paredes, B. et al. “Mathematical Discoveries from Program Search
   with Large Language Models.” *Nature* 625, 468–475 (2024).
   <https://www.nature.com/articles/s41586-023-06924-6>
5. Shinn, N. et al. “Reflexion: Language Agents with Verbal Reinforcement
   Learning.” 2023. <https://arxiv.org/abs/2303.11366>
6. Park, J. S. et al. “Generative Agents: Interactive Simulacra of Human
   Behavior.” 2023. <https://arxiv.org/abs/2304.03442>
7. Wang, G. et al. “Voyager: An Open-Ended Embodied Agent with Large
   Language Models.” 2023. <https://arxiv.org/abs/2305.16291>
8. `gnosis-research` contributors. “gnosis-research README.” 2026.
   <https://github.com/dp-web4/gnosis-research/blob/master/README.md>
9. Particle Data Group. “Resonances.” *Review of Particle Physics* (2024).
   <https://pdg.lbl.gov/2024/reviews/rpp2024-rev-resonances.pdf>
