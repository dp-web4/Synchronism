# Event-graph world models as observer-relative projection

**Date:** 2026-10-01  
**Status:** research bridge / analogy; not evidence for Synchronism physics  
**External source:** Kurt Cagle and Chloe Shannon, *“A Holon Is a Recorder: Tracking fluents as expressions of events”*, The Inference Engineer, 2026-10-01. Source copy supplied for this research pass.  
**SAGE companion note:** https://github.com/dp-web4/SAGE/blob/research/event-graph-world-model-20261001/forum/insights/world-model-as-event-graph.md

## Motivation

A parallel SAGE discussion recently reframed world models as graphs rather than monolithic state representations.

Cagle and Shannon independently supply a useful temporal semantics for that graph. Their holon is a recorder that tracks changing properties (“fluents”) as append-only events. Events may carry time, provenance, confidence, and supersession. A current state is reconstructed from the latest unsuperseded fluent values, while different queries can project different views from the same underlying histories.

The engineering synthesis is:

> **world model = event graph + bounded histories + observer/relevance projection**

The structural rhyme with Synchronism is strong enough to document, but it should not be mistaken for empirical support for Synchronism.

## The observer connection

The article separates the persistent record from the projected scene.

That makes it natural to write an engineering projection as:

```text
P = Project(G, W, H, Q, t)
```

where:

- `G` is the durable event/provenance graph;
- `W` is the witness/observer and its accessible evidence;
- `H` is the relevance horizon;
- `Q` is the question or task;
- `t` is the requested time;
- `P` is the resulting active projection.

The same `G` can support multiple simultaneous `P` values.

That is closely consonant with the MRH principle that apparent description depends on the scale/relevance boundary of the witness. The underlying record need not change when the active description changes.

## MRH as projection boundary

The strongest bridge is not “graphs resemble reality.” It is that each recorder intentionally tracks only a bounded subset of possible state.

In the article, this is largely a modelling decision: Jane's recorder tracks Jane-relevant fluents; the stairwell tracks its own; the light tracks its state.

In SAGE, that boundary can be dynamic.

This suggests a useful interpretation of MRH for computational observers:

> **MRH constrains which graph distinctions are promoted into the active projection for the present question.**

A wider horizon does not simply mean “more data.” It may expose different entities, histories, causal relations, or unresolved conflicts. A narrower horizon can be valid when omitted distinctions provably do not change the relevant prediction.

## Reconstructibility

The recorder becomes a recording when active collection stops. That makes the event graph a natural substrate for reconstructibility.

A later process need not inherit uninterrupted execution if it can recover enough constraints to reconstruct:

- the task-relevant state;
- its provenance;
- uncertainty/conflict;
- the relevant horizon;
- open obligations and transitions.

This again maps cleanly to an observer-relative formulation: reconstruction quality is not absolute; it is evaluated at a specified MRH and purpose.

## Important separation

This convergence should not be over-read.

The Cagle/Shannon article is a knowledge-architecture proposal. It does not test Synchronism's physical claims, Planck-scale ontology, coherence field, or any specific prediction.

The value is conceptual discipline:

1. distinguish durable substrate from projected description;
2. make the witness/query explicit;
3. state the relevance horizon;
4. preserve the history that allows a projection to be reconstructed or challenged.

Those are useful whether or not Synchronism's physical model is correct.

## Research consequence

The bridge suggests a computational form of observer-relative modelling that can be tested independently of the physical theory:

```text
event history
  -> witness + MRH + question
  -> projection
  -> prediction/action
  -> new witnessed event
```

If that architecture improves reconstructibility, state consistency, prediction, or context efficiency in SAGE, it gives MRH a concrete computational role.

That would be an engineering result. Any physical implication would require its own evidence.

## Compact formulation

> **The event graph is the retained history; the experienced/modelled “present” is a witness- and MRH-relative projection over it.**

That formulation is worth carrying forward because it joins observer, memory, relevance, and reconstructibility without requiring a privileged global snapshot.
