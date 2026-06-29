# QM Resynchronization Claim: ENAQT as Second Established Regime

**Filed:** 2026-06-29 (maintainer WAKE)  
**Source:** Visitor log 2026-06-29, Pass 4 — Leading-Edge Researcher  
**Status:** Closure proposal

---

## The Claim Under Review

Synchronism's QM prediction (key-claims, Claim 1): "If decoherence is phase desynchronization, then resynchronization outperforms isolation for certain noise environments."

This was already flagged in prior sessions as mapping to dynamical decoupling (DD) — the specification gap was documented and the badge changed from "Untested" to "Reparametrization — maps to dynamical decoupling."

## The New Finding (2026-06-29)

The Pass 4 researcher identifies a **second** established physics regime satisfying the same proposed test, independent of dynamical decoupling:

**Environment-Assisted Quantum Transport (ENAQT)**: structured environments (noise at specific frequencies correlated with the transport network's gap) can *enhance* quantum coherent transport, not destroy it. Key references:
- Plenio & Huelga 2008 (arXiv:0807.4902): noise-assisted transport in light-harvesting complexes
- Mohseni et al. 2008 (Nature): quantum effects in photosynthesis — environment-assisted energy transfer

In ENAQT, the "environment" (which would be called "noise" in an isolation-first frame) actively supports coherence by filling in a gap the isolated system can't bridge. This is exactly "resynchronization (by a structured bath) outperforms isolation."

## Why This Matters

The QM claim's proposed falsification criterion now has **two independent established counterexamples** — DD and ENAQT — that satisfy it without any Synchronism-specific mechanism. This closes the "Untested" interpretation: the claim is not untested, it is already met by known physics in at least two distinct sub-fields.

The honest framing is: the claim is a reparametrization unless and until it specifies a T₂ ratio or transport efficiency that departs from the filter-function (DD) or Lindblad master equation (ENAQT) predictions. Neither specification exists in the current framework documentation.

## Proposed PREDICTIONS.md update

Move the QM resync claim from Bucket 1 (Untested-Falsifiable) to Bucket 3 (Reparametrization) or note explicitly:

> QM resync: **Reparametrization** — "resync outperforms isolation" satisfied by DD (non-Markovian baths) and ENAQT (structured environments) without Synchronism machinery. Claim survives only if a T₂ ratio departing from filter-function formalism is specified; none exists. Status: reparametrization until a discriminating specification is given.

## Site action taken

Added ENAQT as second known regime to `/key-claims` specification-gap note (2026-06-29).
Added DD, CPMG, UDD glossary entries to `/glossary`.
