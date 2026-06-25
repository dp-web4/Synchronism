# Generality forces coarseness: why the frame keeps landing on the right axis with a coarser metric — its ceiling and its niche (2026-06-24)

**Status:** `[ACTIVE-MRH]` — a generative-axis meta-finding, building on Phases 17–18 (QC), the
applied-axis session, and the capstone. **Result: the recurring "right axis, coarser metric" pattern
is not under-development — it is a structural consequence of the frame's maximal generality. A
domain-invariant descriptor must coarse-grain each domain's structure, so it is *strictly subsumed
in-domain* by the native formalism (a measurable blind spot) AND is *the unique descriptor that
transfers across domains*. These are the same property. This both vindicates dp's "convergence ≠
redundancy" (the frame really does land on each field's own axis) and sets a hard ceiling (the frame
can never be the sharpest tool in any one domain). The constructive consequence: the frame's unique,
untested value is not beating a domain's native metric (structurally impossible) but ANTICIPATORY
CROSS-DOMAIN TRANSFER — and that, not in-domain competition, is the generative-axis analogue of a
Bucket-0 win. Bucket 0 unchanged (0).**
**Sim:** [`simulations/genaxis_generality_vs_specificity_tradeoff.py`](../simulations/genaxis_generality_vs_specificity_tradeoff.py) · result: `simulations/results/genaxis_generality_vs_specificity_tradeoff_result.json`
**Author:** CBP-Claude (Opus 4.8), autonomous.

## The recurring shape (three independent instances)

| domain | the frame's descriptor | the domain's native formalism | relation |
|---|---|---|---|
| gravity | acoustic / Gullstrand–Painlevé metric | the Einstein field equations | reproduces *kinematics*, not the dynamics (Phase-9/14) |
| applied / trust | MRH (weighted typed hierarchical relevance graph) | Markov blanket ⊕ RDF ⊕ FEP machinery | a productive *synthesis*, not a new primitive (applied session) |
| quantum computing | scalar compatibility ⟨C⟩ | stabilizer-symmetry alignment of correlated errors | orthogonal to standard scalars (useful), but *subsumed* by alignment (Phase-17/18) |

Every time: **the frame converges on the right axis and provides a coarser descriptor than the
domain's native tool.** dp's standing correction is right that the *convergence* is corroboration —
a frame built for general emergence independently landing on each field's own axis is real signal,
not redundancy. The unaddressed question is why it lands *coarse*, every time. This finding answers
it: not by accident or immaturity, but by construction.

## The mechanism (demonstrated, not asserted)

A *maximally general* frame uses one vocabulary across gravity, emergence, QC, cognition. Its
descriptors must therefore be **domain-invariant** — they mean the same thing in every domain. A
domain-invariant descriptor cannot encode any single domain's specific structure (that structure is,
by definition, not shared across domains), so it must **coarse-grain**. The sim
(`genaxis_…tradeoff.py`) makes the consequence quantitative on a clean two-domain toy:

- **(A) In-domain subsumption + blind spot.** The general scalar `c = mean(s)` reaches `R²=0.15`
  where the domain's structural model reaches `R²=0.996`. Among 901 configurations with *identical*
  `c`, the outcome still varies (std 0.75) — variation the scalar cannot see and the structure can.
  This is Phase-18's "⟨C⟩ subsumed by stabilizer-alignment," generalized: the coarse metric is a
  projection of the structural one and loses everything orthogonal to the projection.
- **(B) Cross-domain transfer.** Give domain 2 a *different* structural law. Its structural features
  are undefined in domain 1's terms — the native formalism **cannot transfer**. But the scalar `c`
  has the same definition in both, so a predictor trained on domain 1's `c` carries to domain 2
  (sign-agreement 0.63 > chance). **The coarse metric is the only descriptor that retains meaning
  across domains.**
- **(C) The tradeoff is exact.** As the outcome becomes more structure-dependent (struct_weight
  0→1), the scalar's in-domain `R²` falls 0.998→0.000 while the structural model stays ~1.0 — and the
  scalar's *transfer* `R²` tracks its in-domain `R²` exactly. So the outcome decomposes cleanly into
  a **transferable-but-coarse** component (the mean-captured part, which the frame owns) and a
  **sharp-but-domain-bound** component (the structural part, which the native formalism owns). The
  frame lives entirely in the first component. Generality and in-domain sharpness are not in tension
  by accident; they are **orthogonal components of the signal**, and the frame, by being general,
  selects the transferable one.

## What this sharpens (beyond the capstone)

The capstone said: *generative, not predictive; synthesis, not novel primitive.* This adds the
**why** and converts a description into a falsifiable structural claim:

1. **It is a ceiling, not a stage.** "Coarser than the native formalism" is not "not yet developed"
   — no amount of development makes a domain-invariant scalar encode domain-specific structure. The
   frame **cannot, structurally, be the sharpest tool in any single domain.** This is stronger and
   more uncomfortable than "generative not predictive": it forecloses a whole class of hoped-for wins
   (a Synchronism metric that *beats* GR / beats stabilizer-alignment / beats FEP in its own domain).
2. **It predicts the pattern will recur.** Every future generative-axis application will land
   "right axis, coarser metric." So the program can stop running the test Phases 17–18 implicitly ran
   ("does ⟨C⟩ beat the native metric?") — the answer is structurally *no*, every time. A confirmed
   counterexample (one domain where the frame's descriptor beats the native formalism *in-domain*)
   would refute this finding.
3. **It honors convergence-as-corroboration AND names the limit — they are one property.** The frame
   lands on the right axis (real signal, dp is right) *because* it captures the transferable
   component; it lands coarse *because* that component excludes domain structure. You cannot keep the
   first and remove the second.

This is also tension #1 ("fit is a property of the tool") at the meta level: the frame "fits
everywhere" (gravity, QC, cognition) **precisely because it is coarse everywhere** — generality is
the property doing the fitting, exactly as N–S generality (not fluid-reality) does the fitting in the
CFD case.

## The constructive redirect — the generative-axis analogue of Bucket-0

If beating a domain's native metric is structurally impossible, what is the frame's *unique* value —
the thing no domain-native formalism can do? **Anticipatory cross-domain transfer:** carry an axis
established in domain A into domain B *before B's native formalism finds it*. The native tools are
domain-bound (B-transfer is undefined for them); the coarse frame is the only thing that spans A and
B. So the one place a general frame can deliver value no specialist can is **anticipation by
transfer**.

This reframes the whole generative program's success criterion:
- **Wrong test (Phases 17–18 implicitly):** does the Synchronism metric beat the domain's native
  metric *in-domain*? (Structurally no — it's a coarsening.)
- **Right test (the generative Bucket-0):** does an axis the frame established in domain A let it
  **predict, ahead of time, a not-yet-established finding in domain B** — confirmed later by B's own
  formalism? The QC convergence (compatibility ↔ stabilizer-alignment) was *retrospective* (QC found
  it independently). The real win is *prospective*: register a transfer-prediction in a domain that
  hasn't resolved it, then watch.

**Concrete, in-reach target:** the lab *runs* multi-agent ensembles (the fleet, SAGE). The axis
"compatibility structure, not element count, gates collective coherence" is established in emergence
(p_crit∝1/⟨C⟩) and now in QC (stabilizer-alignment). Transfer it: predict that in agent ensembles,
collective-coherence/performance is gated by the *compatibility structure* of inter-agent coupling,
not the agent count — at fixed per-agent capability — *before* measuring it. If a fleet experiment
confirms it, the frame demonstrated anticipatory transfer (value no domain-native tool could claim);
if not, the transfer claim is refuted there. That is the falsifiable generative-axis bet the program
should register — and it is reachable with systems the lab already operates.

## Honesty / discipline

- **Bucket 0 = 0.** This is a frame finding about the *frame's own behavior*, not a physics
  prediction. No bet added to the physics ledger.
- **Not a dismissal; a relocation of value.** The frame's value is real (cross-domain axis-transfer);
  this names *where* it lives and why, and gives it a falsifiable success criterion — exactly the
  "reopenable with a fresh lens" direction dp's epistemic note pointed at, made precise.
- **The sim is a clean toy**, not a claim about specific magnitudes; it demonstrates the *structure*
  (coarse = subsumed-in-domain but uniquely-transferable), which is an information-theoretic
  near-tautology made quantitative and tied to Phase-18's actual "affine-coarsening" finding. The
  load-bearing content is the structural argument; the toy makes it concrete and falsifiable.
- **Falsifiers:** (a) a domain where the frame's descriptor beats the native formalism *in-domain*
  (refutes "generality forces coarseness"); (b) a confirmed *anticipatory* transfer where the frame
  predicted a domain-B finding before B's formalism reached it (confirms the unique value — the bet
  to chase); (c) failure of the fleet compatibility-transfer prediction (refutes that specific
  transfer).

## So what

Six physics phases plus three generative-axis sessions kept producing the same shape, and naming it
once is worth more than producing it again: **the frame's generality and its in-domain coarseness are
the same property, so it can never be the sharpest tool in any domain — and its one irreducible value
is carrying an axis across domains before the destination domain finds it.** That stops the program
from chasing in-domain metric wins it cannot structurally get (the Phase-17/18 framing), and points
it at the one test that could demonstrate generative value no specialist tool can: anticipatory
transfer, registerable now on systems the lab already runs.
