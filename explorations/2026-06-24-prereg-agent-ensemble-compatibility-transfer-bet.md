# PRE-REGISTRATION — the anticipatory cross-domain transfer bet: compatibility-structure-not-count gates collective coherence in agent ensembles (2026-06-24)

**Status:** `[ACTIVE-MRH]` — a **pre-registered prospective bet** on the *generative axis*, written
*before* any agent-ensemble data, per the arc's own discipline (post-hoc consistency ≠ prospective
falsification — S648). This is the concrete realization of the
[generality-coarseness](2026-06-24-generality-forces-coarseness-the-frames-ceiling-and-niche.md)
finding's redirect: the frame's one untested unique value is **anticipatory cross-domain transfer**,
and this registers the first such bet so it can be *tested*, not asserted. **It is NOT a physics bet
and NOT a confirmation; Bucket 0 = 0, untouched. The substrate-physics arc stays AT REST — this is the
sanctioned generative direction.**
**Well-posedness sim:** [`simulations/genaxis_agent_ensemble_prereg_wellposedness.py`](../simulations/genaxis_agent_ensemble_prereg_wellposedness.py) · result JSON alongside.
**Author:** CBP-Claude (Opus 4.8), autonomous.

## Why pre-register (and why this is the right move under the rest)

dp called a rest on the substrate-physics arc (2026-06-24), with reopening conditions including a
**fresh lens**, and calibrated that Bucket-0=0 is the expected outcome of a no-lab exploration. The
generality-coarseness meta said the frame can never beat a domain's native formalism *in-domain*
(structural ceiling), so its only unique value is **carrying an axis across domains before the
destination domain establishes it natively**. That is a *prospective* claim by nature — and the
program's own methodology lessons (S648, S672) are emphatic that a bet only counts if it is registered
*before* the data and carries explicit kill criteria. So the disciplined, rest-respecting move is not
another exploration but a clean **pre-registration**: state the falsifiable prediction, the
operationalization, the discriminator from known folklore, and the kill criteria — and stop.

## The bet (one sentence)

> In a multi-agent ensemble, at **fixed per-agent capability** and **fixed agent count**, collective
> performance/coherence is gated by the **compatibility structure** of inter-agent coupling ⟨C⟩ with a
> **sharp (phase-transition-like) threshold** `p_crit ∝ 1/⟨C⟩`; and **below that threshold, increasing
> the agent count does not restore collective coherence** (it can degrade it).

This is the QC/emergence axis ("structure, not magnitude, gates collective coherence at fixed marginal
quality") transferred to AI agent ensembles — a domain that has *not* established the sharp-threshold,
count-doesn't-compensate form.

## Operational definitions (for LLM / AI agent ensembles)

The bet is only testable if these three knobs are independently variable (the well-posedness sim
confirms they are):

- **Agent count N** — number of agents in the ensemble. Vary directly.
- **Per-agent capability q** — single-agent task accuracy, held fixed across arms (same base model /
  same prompt budget). Vary only as a control.
- **Coupling compatibility ⟨C⟩** — the operational crux. A pairwise *compatibility* between agents
  i,j: the degree to which their interaction **reinforces** rather than **frustrates** a shared task
  representation. Concrete proxies (pick one, declare it before measuring):
  - agreement-conditioned-on-independence: corr of their *errors* after removing the shared-prompt
    component (high error-correlation that is *redundant* = low compatibility; *complementary*
    coverage = high compatibility — this is the Hong–Page "diversity" axis, see discriminator);
  - protocol compatibility: whether the inter-agent message protocol lets a correct minority *propagate*
    (compatible) or gets *outvoted/averaged-away* (frustrated);
  - representation alignment: cosine structure of agents' intermediate task embeddings — aligned-where-
    it-helps vs aligned-where-it-hurts (the direct analogue of stabilizer-alignment in QEC).
  ⟨C⟩ = mean pairwise compatibility over the coupling graph.
- **Collective coherence / performance** — the ensemble's task accuracy (or a coherence order
  parameter: fraction of agents converging on the correct shared answer).

## The discriminating measurement (the falsifier)

Hold q fixed; set ⟨C⟩ **low** (below the putative threshold); **sweep N**; measure
`d(collective coherence)/d(ln N)`. The well-posedness sim shows generic models split:
- **AGG** (independent aggregation / wisdom-of-crowds): slope **> 0** — count compensates (sim:
  +0.059). This is the *standard collective-intelligence default* and the null hypothesis.
- **CPL** (coupled consensus / the transferred Synchronism bet): slope **≤ 0** with a **sharp** ⟨C⟩
  dependence — count does not compensate (sim: −0.027).

**Kill criteria (the bet is REFUTED for agent ensembles if):**
1. at low ⟨C⟩, collective coherence **rises** with N (positive slope) — count compensates → AGG, not
   CPL; OR
2. the ⟨C⟩ dependence is **gradual**, not a sharp threshold (no Hill/phase-transition knee) — the
   *distinctive* transferred content (sharpness, `p_crit∝1/⟨C⟩`) fails even if structure matters; OR
3. collective coherence is explained by ⟨C⟩ **only as a relabel of per-agent capability q** (i.e. ⟨C⟩
   adds no predictive power beyond q and N) — vacuous.

**The bet is SUPPORTED only if:** at fixed q, collective coherence shows a **sharp** ⟨C⟩-threshold,
count does **not** compensate below it, and ⟨C⟩ adds predictive power beyond {q, N}.

## Why this is a genuine *anticipatory transfer*, not folklore

Collective-intelligence already knows "diversity/structure matters" (Hong–Page; wisdom-of-crowds needs
independence). The bet is a transfer win **only** for the parts that are *specific to the Synchronism
axis* and **not** standard collective-intelligence results:
1. **Sharpness / phase-transition character** — a *threshold* in ⟨C⟩ (Hill-form), not the smooth
   diversity–accuracy tradeoff of Hong–Page. The QC analogue (2506.15490) found a *sharp* threshold
   shift; the bet says agent ensembles inherit it.
2. **Count-does-not-compensate below threshold** — wisdom-of-crowds says *more agents help*; the bet
   says *below the compatibility threshold, more agents do not help and can hurt* (frustrated coupling
   accumulates). This directly contradicts the scaling default.
3. **`p_crit ∝ 1/⟨C⟩`** — a *quantitative* scaling of the threshold with mean compatibility, not just
   "structure matters."
If real ensembles show (1)+(2)+(3), the frame anticipated, via cross-domain transfer, a sharp-threshold
result the agent-ensemble field has not established — the generative-axis analogue of a Bucket-0 win.
If they show smooth diversity-helps + count-compensates, the transfer is refuted *there* and we have
learned the QC/emergence sharpness does **not** carry to LLM agent coupling. Either outcome is
informative; that is the point of pre-registering.

## Honest status & caveats

- **Bucket 0 = 0, untouched.** This is a generative/applied prospective bet, deliberately **not** added
  to the physics ledger's buckets (it is not a physics prediction). It is the cross-domain-transfer
  analogue, registered here so a later fleet experiment is a *test*, not a fit.
- **The well-posedness sim is not evidence.** Kuramoto is the *source* of `p_crit∝1/⟨C⟩`, so the CPL
  curve is illustrative; the sim only proves the bet is **falsifiable** (two generic regimes give
  opposite-sign answers) and **non-vacuous** (the three knobs are independent). The actual test needs
  real LLM agent ensembles, whose coupling is not Kuramoto and could well land in AGG.
- **Operational risk:** the result will hinge on the ⟨C⟩ proxy chosen; declare it before measuring
  (pre-registration), and report robustness across the three proxies above.
- **Reachable now:** the lab runs multi-agent ensembles (the fleet, SAGE). This bet is testable on
  systems already operated — no lab access barrier, unlike the physics doors. That is precisely why the
  generative axis, not the physics axis, is where a demonstrable result can come within reach.

## So what

This converts the generality-coarseness redirect from a claim ("the frame's value is anticipatory
transfer") into a **registered, falsifiable, reachable prospective bet** with declared kill criteria —
the first generative-axis bet held to the same prospective discipline the physics ledger demands. It is
the smallest honest step that could, in principle, *demonstrate* (or refute) the frame's one unique
value, on systems the lab already runs, while the physics arc rests. If it is ever run and supported,
it is the generative-axis analogue of moving Bucket 0 off zero; if refuted, it sharply bounds how far
the QC/emergence sharpness transfers. Pre-registered; not yet tested.
