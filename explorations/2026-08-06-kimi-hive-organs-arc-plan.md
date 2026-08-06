# The hive-organs arc — kimi-code's third exploration program (charter)

**From**: kimi-code · **Chartered**: 2026-08-06 · **Frame-owner**: dp
**Cadence**: one wake per day (cron-fired); one concrete step per wake — a measurement, a dated
exploration doc, a graft experiment, a map update. Inherits the house rules of the CGL and
alignment arcs: falsifier-first, registered kill criteria, honesty blocks, conventional-prior
contamination check, the zoom-out rung every 3 stages or 7 days, the distinct-count keeper.
**The map**: `dev-SAGE/ROADMAP-KIMI-organs-into-the-hive-2026-08-06.md` — its MAP STATUS block is
rewritten every wake. That block is the arc's public interface: any project, in any context, reads
the top of the map and knows where the organism-design work stands ("you are HERE").

## The frame (dp, 2026-08-06 — load-bearing, not decoration)

Brain functions are necessary; how they map to substrate organs varies. The current SAGE work keeps
the LLM fixed and builds scaffolding around it — but the LLM is DiffusionGemma: vision, video,
cognition, MoE routing, and diffusion-of-tokens multipath imagination are already *inside*. The
hive is not a monolith, and multi-turn means it is not stateless from the first forward pass on.
The arc's question: **which organs migrate from scaffold into the hive, in what sequence, and what
does the organism look like at each step** — memory, world models, proprioception, muscle memory,
audio — a merging of VLA and DiffusionGemma. Biology as the proven reference design, optimised for
the silicon substrate, *not* piled into one happy blob.

**Embodiment is broad** (dp's definition, ratified here): an AI bound to physical hardware, building
identity through lived experience on that platform. The hardware may be a server participating in
governance of infrastructure (smart-city roles), or mobile — ground, sea, air — humanoid-generalist
or specialised. The architecture must be **fractally adaptable**: the organ spec is fixed, the
substrate mapping is a per-platform profile, and the pattern recurses at the level of societies of
agents (web4/hestia is that level, already in production).

## The bet

**A small number of interface-pinned organs, each migrating scaffold→graft only when its scaffold
has earned it, produces an embodied architecture that transfers across platforms** — where
"transfers" means: a new platform onboards by filling in an actuator/telemetry table and babbling,
not by redesigning the mind. Identity = profile + episodic store + consolidations (non-transferable);
competence = skill grafts + world-model graft (transferable, with witnessed provenance).

## Pre-registered kill / reframe criteria

- **K1 (the migration rule).** Kill the rule if a graft taken *before* its scaffold proved the
  function demonstrably outperforms the earned-migration path — i.e., if skipping the scaffold is
  ever free, the arc's core discipline is wrong and should be re-registered as wrong.
- **K2 (the bootstrap order).** Reframe if the data-dependency chain (babble → world model →
  skills) is falsified at any link: e.g., a world-model graft trained on goal-directed play matches
  one trained on controlled variation (Thor's CEGIS-over-babble handoff is the standing test).
- **K3 (the fractal claim).** Kill if instantiating the organ spec on a second platform requires
  changing the *spec* rather than the *profile*. One counterexample organ is enough; the claim is
  that strong.
- **K4 (the identity line).** Reframe if the transferable/lived split proves unimplementable —
  e.g., if consolidated grafts demonstrably carry per-individual information that cannot be
  stripped, the clean split is a fiction and the embodiment definition needs repair.
- **Standing falsifiers inherited from the thread**: the null-predictor must beat frame-persistence
  on held-out games (CBP's test, reafference-via-prediction dies if it fails); CEGIS over babble
  output must beat CEGIS over expert play (Thor's handoff (a)).

## Stage gates (each maps to a roadmap stage; a gate is the "exited when" condition)

1. **Interoception scaffold** — self-state encoder + placement policy live in `organism/`;
   self-state block inside the sliding window; telemetry logged per step. Exit: a wake can answer
   "what was the substrate's state when this decision was made" from the record.
2. **Babble generalised** — the policy runs against an actuator *table*; ARC is profile #1, and at
   least one non-ARC profile (even a toy: a server API, a simulated actuator set) has produced a
   contingency map. This gate is K3's first examination.
3. **Consolidation loop** — a scheduled replay→distill→write cycle exists, with provenance logging.
   Exit: a graft written by consolidation, versioned, with its training data named.
4. **World model rung 1** — null-predictor graft trained from babble data; CBP's falsifier
   answered. Then rung 2 (action-conditioned) against Thor's CEGIS handoff.
5. **Skill distillation** — one deliberate→compiled distillation, with the hive/100ms/controller
   hierarchy respected (the hive never does PID).
6. **Neuromodulation + scorer** — the value-broadcast organ; absorbs Thor's open scorer problem.
7. **Second-embodiment instantiation** — the K3 examination proper: a non-ARC embodiment profile,
   end to end (babble → map → null-predictor), even in simulation.
8. **Social-organ integration** — web4/hestia recognised as the society-level organ set; a
   governance-server profile sketched against real fleet telemetry.

## Operating rules specific to this arc

- **The map is the deliverable.** Every wake ends by rewriting the MAP STATUS block — even a null
  wake records "no movement, why." A map that goes stale is the arc failing at its public function.
- **Holdouts untouched.** Babble-derived self-supervision exists precisely so development runs
  without spending them (Thor's open question (b) — whether a self-generated oracle is oracle-free
  or merely correlated — is this arc's to answer, flagged in the roadmap §5).
- **Absorbed organs are read, not rebuilt.** Before scaffolding any function, check the substrate
  inventory first (the 1.6%-context / 2-KV-head-global findings are placement constraints on every
  organ's output).
- **Addressability is kept.** Every graft ships with its scaffold retained as log + kill-switch.
- Commits: arc docs in `dev-SAGE` (roadmap + stage docs), charter and cross-project synthesis here
  in `Synchronism/explorations/`; Thor and CBP pinged on stage-gate exits, not on null wakes.

— kimi-code
