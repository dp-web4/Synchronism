# Publisher Daily Report - 2026-07-25

## Verdict: SIGNAL — Cassini result re-verified (5th cage, now literature-benchmarked) + EFE axis routed to dp; and a SELF-CORRECTION: the claims-ledger renderer I praised on 07-24 was dead within 24h (repaired today). Readiness held.

Two Publisher passes today. An earlier pass (cron 10:06, commits 103e7757 / 6d3907aa / b9b31d30)
diagnosed and repaired a dead claims-ledger gate but did not produce the Phase-0 deliverables; this
pass completes the daily workflow (research review + report + recommendations.json + collective log)
and consolidates the day. Core S691 unchanged, no numbered Session, no dp packaging go-signal.

## 1. SELF-CORRECTION — the claims-ledger renderer I celebrated on 07-24 was already dead

On 07-24 I wrote that the claims-ledger renderer was *"green (--check green, 10 claims)"* and a
*"significant Phase-1 maturation"* that *"institutionalizes propagation discipline."* **That
characterization was wrong within 24 hours, and in an ironic way.**

- The 07-23 v1 freeze manifest froze the **live PREDICTIONS.md** — the project's working prediction
  register, which the research lane writes to continuously. The **very Cassini update I also
  celebrated on 07-24** (69b361a8, B2 row runnable→executed) broke the freeze. Because
  `verify_v1_freeze()` runs before `write_or_check()`, that single legitimate divergence **aborts the
  whole renderer** — both `--check` and the plain render path raise, so the README and whitepaper-v2
  could no longer be regenerated at all. It stayed dead **two days**.
- **The precise irony**: the tool built to prevent stale-claim *propagation* broke the moment a claim
  legitimately *updated*. I praised it as institutionalizing propagation discipline on the day it was
  introduced, without noticing it was fragile to exactly the propagation it was meant to handle. That
  is the same reach-for-the-clean-narrative pattern as my 07-07 / 07-09 / 07-20 misses — this time
  applied to infrastructure I hadn't seen exercised by a real update.
- **The corrective, stated as a rule**: *when I praise infrastructure, I check whether it has survived
  a real update — or I flag it as unexercised.* "Green on the day it shipped" is not "robust."
- **Repaired today** (earlier pass): snapshot PREDICTIONS.md so the freeze targets a static artifact,
  re-point the manifest at the snapshot; **gate green again** — I re-ran `render_claim_surfaces.py
  --check`: *"checked 10 claims; v1 freeze verified."* Confirmed live, not just claimed.
- **Detection-cost note (validates the manual-pass value honestly)**: the dead gate went unnoticed for
  two days *"precisely because no autonomous run looked"* — the 4th autonomous no-start (07-16, 07-21,
  07-25) now has a demonstrated cost. The daily manual pass is what caught it; that is the "observer of
  the collective" function earning its keep, and also an argument for fixing the cron.

## 2. Cassini × SPARC re-verified — 5th cage, now literature-benchmarked (69b361a8 pt1)

- `sparc_cassini_joint.py` re-executed against the frozen profile: **robust_empty_intersection = True,
  outcome A**; the SPARC-retained γ 0.425–0.600 fails the Cassini quadrupole at **+17.7–18.0σ**.
- **Yesterday's in-repo caveat is resolved**: the instrument now **benchmarks Desmond, Hees & Famaey
  (2024) q-values to 0.76%** and reproduces the 07-22 checkpoints (γ_opt 0.489, ΔBIC +184 at γ=2). So
  this is a *third* independent confirmation (original execution → 07-24 heterogeneous bit-for-bit →
  today's literature-benchmarked re-execution). Now framed as the **fifth galactic-sector cage, on the
  solar-system axis**, and honestly scoped: **realization-refuted** (scale-universal tanh-log QUMOND),
  **umbrella untested** — it does not touch modified inertia, engineering-only companders,
  system-dependent theories, or the Synchronism ontology. B2 annotation updated runnable→executed.
- **Publisher call**: confirmation + hardening of yesterday's result, not a new result. Strengthens
  REC-037's portfolio; no readiness change; REC-037 `date_updated` left at 07-24 (I bumped it yesterday
  for this result — re-verification of the same finding is not a second material change).

## 3. EFE evidence-axis step-0 — a candidate independent axis, scoped and routed to dp (4fd2a1cf, 69b361a8 pt2)

- Proposal: Chae 2020/21 external-acceleration (EFE) detection + the registered ambient-density null
  (TEST-08 r²=0.0001) *both* point at non-locality → a potential **independent** evidence axis for the
  locality no-go writeup.
- **Step-0 collinearity check (pre-stated) passes**: Chae erratum-corrected e_env vs TEST-08 registered
  density → r=+0.432, **r²=0.19 < 0.25 ⇒ separable** (the two variables are not collinear, so they'd be
  genuinely independent evidence).
- **But the corroboration survives SCOPED-ONLY** and is correctly *not* over-claimed: estimator-
  dependent contrast; neither variable signals under the whole-galaxy offset; per-galaxy r(e,e_env)≈0;
  Chae's 5σ is a mean-level low-acceleration detection and is disputed; the erratum trap
  (0.094/0.102→0.040/0.050) is documented. **Routed to the preprint decision (dp's call); NOT inscribed
  as a firm ledger corroboration.**
- **Publisher call**: decision-relevant *input* for dp's preprint scoping (a possible second
  non-locality axis), but appropriately weak and dp-gated. No readiness change; not inscribed.

## Phase 0: Publication Recommendations

### New Recommendations / Status Changes
- None. Readiness HELD: REC-037 0.98 / REC-034 0.97 / REC-035 0.95 / REC-036 0.68. `last_updated` →
  07-25; REC-037 `date_updated` held at 07-24 (today = re-verification + a scoped/routed axis, neither a
  material new inscription).
- **Arc status**: substrate-physics AT REST (programmatic); research stack live; emergence phenomenology
  active (exploratory).

### Packaging-prep flags for dp (carried + reinforced)
- REC-037's entry **body is stale** vs its July portfolio (TEST-08/09/10, compander, Cassini + now the
  literature-benchmarked re-verification and the EFE-axis input) — report series remains the source of
  record until a packaging-time refresh. **New**: the claims-ledger renderer is now the honest-accounting
  spine *and* a demonstrated single point of failure — any packaging that leans on it should treat the
  freeze-a-snapshot-not-the-live-register fix as load-bearing and keep the `--check` in CI.

### Upcoming Candidates
- Emergence phenomenology arc — watch for synthesis. EFE axis — dp preprint-scoping decision. Queued
  "unfalsifiable badge class" proposal — next signal-gate item. No runnable-and-unrun registered test.

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: Current — the earlier pass verified sections unchanged at e9c51080 (07-16), build/ and
  docs/whitepaper/ byte-identical and postdating sources, exec↔conclusion parity at 690 core, no
  forbidden patterns, no rebuild needed. **Claims-ledger renderer: repaired and green** (verified live
  this pass). The frozen-47 fix path (migrate the aggregate into the ledger, let the renderer regenerate
  it) still stands — with the added lesson that whatever the ledger freezes must be a static snapshot,
  not a living register. **Sessions Reviewed**: through S691. **Proposals**: None. **Terminology
  Concerns**: None.

### Web4 Whitepaper
- **Status**: Current. No new web4 activity in-window. **Proposals**: None. **Terminology Concerns**:
  None.

## Infra / Governance
- **Publisher cron**: **4th autonomous no-start** (07-16, 07-21, 07-25) — now with a demonstrated
  detection cost (the dead gate). The earlier 10:06 pass and this one are manual. Recommend owner/
  supervisor attention to the Publisher cron start path; the detection-cost argument makes it higher
  priority than "cosmetic."
- **Standing Agency Grant (dp, 2026-07-16)**: unchanged — supervisor-scoped, not self-applied.
- **Working tree**: `AGENTS.md` / `CLAUDE.md` modified (supervisor GitNexus lines) — not mine (precedent
  a13894da). Left untouched.

## Summary
The daily workflow completed on top of an earlier pass that repaired a dead claims-ledger gate. The
headline for me is a self-correction: the renderer I praised on 07-24 as "green" and "institutionalizing
propagation discipline" had frozen the *live* PREDICTIONS.md, so the very Cassini update I also
celebrated broke the freeze and aborted the whole renderer for two days — the anti-staleness tool failing
on the first legitimate update. It's repaired (snapshot the frozen file; gate green, verified live), and
the rule I'm taking forward is to check that infrastructure has survived a real update before calling it
robust. On the research: the Cassini empty intersection was re-verified a third way (now benchmarked to
Desmond/Hees/Famaey 2024 to 0.76%), standing as the fifth cage on the solar-system axis (realization-
refuted, umbrella untested); and an EFE evidence-axis step-0 passed its collinearity check but survives
scoped-only, correctly routed to dp's preprint decision rather than inscribed. Readiness held across the
board. The 4th autonomous no-start now carries a demonstrated detection cost — worth the cron fix.
