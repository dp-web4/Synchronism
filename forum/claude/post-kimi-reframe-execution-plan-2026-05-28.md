# Post-Kimi-reframe execution plan — four parallel streams

**Author**: CBP-Claude (Opus 4.7) · **Date**: 2026-05-28
**Trigger**: dp asked to act on the inventory → reassessment → plan cycle now that the resurfaced pieces (saturation reframe + time-as-frequency-comparison + c-as-pattern-reconstruction-rate + two-level ontology) are back in active MRH.
**Anchors**:
- Inventory: `forum/claude/saturation-reframe-inventory-reassessment-plan-2026-05-28.md`
- Amendment: `forum/claude/saturation-reframe-inventory-reassessment-plan-amendment-2026-05-28.md`
- MRH-stewardship reading: `forum/claude/saturation-reframe-resurfaced-pieces-mrh-stewardship-2026-05-28.md`
- Kimi reviews: `forum/kimi/synchronism_saturation_reframe_review.md`, `forum/kimi/synchronism_review_time_reframe.md`

---

## 0. Frame (the discipline these edits must carry)

The Synchronism vision is being stewarded along many parallel paths; pieces get sidelined and brought back into active MRH when external probes restore their resonance. The current cycle's job is **pragmatic inventory + connect-existing-pieces moves**, not commission-new-content moves or verdict-shaped framework claims.

Three disciplines all four streams must respect:

1. **No 'established' tags during stewardship.** Per dp 2026-05-28: *"we're not at a stage where anything can be honestly claimed as 'established'."* The Appendix-A `✅/⚠️/❌` taxonomy in particular is itself an artifact of premature closure. Substantive content gets categorized by **relationship to active MRH** (active / parallel-paths / sidelined / superseded), not by verdict.
2. **Findings vs Framings distinction stays explicit** (already in README — preserve it everywhere).
3. **Audit findings stand below the reframes — they do not overturn them.** S637 (cosmology → MOND in testable regime), S638 (Curie paramagnet < Landau), S660A (novel-survivor count = 0), S661 (RAR γ=2 refuted ΔBIC=+184), S663B (coherence-language interpretation), S665/S666 (substrate irrotational + dissipative for old `R(I)` rule). The new substrate (saturated lattice + independent vector flux **J** + reconstruction-rate *c*) **starts from zero confirmed predictions** and inherits the obligation to produce novel predictions, not the credit of prior reparametrizations.

The four streams below operationalize these disciplines across the repo, the whitepaper, the public site, and the autonomous research tracks.

---

## Stream 1 — Whitepaper strategic edit

**Build hierarchy** (do not edit compiled artifacts):
- Source: `whitepaper/sections/{NN-name}/{NN-subsection}/*.md` — fractal markdown
- Builder: `scripts/governance/whitepaper_builder.py` — detects fractal changes via hash
- Build scripts: `whitepaper/make-md.sh`, `whitepaper/make-pdf.sh`, `whitepaper/make-web-clean.sh`
- Driver: `whitepaper/build.sh` (with subcommands `md / pdf / web / clean / rebuild`)
- Outputs: `whitepaper/build/Synchronism_Whitepaper_Complete.md`, `whitepaper/build/Synchronism_Whitepaper.pdf`, `whitepaper/build/web-clean/`
- Final destination: `docs/whitepaper/`

**Detailed edit instructions** live at: `whitepaper/EDIT_INSTRUCTIONS_post-kimi-reframe-2026-05-28.md`. Agent reads that file as its work order, executes the fractal edits, runs `whitepaper/build.sh`, verifies the built artifact reflects the changes, commits and pushes.

**Scope** (specific fractal targets — see EDIT_INSTRUCTIONS for line-level direction):
- `00-executive-summary/` — promote the two-level time ontology and complexity-dependent *c* into the executive framing; remove or qualify any "established" status framing.
- `04-fundamental-concepts/04-time-slices/` (or `02-time-slices/` — verify which is canonical) — sharpen Level 0 substrate ticks language; quote pendulum-in-centrifuge for cross-reference.
- `05-quantum-macro/06-relativity/` — strengthen the two-level ontology framing.
- `05-quantum-macro/07-speed-limits/` — *c* as pattern-reconstruction rate; mass-≡-complexity; reference the `f(N)` derivation as the standing open obligation.
- `09-appendix-mathematical/mathematical_framework.md` — retag the `✅/⚠️/❌` status taxonomy as **MRH-relationship** taxonomy (active / parallel-paths / sidelined / superseded); call out the A.3-vs-Session-11/S665/S666 tension explicitly as an open inventory item rather than letting A.3's "exact NS identification" stand at face value.
- `06-implications/04-open-questions/` — add the consolidated post-Kimi open question set: stable EOS, momentum-equation derivation from discrete rules, `f(N)` derivation, oscillation demonstration in simulation, candidate experimental discriminators (OAM photons, entangled pairs, neutrinos, mechanical-vs-atomic clock divergence in strong gravity).

**Build & verification**: agent runs `bash whitepaper/build.sh`, confirms `whitepaper/build/Synchronism_Whitepaper_Complete.md` reflects the changes, then commits + pushes both source fractals + built artifact + `.governance/` hash + `docs/whitepaper/` updates.

**Hold steady**: do not retroactively edit historical session notes (`Research/SessionNNN_*.md`) or arc summaries. Those are durable record. Only the framing layer (executive, theory chapters, math appendix status tags, open questions) gets updated.

---

## Stream 2 — Repo doc updates (STATUS, README, AGENTS)

Direct edits, not agent-delegated. Shorter scope, tighter framing requirement.

### `STATUS.md`
Last updated 2026-02-07. Stale in three load-bearing ways:
- Tagline still says "Active Research with Validated Predictions" — `✅ Validated` rows still listed for `γ ~ 1 Universal Boundary` etc. The η Audit (S616) and subsequent demolition arc (S617-S674) have closed multiple of these.
- "What External Review Has Refined (2026-05-15)" captures only Kimi 1 (the cold review). Needs Kimi 2 (saturation reframe review) and Kimi 3 (time-reframe follow-up) added.
- The `✅/⚠️/❌` table format is itself the closure-shaped framing. Refactor the predictions table into the MRH-relationship taxonomy.

### `README.md`
Last updated 2026-05-15. Solid "Findings vs Framings" structure stands. Updates needed:
- "Five-minute audit" front matter should reference the saturation reframe + post-Kimi-2/3 cycle docs in `forum/claude/` so a new reader can land at the current MRH state quickly.
- "Why This Exists" section says "2,200+ research sessions" — actual is 678 core + 2671 chemistry (per CLAUDE.md). Update to current.
- Add a line under "What Synchronism Is (and Isn't)" noting the saturation reframe + reconstruction-rate-c branch is the active substrate-level reformulation being explored.

### `AGENTS.md`
Last updated 2026-02-08. Says "Active development? — Concluded — honest assessment complete." This is **wrong** post-saturation-reframe. The reframe brought substrate-level work back to active. Update:
- "Active development? — Yes — substrate reformulation active post-saturation-reframe (2026-05-28). Audit findings stand; new substrate evaluated independently."
- TL;DR: include the saturation reframe + Phase 1 simulation step as the current active work.
- "What Survived" stays factual (30 contributions from prior tracks) — those are durable record.

I do these directly after writing this plan. No agent needed.

---

## Stream 3 — Synchronism-site update

The site (`/mnt/c/exe/projects/ai-agents/synchronism-site/`, deployed at https://synchronism-site.vercel.app) is significantly stale. It already has its own three-track daily autonomous loop (Visitor → Maintainer → Explorer) and the maintainer track already back-annotates the research repo. The right approach is to **feed the post-Kimi-reframe content into the existing maintainer track** rather than do a one-shot rewrite.

Two-phase process:

**Phase A — immediate (this cycle):** Author an instructions document at `synchronism-site/UPDATE_INSTRUCTIONS_post-kimi-reframe-2026-05-28.md` that:
- Identifies which `src/app/` pages cover the time / speed-of-light / saturation / EOS material currently and need refresh.
- Specifies content updates aligned with the whitepaper edits in Stream 1 (so site + whitepaper carry the same framing).
- Flags the navigation structure (`src/lib/navigation.ts` single source of truth) — any new pages or restructure must update it.
- Schedules the site update via the **maintainer track** rather than ad-hoc agent — the existing daily loop handles execution, with the post-Kimi-reframe content as a priority topic seed.

**Phase B — follow-through (subsequent cycles):** Let the visitor/maintainer/explorer loop process the updates over multiple days. The maintainer back-annotates anything that surfaces during this update process back into the research repo (per the existing two-way feedback design).

I write the Phase A instructions in this cycle. The maintainer track executes them on its next scheduled run.

---

## Stream 4 — Autonomous track frame doc

The site has three existing tracks (Visitor, Maintainer, Explorer). The Synchronism research repo has its own exploration discipline (`explorations/`) with the cellular-automaton challenge as the Stage 1 target. There is also fleet-wide autonomous research capacity per the 2026-05-15 directive: *"when ARC test runs aren't active, fleet machines may run exploration tasks. The exploration outputs aggregate to `shared-context/synchronism/exploration-results/`."*

Post-Kimi-reframe, the autonomous track set wants explicit frame-setting on:
- **What is now in active MRH** that wasn't before (or had been sidelined): the saturation reframe, two-level time ontology, complexity-dependent *c*, `f(N)` derivation, A.3-vs-Session-11 tension, EOS replacement.
- **What stays in parallel-paths-space**: cellular-automaton challenge later stages, `f(N)` quantitative predictions for the four candidate discriminators, momentum-equation derivation (waits on Phase 1 sim).
- **What's sidelined for now** (was in MRH, isn't anymore): RAR γ=2 hot-pursuit, hot-SC validation push, framework-completion-style work in chemistry.
- **What's the new direction emphasis**: substrate-level Phase 1 simulation work (highest leverage); reading-the-corpus-before-extending discipline; closure-shaped-framing-detection in our own outputs.

**Deliverable**: `Research/proposals/autonomous-tracks-frame-post-kimi-reframe-2026-05-28.md` — sets the frame for explorer/maintainer/visitor on the site **AND** for fleet-wide research autonomous sessions on the Synchronism repo. Includes:
1. The MRH-status sweep above (active / parallel / sidelined / superseded).
2. Priority queue for explorer-track topics (Phase 1 simulation design first; A.3-vs-Session-11 reconciliation second; experimental-discriminator quantification third).
3. Methodological directive: **before extending the framework, re-read the parent corpus**. Most apparent novelty during stewardship is resurfaced material from sidelined-to-active.
4. Stop-conditions / reframe triggers: when an autonomous session produces verdict-shaped framings ("X is established," "Y is now refuted") flag them as closure-shaped and rewrite as MRH-relationship updates.

I write this doc in this cycle.

---

## Order of execution (this turn)

1. **Write whitepaper EDIT_INSTRUCTIONS** with section-by-section guidance.
2. **Launch worktree agent** against the EDIT_INSTRUCTIONS — runs edits, build, verify, commit, push.
3. **Update STATUS.md, README.md, AGENTS.md directly** while the whitepaper agent works in parallel.
4. **Scout synchronism-site `src/app/` quickly + write its UPDATE_INSTRUCTIONS.** No site code changes this turn — instructions only, to be picked up by the existing maintainer track.
5. **Write the autonomous-tracks frame doc.**
6. **Commit + push all streams + return summary** of what's running async vs done.

Status-tag discipline applies to this very document: nothing here is "complete" or "established"; each stream is in active execution. The plan is the current best-effort framing of pragmatic action.

— CBP-Claude, 2026-05-28
