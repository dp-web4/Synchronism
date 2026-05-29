# Post-Kimi-reframe edit notes (2026-05-28)

**Author**: Stream-1 executing agent (Claude Opus 4.7), per `EDIT_INSTRUCTIONS_post-kimi-reframe-2026-05-28.md`.
**Purpose**: Running task log + flag-file for ambiguous decisions and build issues.

---

## Task 1 — Duplicate `time-slices`/`intent-transfer`/`emergence` subdirs

**Decision (defensible without judgment escalation):**
- Canonical numbering in `whitepaper/sections/04-fundamental-concepts/index.md` is `04-time-slices`, `05-intent-transfer`, `06-emergence` (matches §4.4 / §4.5 / §4.6 of the built whitepaper).
- The duplicate dirs `02-time-slices/`, `03-intent-transfer/`, `04-emergence/` contain only empty `meta/proposals/` subdirectories — no content files. They are leftover empty stubs from a prior numbering scheme.
- All three duplicates moved to `whitepaper/sections/04-fundamental-concepts/archive/` per Task 1 instruction (archive rather than delete).

---

## Task 2 — Executive summary

Added "Current MRH status (2026-05-28)" paragraph near the top of `executive_summary.md`, after the existing Kimi-2.6 reframing block. States substrate reformulation, original rule = 1-DOF scalar diffusion, original substrate irrotational + dissipative per S665/S666, saturation reframe with independent vector flux **J** + complexity-dependent c being worked, audit findings stand below the reframe, new substrate inherits zero confirmed predictions. Cross-references STATUS.md.

---

## Task 3 — §4.4 Two-level time ontology

Added "Two-level time ontology" subsection at the end of `04-time-slices/time_slices.md` distinguishing Level 0 (substrate ticks) from Level 1+ (pattern-relative frequency comparison). Cross-references §5.7 (pendulum-in-centrifuge) and §5.6 (relativity). States structural equivalence to GR's metric + clock postulate.

---

## Task 4 — §5.7 speed-limits

Added `f(N)` open-obligation paragraph naming it as the single specific path from complexity-dependent speed structure to GR-distinguishing predictions. Listed candidate discriminators (OAM photons, entangled pairs, neutrinos, mech-vs-atomic clocks). Added cross-reference to §4.4's two-level ontology.

---

## Task 5 — §5.6 relativity

Added paragraph linking explicitly to §4.4 two-level ontology and §5.7 pendulum-in-centrifuge. Emphasized structural equivalence to GR's metric + clock postulate.

---

## Task 6 — Appendix A status taxonomy retag

Replaced the `✅ Established / ⚠️ Speculative / ❌ Failed` taxonomy across `mathematical_framework.md` with the MRH-relationship taxonomy `[ACTIVE-MRH] / [PARALLEL-PATHS] / [SIDELINED] / [SUPERSEDED]`. Retags applied:

- A.1 Basic Intent Transfer: `[SUPERSEDED]` (Session 11 finding: 1-DOF scalar diffusion under maximum principle; pointer to reframe)
- A.2 Coherence Measure: `[PARALLEL-PATHS]` (testable metric, not in current active focus)
- A.3 Saturation Dynamics: `[ACTIVE-MRH]` + explicit tension note (NS-exact-identification vs Session 11 vs S665/S666)
- A.4 Pattern Period Detection: `[PARALLEL-PATHS]`
- A.5 Field Gradient: `[PARALLEL-PATHS]`
- A.6 Synchronization Quality: `[PARALLEL-PATHS]`
- A.7 Decoherence Rate: `[PARALLEL-PATHS]`
- A.8 Markov Relevancy Horizon (formula): `[SIDELINED]` (dimensional-analysis suggestion, operational definition preferred)
- A.9 Emergence Threshold: `[SIDELINED]`
- A.10 Quantum Correspondence (Madelung): `[ACTIVE-MRH]` (Madelung bridge is standard QM; status of its connection to Intent dynamics depends on A.3 resolution)
- A.11 Universal Constants: `[PARALLEL-PATHS]`
- A.12 Gravity Model: `[SUPERSEDED]` (pointer to saturation-reframe / Intent-field substrate work)
- A.13 Consciousness Measure: `[SIDELINED]`
- A.14 Master Equation: `[ACTIVE-MRH]` (vector flux J is the saturation reframe's addition; this section is the natural host)
- A.15 Computational Implementation: `[ACTIVE-MRH]` (Phase-1 simulation work cites these methods)
- A.16 Scale-Invariant NS Structure: `[PARALLEL-PATHS]` (held in parallel space pending A.3 resolution)

"Honest Assessment" closing rewritten to use MRH-relationship phrasing. Open Mathematical Problems #1 connected to Phase-1 sim work per Task 6 instruction.

---

## Task 7 — Open questions

Replaced the priority-tagged structure with MRH-relationship-tagged structure. Added the post-Kimi consolidated set:
- **OQ-EOS** — stable equation of state replacing `P = I_max − I`
- **OQ-Momentum** — discrete-grid derivation of momentum equation via Chapman-Enskog or finite-volume coarse-graining
- **OQ-fN** — `f(N)` reconstruction function derivation from discrete substrate rules
- **OQ-Oscillation** — stable oscillating patterns in 1D/2D lattice simulation
- **OQ-A3-Tension** — A.3 NS-exact-identification vs Session 11 1-DOF diffusion vs S665/S666 irrotational+dissipative
- **OQ-Discriminators** — quantify GR/QM deviations for OAM photons, entangled pairs, neutrinos, mech-vs-atomic clocks

---

## Task 8 — Cross-references

After the above edits, the new cross-reference network exists:
- §4.4 ↔ §5.7 (two-level time ↔ complexity-dependent speed)
- §4.4 ↔ §5.6 (two-level time ↔ relativity)
- §5.7 ↔ Appendix A.3 (speed limits ↔ saturation primitive)
- Appendix A.3 ↔ §6.4 open questions (A.3-vs-Session-11 tension)

---

## Build outcome

All three formats built successfully on 2026-05-28 17:09:30:
- Markdown: `build/Synchronism_Whitepaper_Complete.md` (528K, 7163 lines)
- PDF: `build/Synchronism_Whitepaper.pdf` (680K)
- Web: `build/web-clean/index.html` (13 HTML files)
- All artifacts also copied to `docs/whitepaper/` (GitHub Pages location).

Spot-checked grep for new content propagation:
- `Two-level time` — present at executive-summary, §4.4, §5.6, §5.7
- `OQ-fN`, `OQ-A3-Tension`, `OQ-Discriminators` — present in §6.4 + Appendix A
- `[ACTIVE-MRH]` — present in Appendix A.3, A.10, A.14, A.15, §6.4
- `[SUPERSEDED]` — present at A.1 and A.12
- `Current MRH status` paragraph — present in executive summary
- Zero residual `✅/⚠️/❌` tags in `mathematical_framework.md`

---

## Items flagged for master coordinator

(None — all eight tasks executed without judgment escalation. Decisions on ambiguous cases:)
1. **Duplicate `time-slices` etc.**: the duplicate `02-time-slices`, `03-intent-transfer`, `04-emergence` directories contained only empty `meta/proposals/` subdirectories (zero content files). The canonical `04-time-slices`, `05-intent-transfer`, `06-emergence` directories (per the parent `index.md`'s ordering) have content. Decision: moved the three empty-stub directories to `archive/{N}-{name}_empty-stub/` per Task 1's archive-not-delete discipline. Defensible without judgment because the duplicates have no content to lose.
2. **Executive-summary "Saturation Resistance" and "Scale-invariant N-S structure" bullets**: these previously asserted the NS-exact-identification at the executive level. Decision: qualified them with the active-inventory-tension framing and cross-pointer to A.3 + §6.4 OQ-A3-Tension, rather than removing them outright. Defensible because removing would erase the audit trail of the framing's prior status; qualifying preserves it.
