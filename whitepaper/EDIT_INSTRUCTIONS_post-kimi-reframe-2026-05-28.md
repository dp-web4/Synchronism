# Whitepaper edit instructions — post-Kimi-reframe (2026-05-28)

**Author of instructions**: CBP-Claude (Opus 4.7)
**Trigger**: Inventory → reassessment → plan cycle on the saturation reframe (see `forum/claude/saturation-reframe-*-2026-05-28.md` series) → resurfaced pieces from the whitepaper itself need their framing sharpened and lifted into MRH-active position.
**Anchor frame**: `forum/claude/post-kimi-reframe-execution-plan-2026-05-28.md` (master plan, Stream 1).

---

## Frame the executing agent must hold

Three disciplines apply to every edit below. **If a proposed edit violates any, do not make it — flag instead.**

1. **No 'established' tags during stewardship.** Per dp 2026-05-28: *"we're not at a stage where anything can be honestly claimed as 'established'."* The Appendix-A `✅ Established / ⚠️ Speculative / ❌ Failed` taxonomy is itself a closure-shaped framing. Substantive content gets categorized by **relationship to active MRH** (active / parallel-paths / sidelined / superseded), not by verdict on truth-status.
2. **Findings vs Framings distinction stays explicit** (already in README; preserve it everywhere). Quantitative findings ≠ theoretical positions; conflating them is the failure mode external reviewers most often flag.
3. **Audit findings stand below the reframes — they do not overturn them.** S637 (cosmology → MOND in testable regime), S638 (Curie paramagnet < Landau), S660A (novel-survivor count = 0), S661 (RAR γ=2 refuted ΔBIC=+184), S663B (coherence-language interpretation), S665/S666 (substrate irrotational + dissipative for old `R(I)` rule). New substrate starts from **zero confirmed predictions** and inherits the obligation to produce novel predictions, not the credit of prior reparametrizations.

You're not writing new physics here. You're surfacing and connecting pieces already in the whitepaper — `§4.4` two-level time, `§5.7` complexity-dependent speed limits with the pendulum-in-centrifuge analogy, Appendix A.3 saturation primitive — that had been sidelined to mid-document sections, and updating the framing layer (executive, status taxonomy, open questions) to reflect the current MRH state.

**Do not retroactively edit historical session notes** (`Research/SessionNNN_*.md`) or arc summaries. Those are durable record.

---

## Build hierarchy (read before editing)

```
whitepaper/sections/{NN-name}/{NN-subsection}/*.md   ← edit these (fractal sources)
whitepaper/sections/{NN-name}/index.md               ← if subsection structure changes
                ↓
whitepaper/build.sh  (runs make-md.sh / make-pdf.sh / make-web-clean.sh)
                ↓
whitepaper/build/Synchronism_Whitepaper_Complete.md  ← built artifact (do not edit directly)
whitepaper/build/Synchronism_Whitepaper.pdf
whitepaper/build/web-clean/
                ↓ (copied for publishing)
docs/whitepaper/                                      ← GitHub Pages destination
.governance/build_log.json + file_hashes.json        ← updated by builder
```

**Workflow**: edit fractals → `bash whitepaper/build.sh md` → verify `whitepaper/build/Synchronism_Whitepaper_Complete.md` shows your changes → if good, `bash whitepaper/build.sh` for all formats → commit + push fractals + build artifacts + governance state + docs/whitepaper updates.

---

## Edit task list

### Task 1 — Resolve duplicate `time-slices` subdirectory

`whitepaper/sections/04-fundamental-concepts/` has BOTH `02-time-slices/` and `04-time-slices/`. One is current, one is archived/superseded. Determine which:
- Check `whitepaper/sections/04-fundamental-concepts/index.md` for the canonical section ordering.
- Compare last-modified times.
- Read both — keep the one that matches the §4.4 content in the built whitepaper.

**Action**: if one is superseded, move it to the `archive/` subdirectory rather than deleting. Note the move in commit message. Same duplicate pattern exists for `03-intent-transfer / 05-intent-transfer` and `04-emergence / 06-emergence` — investigate and resolve similarly **only if** doing so doesn't require deep judgment; if unclear, leave them and flag.

### Task 2 — `whitepaper/sections/00-executive-summary/executive_summary.md`

Promote the substrate-reformulation status into the executive framing. Specifically:

- Add a short paragraph (3-5 sentences) near the top of the body, after the existing high-level framing, that states:
  > *"Current MRH status (2026-05-28): the substrate layer is being reformulated. The original Intent transfer rule was found to be 1-DOF scalar diffusion (Session 11) and the original substrate was found irrotational and dissipative (Sessions 665/666). A saturation reframe with independent vector flux **J** and complexity-dependent speed-of-light is currently being worked. Audit findings on prior tracks (cosmology, chemistry, hot superconductor) stand below the reframe — the new substrate inherits the obligation to produce novel predictions, not the credit of prior reparametrizations. See `STATUS.md` for the live MRH inventory."*
- Remove or qualify any executive-summary statements that imply a settled/established substrate.

### Task 3 — `whitepaper/sections/04-fundamental-concepts/{NN-time-slices}/` (canonical one)

The pendulum-in-centrifuge analogy lives in §5.7; the substrate-tick framing lives in §4.4. The two-level ontology (Level 0 = substrate ticks absolute, Level 1+ = pattern-relative frequency comparison) is implied by both but not articulated as a unified frame.

- Add a subsection or paragraph titled **"Two-level time ontology"** that explicitly states:
  - **Level 0 (substrate)**: ticks occur globally at Planck frequency, absolute, computational-clock-of-the-universe. Independent of whether there is local intent transfer.
  - **Level 1+ (patterns)**: all measured time is pattern-relative frequency comparison between embedded oscillator and reference oscillator. This is the structure physics has always used (Newton, GR, QM); Synchronism's framing is structurally equivalent to GR's metric + clock postulate, not weaker.
  - Cross-reference §5.7 for the pendulum-in-centrifuge worked example.
  - State that this dissolves Newton-vs-GR by making both correct at their respective levels.

### Task 4 — `whitepaper/sections/05-quantum-macro/07-speed-limits/`

The pendulum-clock-in-centrifuge analogy is already here. The complexity-dependent speed limit / coherence envelope / time-dilation-as-computational-load material is already here. Two updates:

- Add a clearly-labeled paragraph naming `f(N)` as the **standing open obligation**: *"The reconstruction function `f(N)` — the number of ticks required to stabilize a pattern of complexity *N* in an adjacent cell — has not been derived from the discrete substrate rules. Its derivation, with boundary condition `f(N) → 1` as `N → 0` (minimal complexity = photon), is the single specific path from the framework's complexity-dependent speed structure to quantitative predictions distinguishing it from GR. Candidate experimental discriminators: structured light (OAM photons), entangled photon pairs, neutrino speed, mechanical-vs-atomic clock divergence in strong gravity. See `Research/OPEN_QUESTIONS_*` for the active inventory."*
- Cross-reference §4.4's two-level ontology subsection so the time and speed-limit framings are connected.

### Task 5 — `whitepaper/sections/05-quantum-macro/06-relativity/`

Strengthen the two-level ontology framing in the relativity section. The existing content is structurally compatible; this is a connect-existing-pieces edit, not a rewrite. Add a paragraph linking explicitly to §4.4's two-level ontology and to §5.7's pendulum-in-centrifuge analogy. Emphasize that the framing is structurally equivalent to GR's metric + clock postulate, not weaker.

### Task 6 — `whitepaper/sections/09-appendix-mathematical/mathematical_framework.md` — **status taxonomy retag**

This is the most load-bearing edit. The Appendix A subsections currently carry status tags `✅ Established`, `⚠️ Speculative`, `❌ Failed`. **Replace the status taxonomy with an MRH-relationship taxonomy.**

New taxonomy:
- **`[ACTIVE-MRH]`** — currently in active research focus, content being extended or revised
- **`[PARALLEL-PATHS]`** — in the framework's parallel hypothesis space, not currently in active focus but not abandoned
- **`[SIDELINED]`** — was in active focus, currently not pursued; reasons documented
- **`[SUPERSEDED]`** — replaced by a later formulation in the active or parallel space; pointer to successor

**Concrete edits**:
- Retag A.1, A.2 etc. by their actual MRH-relationship status, not by verdict-on-truth.
- **A.3 Saturation Dynamics**: currently tagged `✅ Established` with the claim *"the full Intent transfer equation in continuum form IS the incompressible Navier-Stokes equation with this variable viscosity — not an analogy, but an exact identification."* Two changes:
  1. Retag as `[ACTIVE-MRH]` (it's the substrate-reformulation focus right now).
  2. Add an inline note flagging the tension: *"This identification is in active inventory tension with Session 11's finding (the rule reduces to 1-DOF scalar diffusion under the maximum principle for parabolic PDEs) and with S665/S666's findings (the substrate is irrotational, curl(v) ≡ 0 for any R(I); and dissipative, first-order ∂I/∂t with decreasing Lyapunov functional). At least one of these three statements requires qualification. Phase 1 simulation work (1D/2D lattice with R(I) = [1 − (I/I_max)^n] sweeping n, plus independent vector flux **J**) will determine which. See `Research/OPEN_QUESTIONS_*` and `forum/claude/saturation-reframe-resurfaced-pieces-mrh-stewardship-2026-05-28.md` §2."*
- A.12 (gravity model) currently tagged `❌ Failed`: retag as `[SIDELINED]` or `[SUPERSEDED]` with pointer to the saturation-reframe / Intent-field substrate work that's now active.
- Update the "Honest Assessment" closing section to remove `✅/⚠️/❌` language. Replace with the MRH-relationship phrasing: *"Sections marked `[ACTIVE-MRH]` are in current research focus. Sections marked `[PARALLEL-PATHS]` are alternative formulations carried in the parallel space. `[SIDELINED]` content is not currently pursued but not abandoned. `[SUPERSEDED]` content points to its successor formulation. The mathematics is a work in progress through stewardship, not a completed foundation."*
- "Open Mathematical Problems" #1 (`What transfer rules generate stable patterns?`) — add a note: *"This is the same question Phase 1 simulation work directly addresses. See post-Kimi-reframe execution plan."*

### Task 7 — `whitepaper/sections/06-implications/04-open-questions/`

Add (or extend) the open-questions content with the post-Kimi consolidated set:
- **OQ-EOS**: stable equation of state replacing `P = I_max − I` (which gives `c_s² = −I_max < 0`, breaks under stability analysis). Polytropic `P ∝ ρ^γ` or similar with `dP/dρ > 0` needed.
- **OQ-Momentum**: discrete-grid derivation of the momentum equation via Chapman-Enskog or finite-volume coarse-graining. Currently asserted at continuum level; required obligation per Kimi 2026-05-28.
- **OQ-fN**: `f(N)` reconstruction function — derive from discrete substrate rules. Boundary condition `f(N) → 1` as `N → 0`.
- **OQ-Oscillation**: demonstrate stable oscillating patterns in 1D/2D lattice simulation. Sweep `(n, I_max, T_ij)`. Test Mechanism A (conservative `J`) vs Mechanism B (CFL violation + saturation feedback → limit cycle).
- **OQ-A3-Tension**: reconcile Appendix A.3's "exact NS identification" claim with Session 11's 1-DOF scalar diffusion finding and S665/S666's irrotational + dissipative proofs.
- **OQ-Discriminators**: quantify predicted deviations from GR/QM for OAM photons (function of ℓ), entangled photon pairs, neutrino propagation, mechanical-vs-atomic clock divergence in strong gravitational fields. All gated on `f(N)` derivation.

Use the MRH-relationship taxonomy (`[ACTIVE-MRH]`, `[PARALLEL-PATHS]`, etc.) on each open question, not `priority: high/medium/low`.

### Task 8 — Cross-references

After the above edits, verify cross-references work:
- §4.4 ↔ §5.7 (two-level time ↔ complexity-dependent speed)
- §4.4 ↔ §5.6 (two-level time ↔ relativity)
- §5.7 ↔ Appendix A.3 (speed limits ↔ saturation primitive)
- Appendix A.3 ↔ §6.4 open questions (A.3-vs-Session-11 tension)

The whitepaper already has a cross-reference style — preserve it. Add new cross-refs where the edits create new connections.

---

## Build + verify

After fractal edits:

```bash
cd whitepaper
bash build.sh md    # build markdown first (cheapest)
```

- Read `whitepaper/build/Synchronism_Whitepaper_Complete.md` and grep for the new content (`grep -n "Two-level time" build/Synchronism_Whitepaper_Complete.md`, etc.) to confirm changes propagated.
- If markdown looks good:
  ```bash
  bash build.sh  # full rebuild (md + pdf + web)
  ```
- If any format fails, investigate before committing. PDF builds may fail on LaTeX issues unrelated to content — if so, commit md + web and flag pdf separately.

## Commit + push

```bash
cd /mnt/c/exe/projects/ai-agents/Synchronism
git add whitepaper/sections/ whitepaper/build/ docs/whitepaper/ .governance/
git status --short   # verify scope
git commit -m "whitepaper: post-Kimi-reframe strategic edit (2026-05-28)

Surfaces and connects pieces already in the whitepaper that had been
sidelined relative to the saturation-reframe arc. Lifts §4.4 two-level
time ontology and §5.7 complexity-dependent speed-limit/c-as-reconstruction-rate
into active MRH framing. Retags Appendix A status taxonomy from
verdict-shaped (✅/⚠️/❌) to MRH-relationship (active/parallel/sidelined/
superseded) per dp's 'nothing established during stewardship' discipline.
Flags Appendix A.3 vs Session 11 vs S665/S666 tension as inventory item.
Adds post-Kimi consolidated open-question set (OQ-EOS, OQ-Momentum, OQ-fN,
OQ-Oscillation, OQ-A3-Tension, OQ-Discriminators).

Audit findings on prior tracks (S637, S638, S660A, S661, S663B, S665/S666)
stand below the reframes; new substrate inherits zero confirmed predictions.

Cycle docs: forum/claude/saturation-reframe-*-2026-05-28.md series.
Build via: bash whitepaper/build.sh.

Co-Authored-By: Claude Opus 4.7 <noreply@anthropic.com>"
git pull --no-rebase
git push
```

## What to flag back rather than execute silently

- Any fractal edit where the existing content's framing is in tension with the post-Kimi discipline AND the right resolution is not obvious — flag in a `whitepaper/POST_KIMI_REFRAME_EDIT_NOTES_2026-05-28.md` file, do not guess.
- Any build failure that you cannot resolve (LaTeX env issues, missing pandoc filter, etc.) — commit what builds, flag the rest.
- The duplicate `time-slices / intent-transfer / emergence` subdirectories: if disambiguating one of them requires reading-comprehension judgment about which version is canonical, flag it rather than guess.

## Out of scope

- Editing the Synchronism site (`/mnt/c/exe/projects/ai-agents/synchronism-site/`) — different repo, handled in Stream 3 of master plan.
- Editing STATUS.md, README.md, AGENTS.md — handled in Stream 2 by the master agent.
- Writing the autonomous-tracks frame doc — handled in Stream 4 by the master agent.
- Editing chemistry validation arcs or hot-SC arcs — those audit findings stand and stay where they are.

— Instructions written 2026-05-28 by CBP-Claude
