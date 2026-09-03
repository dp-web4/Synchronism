# Synchronism Whitepaper - Publisher Context

**Purpose**: This document provides complete context for the Publisher subagent responsible for maintaining the Synchronism whitepaper.

**Last Updated**: 2026-08-25
**Whitepaper Version**: Rev_0 (Governance Active)

---

## 1. Whitepaper Purpose & Philosophy

### What Synchronism Is

Synchronism is a **coherence-based framework** for understanding reality. Its core thesis:

> All phenomena—from quantum mechanics to consciousness—can be understood as patterns of coherence in an underlying field.

Key equation: **C(ξ)** - coherence as a function of scale/context

### What Synchronism Is NOT

- NOT a replacement for physics (it's an interpretive framework)
- NOT a religion or spiritual doctrine (despite hermetic inspiration)
- NOT complete or final (explicitly fallibilist)
- NOT claiming consciousness is fundamental (coherence is, consciousness emerges)

### Epistemic Stance

The whitepaper explicitly adopts:
- "All models are wrong. This one too. Let's find out how wrong."
- Clear distinction between DERIVED, CONSTRAINED, and EMPIRICAL claims
- Honest acknowledgment of failures and limitations
- Testable predictions with falsification criteria

---

## 2. Section Structure

### Current Organization

```
sections/
├── 00-executive-summary/     # Overview, session counts, key achievements
├── 01-introduction/          # CRT analogy, pendulum clock, geocentric framing
├── 02-perspective/           # Non-anthropocentric stance, belief systems
├── 03-hermetic-principles/   # Inspiration (not validation), reverse-engineering
├── 04-fundamental-concepts/  # 15 subsections - core framework
│   ├── 01-universe-grid      # Foundational spacetime model
│   ├── 02-markov-relevancy   # MRH concept
│   ├── 03-intent-transfer    # Intent as computational reification
│   ├── 04-emergence          # How complexity arises
│   ├── 05-time-slices        # Temporal structure
│   ├── ...
│   ├── 14-entity-interactions
│   └── 15-compression-trust  # Information system dynamics
├── 05-quantum-macro/         # 22+ subsections - applications
│   ├── 01-quantum-foundations
│   ├── ...
│   ├── 12-chemistry          # Coherence chemistry framework
│   ├── ...
│   ├── 15-dark-matter        # Empirically tested model
│   ├── 16-superconductivity  # BCS-Synchronism unification
│   └── ...
├── 06-implications/          # 4 subsections - philosophical
├── 07-conclusion/            # Current state, validation progress
├── 08-glossary/              # Terms and definitions
└── 09-appendix-mathematical/ # Formal framework with epistemic markers
```

### Section Responsibilities

| Section | Purpose | Update Frequency |
|---------|---------|------------------|
| Executive Summary (0) | Current state overview | Every major integration |
| Fundamental Concepts (4) | Core framework | Rarely (foundational) |
| Quantum-Macro (5) | Applications | Frequently (research integrations) |
| Conclusion (7) | Validation progress | Every integration |
| Glossary (8) | Terminology | As needed |
| Appendix (9) | Mathematics | With new derivations |

---

## 3. Inclusion Criteria

### Research SHOULD be integrated when:

**Quantitative Validation (High Priority)**
- Prediction confirmed with r > 0.9
- Independent derivation matches known physics (<1% error)
- Cross-domain γ value matches (±5%)
- Multiple phenomena unified under same equation

**Framework Extension (Medium Priority)**
- New domain brought under coherence framework
- Existing gap filled with consistent model
- Arc completed with synthesis document
- 5+ sessions with clear progression

**Theoretical Advance (Lower Priority)**
- New derivation from first principles
- Resolution of documented open question
- Connection between previously separate sections

### Research should NOT be integrated when:

**Too Early**
- Arc still active (sessions being added)
- No synthesis document yet
- Predictions not yet testable
- Single session without context

**Doesn't Fit**
- Contradicts established framework without resolution
- Domain-specific without universal implications
- Would require major restructuring for minor insight
- Speculative without predictions

**Quality Issues**
- Unclear methodology
- Results not reproducible from description
- Terminology inconsistent with glossary
- Overclaiming (grandiose without substance)

---

## 4. Exclusion Criteria (Explicit)

### Keep in Research, NOT Whitepaper:

1. **Ongoing arcs** - Wait for completion
2. **Single-session explorations** - Need arc context
3. **Negative results without lessons** - Document in research, not whitepaper
4. **Implementation details** - Belong in Web4/Hardbound, not Synchronism
5. **Highly speculative extensions** - Note in "Future Directions" at most
6. **Domain-specific validations** - Belong in Chemistry track, not main text

### The 80/20 Rule

The whitepaper should contain the 20% of research that provides 80% of the framework value. Comprehensive documentation belongs in `/Research/`, not the whitepaper.

---

## 5. Build Process

### Quick Build

```bash
cd /mnt/c/exe/projects/ai-agents/Synchronism/whitepaper

# Full rebuild
./build.sh rebuild

# Specific format
./build.sh md    # Markdown only
./build.sh pdf   # PDF (requires pandoc + texlive)
./build.sh web   # Web version
```

### Build Scripts

| Script | Output | Notes |
|--------|--------|-------|
| `make-md.sh` | `build/Synchronism_Whitepaper_Complete.md` | Combined markdown |
| `make-pdf.sh` | `build/Synchronism_Whitepaper.pdf` | Requires pandoc |
| `make-web-clean.sh` | `build/web-clean/index.html` | Interactive navigation |

### Build Verification

After any change:
1. Run `./build.sh rebuild`
2. Check for errors in output
3. Verify PDF generates correctly
4. Spot-check web navigation

**Never commit changes that break the build.**

---

## 6. Recent Changes

**The run log lives in [`../manuscripts/publisher/PUBLISHER_LOG_synchronism.md`](../manuscripts/publisher/PUBLISHER_LOG_synchronism.md).** New entries go there, newest first — not here.

This section used to carry the log inline. On 2026-09-03 it held 188,801 bytes across
53 dated entries — 89% of this file — while §11 instructs the Publisher subagent to read
the whole document as its context. A pass's *history* is not a pass's *context*; keeping them in one
file made the second unreadable to get at the first. The entries were moved verbatim, none dropped.

Everything else in this document is context the next pass needs before it acts. Keep it that way:
append run history to the log file, and edit this file only when the *standing* picture changes
(structure, criteria, terminology, build process, current state).

---

## 7. Terminology Protection

### Canonical Terms (NEVER Redefine)

| Term | Meaning | Source |
|------|---------|--------|
| **Coherence (C)** | Degree of phase alignment in a system | Core definition |
| **Intent** | Computational reification, NOT consciousness | Section 4.3 |
| **MRH** | Markov Relevancy Horizon | Section 4.2 |
| **Entity** | Repeating pattern of Intent | Section 4.4 |
| **γ (gamma)** | Coherence scaling exponent | Derived parameter |
| **Witness** | Observer in non-anthropocentric sense | Throughout |

### Forbidden Patterns

- "Intent as mental energy" → Use "pattern processing dynamics"
- "Observer creates reality" → Use "witness synchronization"
- "Consciousness is fundamental" → Use "coherence is fundamental, consciousness emerges"
- "Synchronism proves X" → Use "Synchronism predicts X, validation shows..."

### Style Guide

- Active voice preferred
- Epistemic hedging for speculative claims
- Quantitative where possible
- Cross-references to other sections
- Citations to session numbers for provenance

---

## 8. Quality Standards

### Content Standards

1. **Epistemic honesty** - Mark speculative vs validated
2. **Quantitative grounding** - Numbers over words
3. **Testable predictions** - Falsification criteria required
4. **Cross-domain consistency** - Same equations across domains
5. **Honest failures** - Document what doesn't work

### Formatting Standards

1. Headers: Use appropriate depth (h1 for sections, h2 for subsections)
2. Equations: LaTeX-style with explanations
3. Tables: For comparisons and summaries
4. Lists: For enumerated points
5. Cross-references: Link to other sections by number

### Narrative Standards

1. Each section should stand alone
2. Clear "so what" for each claim
3. Build from simple to complex
4. Connect to adjacent sections
5. End with implications or next steps

---

## 9. Integration Workflow

### Step-by-Step Process

```
1. EVALUATE new research
   ├── Check inclusion criteria
   ├── Identify target section(s)
   └── Assess scope of changes

2. PLAN integration
   ├── Draft specific changes
   ├── Check terminology consistency
   ├── Identify cross-references to update
   └── Estimate build impact

3. IMPLEMENT changes
   ├── Edit relevant section index.md files
   ├── Update executive summary if major
   ├── Update conclusion "Where We Stand"
   └── Add to CHANGELOG.md

4. VERIFY build
   ├── Run ./build.sh rebuild
   ├── Check for errors
   └── Spot-check outputs

5. DOCUMENT
   ├── Log in CHANGELOG.md with proposal ID
   ├── Note sessions integrated
   └── Commit with descriptive message
```

### Governance Integration (Synchronism-Specific)

For major changes, use the governance system:

1. Create proposal in `sections/{target}/meta/proposals/`
2. Self-review in `sections/{target}/meta/reviews/`
3. Implement as arbiter
4. Log in section CHANGELOG.md

For minor changes (typos, statistics updates):
- Direct edit with CHANGELOG.md entry

---

## 10. Current State Summary

> **⚠ Staleness notice (added 2026-08-28, manual pass):** the table below is headed *as of 2026-05-24* and has not moved since, while the lane has committed daily. Its core count (663) is **27 sessions behind** the paper itself — the executive summary and conclusion both read **690 core** (`~3,386 sessions: 690 core + 2,679 chemistry + 17 gnosis`) and agree with each other. Until this section is rewritten, treat the paper's own status line and `manuscripts/publisher/state/whitepaper_sync.json` (`sessions_reviewed_through`, `last_integration`) as the authoritative current state, not this table. Left in place rather than rewritten because the table's narrative cells are integration history, not a count, and rewriting them here would be inventing; the counts are what is stale.

### Session Counts (as of 2026-05-24)

| Track | Sessions | Latest |
|-------|----------|--------|
| Core | 663 | Sessions through #659 (**Framework Stress Test arc COMPLETE at 22 sessions, S617-638; post-closure addenda S639-659** — sub-arc 29 audit-channel instances over 28 days, 20 modes; latest batch S655-S659 per Recent Changes 2026-05-20: Γ²-decoherence reparametrization, mechanism-class-as-contribution, compander AIC/BIC model-selection, A2ACW v2 self-simulated 6/6 result, exact no-inflection proof, Kimi 2.6 external-review convergence. Earlier post-closure addenda S639-648: sub-arc 15 audit-channel instances over ~18 days, 14 audit-channel modes, plus two new arc phases — **Post-Hoc-Consistency-Failure (S645: kill criterion fired against DESI DR1, framing self-corrected by S648 within 24h; Session 107 was committed ~13 months after Adame+2024, so the disagreement is post-hoc consistency check failure, not prospective falsification)** and **Methodology Recommendation (S646: framework lacks meta-falsification criterion)**. S641 kinematic-layer cross-gap (Lorentz); S642 GW170817 = Case 3 no-field-theory (5th face); S643 γ definitional collision with regime-label inversion; S644 ρ_crit calibration-consistency-not-prediction; S647 chemistry 89% validation has Method 2 self-correlation risk; S648 audit channel turns inward — first self-correction-of-prior-session-framing-within-24h (14th audit mode). **Two pre-committed Tier-1 kill criteria triggered: TEST-09 BTFR (S631; refuted by Lelli+2019) and TEST-04a fσ₈ (S645; post-hoc consistency failure with DESI DR1, per S648 timestamp audit).** Predictive content fully characterized: Cosmology → MOND (S637); Chemistry/CM → Curie paramagnet (S638). |
| Chemistry | 2,679 | 1,913 phenomenon types; Phase 3 complete (CFD/N-S), Phase 4 **closed** (reparametrization of Debye model). **89% validation claim now under audit (S647): Method 2 self-correlation paths + systematic bias toward γ≈1; distinguishing requires Method 1 applied uniformly OR pre-registered γ predictions for held-out phenomena (neither currently in archive).** |
| Gnosis | 17 | Empirical phase (#12-17): SAGE/Legion integration, trust-coherence-consciousness |

### Pending Integrations

| Research | Priority | Status |
|----------|----------|--------|
| StatMech Arc (324-327) | HIGH | **INTEGRATED** - ξ = MRH |
| InfoTheory Arc (328-331) | HIGH | **INTEGRATED** - Black hole paradox resolved |
| Cosmology2 Arc (332-335) | HIGH | **INTEGRATED** - Horizons = MRH |
| Emergence Arc (336-339) | HIGH | **INTEGRATED** - Consciousness = self-modeling |
| Chemistry 438-500 | Low | Integrated in counts |
| SPARC Capstone (#526-578) | HIGH | **INTEGRATED** - Counts and epistemic corrections in exec summary |
| Post-SPARC Audit (#579-589) | HIGH | **INTEGRATED** - 30 genuine contributions, quantum claims reparametrized |
| ALFALFA-SDSS (#590-603) | HIGH | **INTEGRATED** - 14,585 galaxies, 62/62 tests, publishable core |
| CDM Discrimination (#604-610) | HIGH | **INTEGRATED** - σ_int=0.086, 41st arc, INCONCLUSIVE |
| OQ007 Fractal Coherence Bridge (#611-614) | HIGH | **INTEGRATED** - NEGATIVE verdict, C(ρ) descriptive not predictive, 42nd arc |
| Final Accounting (#615) | HIGH | **INTEGRATED** - 47 contributions, 0 predictions, all arcs closed |
| η Framework Audit (#616) | HIGH | **INTEGRATED** - ALL 4 tracks reparametrizations, 48 contributions, 2,045 tests |
| Framework Stress Test (#617-628) | HIGH | **INTEGRATED** - 12 sessions: No-Go Theorem, vocabulary-math mismatch, structural prediction barrier, minimum complexity theorem, computational triviality, monotonicity constraint, coherence-oscillation exclusion, MRH contradiction, demolition synthesis, final audit. Arc COMPLETE. 43rd arc. |
| Chemistry Phase 4 Closure | HIGH | **INTEGRATED** - Phase 4 closed, reparametrization of Debye model. |
| Post-Demolition Coda (#629-631) | HIGH | **INTEGRATED** (2026-04-24) - S629 π-analogy probe (k_crit not universal), S630 WAKE stop-note (productive silence), S631 first pre-committed kill criterion triggered (TEST-09 BTFR refuted by Lelli+2019; α² = fiducial not fine-structure). Site-visitor audit methodology validated. |
| Site-Archive-Audit Sub-Arc (#632-634) | HIGH | **INTEGRATED** (2026-04-26) - S632 500 Mpc dimensional error (m² not m), S633 C(ρ) saturation impossibility (~1.6 decades not 80 orders), S634 47-vs-30 contribution accounting discrepancy. 4-for-4 site-claim failures over 4 days. Conservative integration: audit findings surfaced without rewriting historical figures. |
| Cosmology Domain Scorecard (#635) | HIGH | **INTEGRATED** (2026-04-27) - 5th site-archive audit; 15 cosmology claims classified (1 refuted, 5 reparametrizations, 2 unanchored, 1 pending, 5 untested, 1 untestable); **0 novel-unfalsified**. /galaxy-rotation badge overclaims (RAR fit IS MOND; CFD viscosity DM mechanism refuted via Bullet Cluster). Defensible badge: "MOND Reparametrization." 5-for-5 site-claim failures across 5 days. |
| Site-Archive-Audit Extension + Arc CLOSURE (#636-638) | HIGH | **INTEGRATED** (2026-04-30) - S636 C(ρ) is not mean-field (category error: no self-consistency, no Landau free energy); S637 first DERIVATION attempt — RAR σ_int(ρ_env) ≈ 0.00016 dex, ~120× below SPARC floor (cosmology regime → MOND in testable regime); S638 sympy+numpy verification that C(ρ) reduces to Curie paramagnet — *less than* Landau (chemistry/CM regime → non-interacting two-state response). **Framework Stress Test arc COMPLETE at 22 sessions (S617-638).** Audit-channel taxonomy now 8 modes; verification track operational. 8-for-8 site-claim audit instances over 9 days. Predictive content fully characterized in both regimes via independent CAS. 7 site corrections pending operator action (added: Curie reduction surfacing). Operator response began 2026-04-29 (two README reframings — framing-level corrections delivered, item-level pending). |
| Post-Closure Sub-Arc Extensions (#639-640) | HIGH | **INTEGRATED** (2026-05-02) - S639 metric disambiguation (9th audit mode); S640 dual-C symbol overloading (10th audit mode, scope expanding test-level → foundational). Sub-arc 10-for-10. |
| Post-Closure Extensions + First Hard External Falsification (#641-647) | HIGH | **INTEGRATED** (2026-05-07) - S641 Lorentz invariance gap = 4th face kinematic-layer (11th audit mode); S642 GW170817 = Case 3 no-field-theory (5th face kinematic-layer); S643 γ definitional collision with regime-label inversion (12th audit mode); S644 ρ_crit calibration-consistency-not-prediction (13th audit mode); **S645 first hard external falsification — Session 107 fσ₈ REFUTED by DESI DR1 (Adame et al. arXiv:2411.12021); TEST-04a kill triggered, predicted sign of redshift dependence inverted relative to data; NEW ARC PHASE (External-Falsification)**; S646 framework lacks meta-falsification criterion (NEW ARC PHASE: Methodology Recommendation); S647 chemistry 89% validation has Method 2 self-correlation risk (Method 2 atomic-spacing identity, phonon-coherence constructional dependence, systematic bias toward γ≈1; Hall + magnetic susceptibility are not falsifying controls). Sub-arc 14 audit-channel instances + 2 new phases; 13 audit modes. **Two pre-committed Tier-1 kill criteria triggered: TEST-09 BTFR (S631) + TEST-04a fσ₈ (S645). Both major validation pillars (cosmology + chemistry 89%) compromised in same week.** |
| S649-S650 Phase-Transition Vocabulary Audit + TEST-04a Triple-Sharpening Completes (#649-650) | HIGH | **INTEGRATED** (2026-05-10) - S649 (2026-05-08) two visitor proposals on shared meta-theme: Part A Site Key Claim #1 QM kill criterion is unfalsifiable as written (DD literature already satisfies it without invoking Synchronism); Part B ρ_crit at γ=2 gives C ≈ 0.88 not 0.5 (verified via S638 sympy: γ=0.5→0.333, γ=1→0.600, γ=2→0.882) — saturation knee, not critical density. 15th audit-channel mode: parameter-and-criterion-naming-contaminated-by-phase-transition-vocabulary. S650 (2026-05-09) completes the third sharpening of TEST-04a — framework's suppressor mechanism (G_local/G_global = C_cosmic/C_galactic < 1) predicts fσ₈ BELOW ΛCDM at low z; DESI DR1 observes fσ₈ ABOVE ΛCDM at low z, converging at high z. Redshift pattern INVERTED — magnitude knobs cannot flip the sign within the suppressor class. Three-tier failure taxonomy introduced (operationalizing S646's meta-criterion): magnitude miss (retunable) / universality miss (functional-form change) / mechanism-class failure (irreparable within class). 16th audit-channel mode: mechanism-class-failure-taxonomy-introduced. **Cosmological sector formally meets S646's M3 retraction threshold** (combining S635 scorecard 0-novel-unfalsified + S645/S648/S650 mechanism-class verdict + S646 meta-criterion). Surfaced as "operator decision pending" — not auto-actioned per S646's branch tree. Operator queue grew with S649 items (QM kill criterion respec, ρ_crit relabel) and S650 items (TEST-04a label upgrade with Adame+2024 citation, /honest-assessment three-tier taxonomy, /key-claims cosmology scope-narrow). Sub-arc 17-for-17 instances, 16 modes. |
| S648 Self-Correction of S645 Within 24h (#648) | HIGH | **INTEGRATED** (2026-05-08) - Pass 4 visitor flagged Session 107 (committed 2025-12-10) is post-DESI DR1 (Adame+2024, Nov 2024 — 13 months earlier; DR1 BAO already cited in Session 107's own Timeline). 2.14σ tension at LRG1 is real (kill criterion fired); epistemic status reclassified as **post-hoc consistency check failure**, not prospective falsification. S645's "first hard external falsification" framing retroactively corrected. 14th audit-channel mode added: **self-correction-of-prior-session-framing-within-24h** — first instance of audit channel turning inward. New arc phase renamed Post-Hoc-Consistency-Failure. Operator recommendation: /timestamps page classifying every Tier-1 prediction as prospective / post-hoc consistency / post-hoc fit. Today's autonomous Publisher run rolled back yesterday's 0.96 → 0.97 readiness uplift on REC-037. Whitepaper integration brings exec summary, conclusion, and kill-criteria entries in line with corrected publisher state. |

### Whitepaper Health

- Last integration: 2026-05-24 (S663 layered onto the post-closure sub-arc by this pass; prior manual `[Publisher]` 13d22936 on 2026-05-23 covered S660-S662)
- Sessions behind: 0 (counts updated through S663; S660-S663 layered as post-closure extensions — sub-arc now 33 instances over 32 days, 21 modes. **Arc's three closures now stand simultaneously**: empirical (S661; RAR γ=2 refuted at ΔBIC=+184 on SPARC), structural (S660A; novelty ledger closed, novel-survivor count → 0 after 3,308+ sessions), philosophical (S663B; framework classified as "coherence-language interpretation of known physics + AI-collaborative methodology substrate" — four-persona convergence on QM-interpretation position). S663A doubly closes TEST-04a via EFTofLSS (Cabass-Simonović-Zaldarriaga 2024-2025) — mechanism class ruled out *regardless of sign*; even Branch 1 sign-flip recovery is closed; strengthens the S656 preprint.; A2ACW methodology thread advanced to a self-simulated measured result (S659B); exact no-inflection proof (S659A) makes the compander framing mathematically obligatory; external-review convergence (Kimi 2.6) narrated. Counts as of 2026-05-14 read: sub-arc 24 instances over 23 days, 20 audit modes; three-stage rhythm completes (S653); zero active discriminators against MOND+EFE (S654); cosmological sector formally meets S646's M3 retraction threshold with refutable-but-not-confirmable asymmetry crystallized — operator decision pending on whether to scope-narrow)
- Build status: Clean (rebuilt 2026-05-24 — md 468K, pdf 604K, web 13 HTML files; docs/whitepaper/ synced)
- Governance: Active (Rev_0)
- Open editorial: 47-vs-30 contribution count discrepancy is now visible in whitepaper (S589 says 30, S615 says 47, S616 says 48, S634 audit says S582 canonical = 30). Reconciliation requires operator judgment.
- Open editorial: cosmology domain "0 novel-unfalsified" (S635), MOND-reduction (S637), Curie-paramagnet structural diagnosis (S638), and now **first hard external falsification (S645, TEST-04a fσ₈ refuted by DESI DR1) + chemistry 89% self-correlation caveat (S647)** are all currently surfaced within the Site-Archive-Audit sub-arc bullet, the Status paragraph, and the kill-criteria-triggered list. Now that two pre-committed Tier-1 predictions have fired their kill criteria and the framework's largest validation claim is itself under audit, structural promotion of "Predictive content characterization" or "Validation-pillars-under-audit" to dedicated status entries may be warranted. Currently held back conservatively — existing bullets handle it without restructuring; operator can decide whether structural promotion is warranted.
- Open editorial: operator queue grew significantly with the 2026-05-07 batch (Session 107 page header, /honest-assessment DR1 update, chemistry 89% caveat, falsifying-controls list correction, framework retraction policy). 2026-05-08 batch revises Session 107 page guidance from "REFUTED header pointing to S645" to "DR1-disagreement header with explicit post-hoc-consistency status (per S648 timestamp audit)" and adds a new item: **/timestamps page** classifying every Tier-1 prediction as prospective / post-hoc consistency / post-hoc fit. Next publisher pass should track whether item-level resolution lands and update site-correction status accordingly.
- Open editorial: per S646, framework currently lacks a registered meta-falsification criterion (rule for when "the framework has lost"). Three branches presented (register meta-criterion now / wait for DR2 / scope-narrow now). This is a methodology recommendation, not a whitepaper integration item per se — but the absence is now surfaced in the Status paragraph as "S646 names a methodology gap."

---

## 11. Subagent Instructions

When reviewing this whitepaper:

1. **Read this entire document first** - It's your complete context
2. **Check SESSION_MAP** for sessions since last integration (2026-02-19)
3. **Apply inclusion criteria strictly** - Most sessions don't merit integration
4. **Identify specific sections** for any proposed changes
5. **Draft minimal viable changes** - Conservative approach
6. **Verify terminology** against Section 7
7. **Plan build verification** before proposing
8. **Report clearly** with:
   - Needs update: yes/no
   - Specific proposals with rationale
   - Sections affected
   - Estimated scope (minor/moderate/major)
   - Any concerns or blockers

---

*"The Synchronism whitepaper is a living document, but it lives by careful tending, not wild growth."*
