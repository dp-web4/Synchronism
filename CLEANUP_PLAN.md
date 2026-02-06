# Synchronism Repository Cleanup Plan

**Status**: Approved with modifications
**Decisions**: See [User Decisions](#user-decisions) section

---

## User Decisions

1. **Chemistry track** - Keep as subdirectory of Research/ (it's a phenomenon in itself)
2. **Simulations/** - Add to .gitignore, track scripts only
3. **STATUS.md** - Yes, create honest assessment file
4. **Session organization** - Organize by arc, ensure future sessions follow structure
5. **memory/** - Keep in Synchronism (repo-specific), separate from private-context epistemic memory

---

## Current State Analysis

### Root Files (11 total)
| File | Keep/Move | Destination |
|------|-----------|-------------|
| README.md | Keep | Root |
| CLAUDE.md | Keep | Root |
| COHERENCE_AS_COMPRESSION.md | Move | docs/theory/ |
| COMPRESSION_ACTION_THRESHOLD.md | Move | docs/theory/ |
| PARAMETER_DEFINITIONS_AND_DERIVATIONS.md | Move | docs/theory/ |
| RESEARCH_PHILOSOPHY.md | Move | docs/why/ |
| RESEARCH_DIRECTION_2025_11.md | Move | archive/planning/ |
| RESEARCH_PIVOT_QUANTUM_COMPUTING.md | Move | archive/planning/ |
| MIGRATION_PLAN.md | Move | archive/planning/ |
| REORGANIZATION_SUMMARY_FOR_REVIEW.md | Move | archive/planning/ |
| SESSION_29_CONTEXT.md | Move | archive/legacy/ |

**Target: 2 root files** (README.md, CLAUDE.md)

### Directory Analysis

| Directory | Size | Status | Action |
|-----------|------|--------|--------|
| simulations/ | 725M | Active (session code) | Keep |
| data/ | 19M | Active | Keep |
| Research/ | 8.6M | Core content (378 sessions) | Keep, add discoveries/ |
| Documentation/ | 8.1M | Historical PDFs | Consolidate with docs/ |
| manuscripts/ | 3.9M | Publications | Keep |
| whitepaper/ | 3.2M | Publication | Keep |
| scripts/ | 2.6M | Active code | Keep |
| figures/ | 2.4M | Visualizations | Keep |
| forum/ | 1.9M | A2ACW sessions | Keep |
| docs/ | ~1M | Technical docs | Expand as main docs/ |
| explorations/ | ~500K | Oct 2025 experiments | Archive |
| Mathematical_Frameworks/ | ~500K | Mixed code/PDFs | Consolidate |
| Experimental/ | ~100K | MRH experiments | Archive or merge |
| governance/ | ~100K | Governance system | Keep |
| memory/ | ~100K | Knowledge DB | Keep |
| articles/ | ~50K | News items | Keep |
| supplementary/ | ~50K | Validation code | Merge to scripts/ |
| web-version/ | 1.6M | Website | Keep |

### Overlapping Directories
- `Documentation/` and `docs/` both exist
- `Mathematical_Frameworks/` has code that belongs in scripts/
- `supplementary/` has validation code that belongs in scripts/
- `Experimental/` and `explorations/` both hold old experiments

---

## Proposed Structure

```
Synchronism/
├── README.md                    # Navigation hub (rewrite)
├── CLAUDE.md                    # AI session context
├── Research/                    # 378+ sessions (KEEP AS-IS)
│   ├── SESSION_MAP.md
│   ├── Chemistry/
│   ├── Gnosis/
│   ├── Open_Questions/
│   └── discoveries/            # NEW: Major validated findings
├── docs/                        # EXPANDED
│   ├── README.md               # Documentation navigation
│   ├── theory/                 # Core theory documents
│   ├── why/                    # Philosophy and motivation
│   ├── how/                    # Technical guides
│   └── reference/              # PDFs and formal documents
├── archive/                     # NEW
│   ├── README.md               # Archive navigation
│   ├── planning/               # Old planning docs
│   ├── legacy/                 # Obsolete files
│   └── explorations/           # Oct 2025 experiments
├── simulations/                 # Session simulation code
├── scripts/                     # Analysis and utility code
├── manuscripts/                 # Publications
├── whitepaper/                  # Whitepaper versions
├── forum/                       # A2ACW adversarial sessions
├── governance/                  # LRC governance system
├── data/                        # Datasets
├── figures/                     # Visualizations
├── articles/                    # News/articles
├── memory/                      # Knowledge DB
└── web-version/                 # Website
```

---

## Major Discoveries to Surface

Like HRM's `docs/what/discoveries/`, create `Research/discoveries/`:

### 1. γ ~ 1 Universal Boundary
- **Status**: Validated across 1703 phenomenon types
- **Sessions**: Chemistry #1-1840
- **Evidence**: 89% prediction validation rate
- **Key insight**: γ = 2/√N_corr unifies quantum-classical transition

### 2. NP2 RAR Scatter Validation
- **Status**: Strongly supported (p = 5×10⁻⁶)
- **Sessions**: #376-378 (Gas Fraction Control Arc)
- **Evidence**: Bootstrap, permutation, non-parametric all confirm
- **Key insight**: Morphology-dependent scatter supports environment effect

### 3. Gnosis Consciousness Threshold
- **Status**: Theory complete
- **Sessions**: Gnosis #1-11
- **Evidence**: 8-way convergence at C ≈ 0.50
- **Key insight**: Self-awareness emerges at coherence threshold

### 4. η Formalism for Superconductivity
- **Status**: Validated (cuprate predictions match)
- **Sessions**: #292, #297
- **Evidence**: YBCO η = 0.38 extracted from data
- **Key insight**: Reachability factor determines T_c potential

### 5. ONE EQUATION Unification
- **Status**: Theoretical + partial validation
- **Sessions**: Integration Arc #360-363
- **Evidence**: Spans 80 orders of magnitude
- **Key insight**: γ = 2/√N_corr connects all scales

---

## Implementation Phases

### Phase 1: Create New Structure
1. Create `docs/theory/`, `docs/why/`, `docs/how/`, `docs/reference/`
2. Create `archive/`, `archive/planning/`, `archive/legacy/`, `archive/explorations/`
3. Create `Research/discoveries/`
4. Create `Research/arcs/` with arc index files

### Phase 2: Move Root Files (9 files)
Using `git mv` to preserve history:

| File | Destination |
|------|-------------|
| COHERENCE_AS_COMPRESSION.md | docs/theory/ |
| COMPRESSION_ACTION_THRESHOLD.md | docs/theory/ |
| PARAMETER_DEFINITIONS_AND_DERIVATIONS.md | docs/theory/ |
| RESEARCH_PHILOSOPHY.md | docs/why/ |
| RESEARCH_DIRECTION_2025_11.md | archive/planning/ |
| RESEARCH_PIVOT_QUANTUM_COMPUTING.md | archive/planning/ |
| MIGRATION_PLAN.md | archive/planning/ |
| REORGANIZATION_SUMMARY_FOR_REVIEW.md | archive/planning/ |
| SESSION_29_CONTEXT.md | archive/legacy/ |

### Phase 3: Consolidate Directories
1. Move `Documentation/` PDFs → `docs/reference/`
2. Move `Documentation/governance/` → `docs/governance/` (or keep separate)
3. Move `explorations/` → `archive/explorations/`
4. Move `Mathematical_Frameworks/*.py` → `scripts/math/`
5. Move `Mathematical_Frameworks/*.pdf` → `docs/reference/`
6. Move `supplementary/*.py` → `scripts/validation/`
7. Move `Experimental/*.py` → `scripts/experimental/`

### Phase 4: Create STATUS.md
Write honest assessment file:
- What works (validated predictions, γ ~ 1 universality)
- What doesn't (melting point 53% error, critical exponents 2×)
- What's untested (34 Gnosis predictions, OQ005/OQ006)
- How to evaluate this work

### Phase 5: Update .gitignore
Add simulation outputs:
```
simulations/*.csv
simulations/*.json
simulations/*.pkl
simulations/*.npy
simulations/output/
simulations/results/
```

### Phase 6: Write Discovery Documentation
Create summaries in `Research/discoveries/`:
- `gamma-universal-boundary.md` - 1703 phenomenon types, 89% validated
- `np2-rar-scatter-validation.md` - p = 5×10⁻⁶, Gas Fraction Arc
- `gnosis-consciousness-threshold.md` - C ≈ 0.50, 8-way convergence
- `eta-superconductivity-formalism.md` - Cuprate predictions match
- `one-equation-unification.md` - 80 orders of magnitude

### Phase 7: Create Arc Index Files
In `Research/arcs/`:
- `README.md` - Arc navigation guide
- One `.md` file per major arc with session links
- Template for future arc documentation

### Phase 8: Update README
Rewrite README.md as navigation hub:
- Lead with achievements (like HRM)
- Current status snapshot
- Audience-based navigation (Newcomer/Researcher/Developer/AI)
- Link to SESSION_MAP.md and discoveries/

### Phase 9: Update CLAUDE.md
- Add arc structure guidance for future sessions
- Update any changed paths
- Add STATUS.md reference

### Phase 10: Verification
- Check all internal links
- Test navigation from fresh perspective
- Verify autonomous session scripts still work

---

## Arc Organization Strategy

SESSION_MAP.md already documents 14+ arcs clearly. Options for physical organization:

### Option A: Full Reorganization (Disruptive)
Move all 377 session files into arc directories:
```
Research/arcs/
├── active/
│   ├── experimental-validation/   # 368-378+
│   └── hot-superconductor/        # 292, 297-300
├── complete/
│   ├── technology/                # 364-367
│   ├── integration/               # 360-363
│   ├── consciousness-2.0/         # 356-359
│   └── ... (14+ arc directories)
└── early/                         # Pre-arc sessions (1-~250)
```
**Pros**: Clean organization, easy navigation
**Cons**: Breaks all existing links, disrupts autonomous sessions

### Option B: Hybrid (Recommended)
1. Keep existing session files in place
2. Create `Research/arcs/` with index files pointing to sessions
3. New sessions go directly into arc directories
4. Eventually migrate as sessions are touched

```
Research/
├── arcs/
│   ├── README.md                  # Arc navigation
│   ├── experimental-validation.md # Index for #368-378
│   ├── hot-superconductor.md      # Index for #292, 297-300
│   └── ... (index files)
├── Session368_*.md                # Existing files stay
├── Session369_*.md
└── ...
```

### Option C: Future-Only
1. Keep all existing files as-is
2. Create arc directories only for NEW sessions
3. SESSION_MAP.md remains the primary navigation

**Selected: Option B** - Preserves existing structure while enabling arc-based navigation for new work.

---

## Simulations .gitignore

Add to `.gitignore`:
```
# Simulation outputs (scripts tracked, data regenerable)
simulations/*.csv
simulations/*.json
simulations/*.pkl
simulations/*.npy
simulations/output/
simulations/results/
```

Keep tracked:
- `simulations/*.py` (session scripts)
- `simulations/README.md`
- `simulations/requirements.txt` (if exists)

---

## Comparison to HRM Cleanup

| Metric | HRM Before | HRM After | Synchronism Before | Synchronism Target |
|--------|------------|-----------|--------------------|--------------------|
| Root files | 73 | 7 | 11 | 2-3 |
| Discovery docs | 0 | 5 | 0 | 5 |
| Overlapping dirs | 3+ | 0 | 4 | 0 |
| Archive structure | None | Organized | Scattered | Organized |
| Newcomer orientation | Poor | Good | Medium | Good |

---

## Success Criteria

After implementation:
- [ ] Root has ≤3 essential files (README, CLAUDE.md, STATUS.md)
- [ ] README highlights achievements on first screen
- [ ] STATUS.md provides honest assessment
- [ ] All 378+ sessions navigable via SESSION_MAP
- [ ] Major discoveries have dedicated documentation (5 files)
- [ ] Historical/planning docs archived (not deleted)
- [ ] No overlapping directories (Documentation/, docs/ consolidated)
- [ ] Arc structure in place for future sessions
- [ ] Simulation outputs gitignored
- [ ] All internal links work
- [ ] Autonomous session scripts unaffected

---

## Estimated Effort

| Phase | Tasks | Priority |
|-------|-------|----------|
| 1. Create structure | mkdir commands | HIGH |
| 2. Move root files | 9 git mv | HIGH |
| 3. Consolidate dirs | ~20 git mv | MEDIUM |
| 4. STATUS.md | Write new file | HIGH |
| 5. .gitignore | Update file | HIGH |
| 6. Discovery docs | 5 new files | HIGH |
| 7. Arc indexes | ~15 index files | MEDIUM |
| 8. README rewrite | 1 file | HIGH |
| 9. CLAUDE.md update | 1 file | MEDIUM |
| 10. Verification | Link checks | HIGH |

Can be done in one focused session or spread across multiple.

---

*Plan created: February 5, 2026*
*Approved: February 5, 2026*
*Based on HRM reorganization pattern*
