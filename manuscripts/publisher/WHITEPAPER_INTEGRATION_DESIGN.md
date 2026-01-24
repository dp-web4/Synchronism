# Publisher Whitepaper Integration Design

**Version**: 2.0
**Date**: 2026-01-23
**Status**: Design Document

---

## Overview

The Publisher track is expanding from "catalog and recommend" to include **active whitepaper maintenance**. This document outlines the architecture for managing two whitepapers:

1. **Synchronism Whitepaper** - Core coherence physics framework
2. **Web4 Whitepaper** - Trust-native distributed intelligence protocol

---

## Design Principles

### 1. Meaning Over Mechanics

Whitepapers are the primary interface between research and the world. Changes must:
- Preserve conceptual integrity
- Maintain narrative coherence
- Respect established terminology
- Enhance rather than dilute

### 2. Comprehensive Review, Conservative Changes

Each review cycle should:
- Thoroughly evaluate new research against whitepaper scope
- Identify genuine improvements (not just additions)
- Make changes only when truly merited
- Document rationale for all decisions

### 3. Subagent Isolation

Each whitepaper gets its own subagent with:
- Full context on structure, build process, history
- Clear inclusion criteria
- Authority to propose changes
- No cross-contamination between whitepapers

### 4. Respect Existing Governance

Synchronism whitepaper has an existing governance structure (proposals/reviews/arbiter).
Web4 whitepaper has a simpler model.
Publisher adapts to each, doesn't override.

---

## Architecture

### Publisher Orchestration Flow

```
Publisher Daily Run (02:30 UTC)
│
├── Phase 1: Catalog & Recommend (existing)
│   ├── Scan SESSION_MAP for new sessions
│   ├── Update recommendations.json
│   └── Generate daily report
│
├── Phase 2: Whitepaper Review (NEW)
│   ├── Launch Synchronism Whitepaper Subagent
│   │   └── Returns: {needs_update: bool, proposals: [], report: str}
│   │
│   ├── Launch Web4 Whitepaper Subagent
│   │   └── Returns: {needs_update: bool, proposals: [], report: str}
│   │
│   └── Merge subagent reports into daily report
│
└── Phase 3: Commit & Push
    ├── Commit all changes
    └── Push to remotes
```

### Subagent Context Structure

Each whitepaper directory gets a `PUBLISHER_CONTEXT.md` providing:

```
PUBLISHER_CONTEXT.md
├── Whitepaper Purpose & Philosophy
├── Section Structure (what each section covers)
├── Inclusion Criteria (what merits whitepaper inclusion)
├── Exclusion Criteria (what stays in research, not whitepaper)
├── Build Process (how to regenerate outputs)
├── Recent Changes (last 5 integrations with dates)
├── Terminology Protection (canonical terms, forbidden redefinitions)
├── Quality Standards (style, rigor, formatting)
└── Integration Workflow (how to propose and implement changes)
```

---

## Inclusion Criteria Framework

### For Synchronism Whitepaper

**INCLUDE when research:**
- Derives new results from coherence framework
- Validates existing predictions quantitatively
- Provides new cross-domain connections
- Fills documented gaps in the framework
- Has reached natural terminus (complete arc)

**EXCLUDE when research:**
- Is still actively evolving
- Contradicts established framework without resolution
- Adds complexity without proportional insight
- Is domain-specific without universal implications
- Would require restructuring multiple sections

### For Web4 Whitepaper

**INCLUDE when research:**
- Clarifies or extends LCT/T3/R6 specifications
- Provides implementation patterns for adoption
- Addresses documented gaps in protocol
- Has been validated in hardbound/HRM implementations

**EXCLUDE when research:**
- Is theoretical without implementation path
- Belongs in Synchronism (physics) not Web4 (protocol)
- Would break backward compatibility without justification
- Is speculative future direction

---

## Subagent Responsibilities

### Synchronism Whitepaper Subagent

**Inputs:**
- Full PUBLISHER_CONTEXT.md
- Recent SESSION_MAP entries
- Publisher recommendations (relevant ones)
- Current section CHANGELOGs

**Process:**
1. Load whitepaper context
2. Scan new sessions since last review
3. Evaluate each against inclusion criteria
4. For included sessions:
   - Identify target section(s)
   - Draft integration approach
   - Check for terminology consistency
   - Estimate scope of changes
5. Generate review report

**Outputs:**
- Needs update: yes/no
- Proposed integrations: list with rationale
- Recommended section changes
- Concerns or blockers
- Summary for Publisher daily report

### Web4 Whitepaper Subagent

**Inputs:**
- Full PUBLISHER_CONTEXT.md
- Recent web4 developments (from hardbound, HRM)
- Current ARCHITECTURE.md files
- Protocol specifications

**Process:**
1. Load whitepaper context
2. Check for new implementations (hardbound-core, web4-core)
3. Check for specification updates
4. Evaluate documentation currency
5. Identify gaps or inconsistencies
6. Generate review report

**Outputs:**
- Needs update: yes/no
- Proposed sections to update
- New content to add
- Deprecated content to remove
- Summary for Publisher daily report

---

## Change Workflow

### Conservative Change Process

```
1. Subagent identifies potential update
   │
2. Evaluates against inclusion criteria
   │ (Most stop here - criteria not met)
   │
3. Drafts specific changes with rationale
   │
4. Checks for unintended consequences
   │ - Terminology drift
   │ - Narrative disruption
   │ - Cross-reference breakage
   │
5. If all checks pass → Propose change
   │
6. For Synchronism: Create proposal in meta/proposals/
   For Web4: Direct edit (simpler governance)
   │
7. Log in CHANGELOG.md
   │
8. Rebuild outputs (make-*.sh)
   │
9. Commit with detailed message
```

### Change Types

| Type | Description | Threshold |
|------|-------------|-----------|
| **Minor** | Typos, formatting, clarifications | Low - direct edit |
| **Moderate** | New examples, updated statistics | Medium - log rationale |
| **Major** | New sections, restructuring | High - full review cycle |
| **Critical** | Core framework changes | Very High - human review required |

---

## File Structure

### Synchronism Whitepaper Additions

```
whitepaper/
├── PUBLISHER_CONTEXT.md    # NEW - Subagent context
├── sections/
│   └── (existing structure)
├── build/
└── (existing files)
```

### Web4 Whitepaper Additions

```
whitepaper/
├── PUBLISHER_CONTEXT.md    # NEW - Subagent context
├── sections/
│   └── (existing structure)
├── build/
└── (existing files)
```

### Publisher Directory Updates

```
publisher/
├── CLAUDE.md               # Updated with Phase 2
├── state/
│   ├── recommendations.json
│   ├── published.json
│   ├── whitepaper_sync.json    # NEW - Last sync state
│   └── whitepaper_web4.json    # NEW - Last sync state
├── reports/
│   └── (includes whitepaper sections)
└── (existing files)
```

---

## Safety Constraints

1. **Never auto-commit critical changes** - Flag for human review
2. **Always preserve existing content** - Archive before major changes
3. **Respect terminology protection** - Canonical terms are immutable
4. **Build must succeed** - No commits if build fails
5. **Log everything** - Full audit trail of decisions

---

## Implementation Phases

### Phase 1: Context Documents (Today)
- Create PUBLISHER_CONTEXT.md for Synchronism whitepaper
- Create PUBLISHER_CONTEXT.md for Web4 whitepaper
- Document inclusion criteria, structure, process

### Phase 2: Publisher CLAUDE.md Update (Today)
- Add Phase 2 workflow to Publisher
- Define subagent launch protocol
- Update daily report format

### Phase 3: State Tracking (Today)
- Add whitepaper_sync.json and whitepaper_web4.json
- Track last review date, sessions reviewed, pending proposals

### Phase 4: Testing (Subsequent Session)
- Run full cycle manually
- Verify subagent isolation
- Validate change workflow
- Confirm build integration

---

## Success Criteria

1. **Whitepapers stay current** - Research integrated within 7 days of arc completion
2. **No meaning drift** - Terminology and framework preserved
3. **Quality maintained** - No degradation of narrative coherence
4. **Audit trail complete** - All changes documented and attributable
5. **Build stability** - Zero broken builds from Publisher changes

---

*"The whitepaper is the face of the research. The Publisher ensures that face reflects truth, not just activity."*
