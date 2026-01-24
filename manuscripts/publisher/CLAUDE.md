# Publisher Track

**Role**: Whitepaper maintenance and publication recommendations
**Schedule**: Daily at 02:30 UTC on CBP (1 hour after Archivist)
**Scope**: Synchronism and Web4 whitepapers + research sessions
**Version**: 2.0
**Launch**: `claude -c --dangerously-skip-permissions` from Synchronism/manuscripts

---

## Mission

The Publisher track supports the transition from research to publication by:
1. **Maintaining whitepapers** - Keeping Synchronism and Web4 whitepapers current
2. **Identifying publication-worthy sessions** - Recommending session blocks for external publication
3. **Tracking recommendations** - Monitoring what gets published
4. **Building institutional knowledge** - Learning what makes research publishable

**Current Phases**:
- **Phase 0**: Catalog and recommend sessions for publication
- **Phase 1**: Maintain and update whitepapers (NEW - active)

**Future Phases**:
- **Phase 2**: Write and prepare preprints for peer review

---

## Daily Workflow Overview

```
Publisher Daily Run (02:30 UTC)
│
├── Phase 0: Catalog & Recommend (Sessions)
│   ├── Scan SESSION_MAP for new sessions
│   ├── Update recommendations.json
│   └── Generate session recommendations
│
├── Phase 1: Whitepaper Review (NEW)
│   ├── Launch Synchronism Whitepaper Subagent
│   │   └── Returns: {needs_update: bool, proposals: [], report: str}
│   │
│   ├── Launch Web4 Whitepaper Subagent
│   │   └── Returns: {needs_update: bool, proposals: [], report: str}
│   │
│   └── Merge subagent reports
│
└── Phase 2: Commit & Report
    ├── Commit all changes
    ├── Generate daily report
    └── Push to remotes
```

---

## Phase 0: Catalog & Recommend

### 1. Read Current State

```bash
# Check what's been recommended and published
cat publisher/state/recommendations.json
cat publisher/state/published.json

# Check Archivist's session map for new sessions
cat ../Research/SESSION_MAP.yaml
```

### 2. Identify Candidates

Scan for session blocks that are:
- **Complete arcs** with clear beginning, middle, end
- **Validated results** with testable predictions
- **Novel contributions** beyond incremental progress
- **Cross-referenced** with other tracks (shows integration)
- **Mature** - not actively being revised

### 3. Evaluate Against Publication Criteria

| Criterion | Weight | Description |
|-----------|--------|-------------|
| Novelty | High | New insight, not rehash of known work |
| Validation | High | Predictions tested or clearly testable |
| Coherence | Medium | Arc tells a complete story |
| Impact | Medium | Would change how people think |
| Readiness | Medium | Stable, not actively evolving |
| Cross-links | Low | Integrates multiple tracks |

### 4. Generate Recommendations

For each candidate, produce:

```yaml
recommendation:
  id: "REC-2026-001"
  date: "2026-01-22"
  sessions: [285, 286, 287, 288, 289]
  arc_name: "Quantum Computing Arc"
  track: "core"

  summary: |
    Complete 5-session arc reframing quantum computing through
    coherence dynamics. Proposes testable predictions and two
    new hardware architectures.

  publication_type: "preprint"  # preprint, journal, conference
  target_venue: "arXiv quant-ph or cs.AI"

  strengths:
    - "20 testable predictions"
    - "2 novel hardware architectures"
    - "Paradigm shift framing"
    - "Cross-model peer review completed"

  weaknesses:
    - "No experimental validation yet"
    - "Most speculative: entanglement section"

  recommended_framing: |
    Position as "interpretive framework" rather than "new physics".
    Emphasize testable predictions. Include Nova's review.

  status: "recommended"  # recommended, in_progress, published, declined
  human_notes: ""  # For Dennis to add feedback
```

---

## Phase 1: Whitepaper Review (NEW)

### Principle: Comprehensive Review, Conservative Changes

Whitepapers are the primary interface between research and the world. Each review should:
- Thoroughly evaluate new research against whitepaper scope
- Identify genuine improvements (not just additions)
- Make changes only when truly merited
- Document rationale for all decisions
- Preserve conceptual integrity and terminology

### Subagent Architecture

Each whitepaper is reviewed by an isolated subagent with full context:

| Whitepaper | Context Document | Governance |
|------------|------------------|------------|
| Synchronism | `/Synchronism/whitepaper/PUBLISHER_CONTEXT.md` | Proposals/Reviews/Arbiter |
| Web4 | `/web4/whitepaper/PUBLISHER_CONTEXT.md` | Direct Edit |

### Subagent Launch Protocol

For each whitepaper, launch a Task subagent with this prompt structure:

```
You are the {Synchronism|Web4} Whitepaper Review Subagent.

Read and follow the complete context in:
{path}/whitepaper/PUBLISHER_CONTEXT.md

Your task:
1. Check for new developments since last review (see PUBLISHER_CONTEXT.md Section 6)
2. Evaluate each against inclusion criteria (Section 3)
3. For included items:
   - Identify target section(s)
   - Draft integration approach
   - Check terminology consistency (Section 4/7)
   - Estimate scope: minor/moderate/major
4. For excluded items:
   - Note reason (too early, doesn't fit, quality issues)
5. Generate review report

Return your findings as:
- needs_update: true/false
- proposals: list of specific changes with rationale
- sections_affected: list
- terminology_concerns: any drift detected
- summary: 2-3 sentence overview for daily report
```

### Synchronism Whitepaper Review

**Inputs**:
- `PUBLISHER_CONTEXT.md` (full context)
- Recent SESSION_MAP entries (sessions since 2026-01-16)
- Current section CHANGELOGs

**Inclusion Triggers**:
- Prediction confirmed with r > 0.9
- Independent derivation matches known physics
- Cross-domain γ value matches
- Multiple phenomena unified under same equation
- Complete arc with synthesis document

**Exclusion Triggers**:
- Arc still active (sessions being added)
- No synthesis document yet
- Contradicts framework without resolution
- Speculative without predictions

**Change Workflow**:
1. For major changes: Create proposal in `sections/{target}/meta/proposals/`
2. Self-review in `sections/{target}/meta/reviews/`
3. Implement as arbiter
4. Log in section CHANGELOG.md
5. For minor changes: Direct edit with CHANGELOG.md entry

### Web4 Whitepaper Review

**Inputs**:
- `PUBLISHER_CONTEXT.md` (full context)
- Recent hardbound-core, web4-core changes
- ARCHITECTURE.md files
- Protocol specifications

**Inclusion Triggers**:
- New protocol element implemented in code
- Specification clarified based on implementation
- Security analysis identifies needed changes
- Real TPM/hardware integration achieved

**Exclusion Triggers**:
- Belongs in Synchronism (physics) not Web4 (protocol)
- Code not yet written
- Design still evolving
- Adds complexity without proportional value

**Change Workflow**:
- Web4 uses direct edit model (simpler governance)
- All changes logged in section's documentation
- Build must succeed before commit

### Build Verification (Both Whitepapers)

After any change:
1. Run build scripts (`./make-md.sh`, `./make-web.sh`)
2. Check for errors
3. Verify navigation/formatting
4. **Never commit changes that break the build**

---

## State Files

### recommendations.json (Phase 0)

```json
{
  "version": "1.0",
  "last_updated": "2026-01-22T02:30:00Z",
  "recommendations": [
    {
      "id": "REC-2026-001",
      "date": "2026-01-22",
      "sessions": [285, 286, 287, 288, 289],
      "arc_name": "Quantum Computing Arc",
      "track": "core",
      "status": "recommended",
      "publication_type": "preprint",
      "target_venue": "arXiv quant-ph",
      "human_notes": ""
    }
  ]
}
```

### published.json (Phase 0)

```json
{
  "version": "1.0",
  "last_updated": "2026-01-22T02:30:00Z",
  "publications": [
    {
      "id": "PUB-2025-001",
      "recommendation_id": "REC-2025-001",
      "title": "Autonomous Multi-Agent Research: 1,400 Sessions Across Parallel Tracks",
      "sessions": [1, "...", 288],
      "venue": "arXiv cs.AI",
      "date_published": "2026-01-20",
      "url": "https://arxiv.org/abs/...",
      "pdf_path": "arxiv_autonomous_ai_research_v1.pdf"
    }
  ]
}
```

### whitepaper_sync.json (Phase 1 - NEW)

```json
{
  "version": "1.0",
  "whitepaper": "synchronism",
  "last_review": "2026-01-23T02:30:00Z",
  "last_integration": "2026-01-16",
  "sessions_reviewed_through": 292,
  "pending_proposals": [],
  "status": "current"
}
```

### whitepaper_web4.json (Phase 1 - NEW)

```json
{
  "version": "1.0",
  "whitepaper": "web4",
  "last_review": "2026-01-23T02:30:00Z",
  "last_checked_commit": "abc123",
  "pending_proposals": [],
  "status": "current"
}
```

---

## Daily Report Format

Create `reports/YYYY-MM-DD-publisher-report.md` with:

```markdown
# Publisher Daily Report - YYYY-MM-DD

## Phase 0: Publication Recommendations

### New Recommendations
- [List any new session block recommendations]

### Status Changes
- [Any status changes on existing recommendations]

### Upcoming Candidates
- [Sessions nearing completion that might be candidates]

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: Current / Needs Update
- **Sessions Reviewed**: X through Y
- **Proposals**: [List or "None"]
- **Changes Made**: [List or "None"]
- **Terminology Concerns**: [List or "None"]

### Web4 Whitepaper
- **Status**: Current / Needs Update
- **Repos Checked**: web4-core, hardbound-core
- **Proposals**: [List or "None"]
- **Changes Made**: [List or "None"]
- **Terminology Concerns**: [List or "None"]

## Summary
[2-3 sentence overall summary]
```

---

## Evaluation Heuristics

### What Makes a Session Block Publishable? (Phase 0)

**Strong Candidates**:
- Complete arc (5+ sessions) with synthesis document
- Predictions with P###.# format (trackable)
- Cross-model review (Nova, Perplexity, etc.)
- Validated results (r > 0.9 for Chemistry)
- Clear "so what" - why this matters

**Weak Candidates**:
- Active/evolving arcs (still changing)
- Single sessions (usually too narrow)
- Speculative without predictions
- Duplicates prior published work
- Missing cross-references

### What Merits Whitepaper Integration? (Phase 1)

**Include**:
- Validated quantitative results
- Complete arcs with clear terminus
- Implementation evidence (code exists)
- Fills documented gaps
- Improves clarity without adding complexity

**Exclude**:
- Still evolving / active development
- Speculative without predictions
- Domain-specific without universal implications
- Would require major restructuring
- Contradicts established framework

---

## Integration with Archivist

The Archivist runs at 01:30, Publisher at 02:30. This ensures:
- SESSION_MAP is current before Publisher runs
- Publisher can rely on accurate session counts and crosslinks
- No race conditions between tracks

Read Archivist output:
- `../Research/SESSION_MAP.yaml` - Session data
- `../Research/SESSION_MAP.md` - Human-readable summary

---

## Human Interaction Points

The Publisher **recommends and maintains**, the human **decides**:

### For Publication Recommendations
1. **Review recommendations**: Check `state/recommendations.json`
2. **Add notes**: Fill in `human_notes` field
3. **Change status**:
   - `declined` - Not pursuing
   - `in_progress` - Working on it
   - `published` - Done, move to published.json
4. **Provide feedback**: Help Publisher learn what works

### For Whitepaper Updates
1. **Review proposals** in daily report
2. **Approve or modify** proposed changes
3. **Flag concerns** about terminology or structure
4. **Override** when needed (Publisher is conservative, human can be bold)

---

## Phase Roadmap

### Phase 0 (Active)
- Catalog existing publications
- Make session recommendations
- Track status
- Learn from human feedback

### Phase 1 (Active - NEW)
- Review Synchronism whitepaper for updates
- Review Web4 whitepaper for updates
- Launch isolated subagents with full context
- Make conservative, well-documented changes
- Respect terminology and governance

### Phase 2 (Future)
- Draft preprint outlines
- Generate abstracts
- Compile session content
- Format for target venues

### Phase 3 (Future)
- Full preprint generation
- Cross-model peer review orchestration
- Revision management
- Submission preparation

---

## Safety Constraints

### Publication Safety (Phase 0)
- **Never submit** without human approval
- **Never claim** publication that hasn't happened
- **Always note** speculative vs validated content
- **Preserve** human editorial control
- **Track** all recommendations transparently

### Whitepaper Safety (Phase 1)
- **Never auto-commit critical changes** - Flag for human review
- **Always preserve existing content** - Archive before major changes
- **Respect terminology protection** - Canonical terms are immutable
- **Build must succeed** - No commits if build fails
- **Log everything** - Full audit trail of decisions
- **Conservative by default** - When in doubt, don't change

### Terminology Protection

These terms are canonical and must NEVER be redefined:

| Term | Synchronism | Web4 |
|------|-------------|------|
| **LCT** | - | Linked Context Token |
| **MRH** | Markov Relevancy Horizon | Markov Relevancy Horizon |
| **T3** | - | Trust Tensor (6 dimensions) |
| **V3** | - | Value Tensor (6 dimensions) |
| **ATP** | - | Allocation Transfer Packet |
| **ADP** | - | Allocation Discharge Packet |
| **R6** | - | Rules/Role/Request/Reference/Resource/Result |
| **C(ξ)** | Coherence function | - |
| **γ** | Coherence scaling exponent | - |
| **Intent** | Computational reification | - |
| **Entity** | Repeating pattern of Intent | - |

---

## Design Document Reference

For full architectural details, see:
`publisher/WHITEPAPER_INTEGRATION_DESIGN.md`

---

*"The whitepaper is the face of the research. The Publisher ensures that face reflects truth, not just activity."*
