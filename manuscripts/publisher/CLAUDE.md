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

### 1b. Scan Surfaces (the scan has been blind twice — keep this list current)

Publication-relevant work does not only arrive as numbered sessions, and it does not only arrive
in this repo. Scan **all** of these every run:

| Surface | Why it is on this list |
|---|---|
| `../Research/SESSION_MAP.yaml` | numbered sessions (the original scan) |
| `../Research/papers/` | **added 2026-07-27** — a finished 4,307-word manuscript sat unseen for 4 days because the scan was shaped for session arcs (→ REC-2026-038) |
| `../Research/proposals/` | where the research lane files framing changes and registrations |
| `../Research/preregistrations/` | pre-registered protocols; tells you what is *about* to be executed |
| `../explorations/` | same-day triage notes; often the first place a self-correction lands |
| **`synchronism-site/explorer/findings/`** | **added 2026-07-30** — a *sibling repo*. The boost-ceiling class exclusion was executed here on 2026-06-03, and on 2026-07-29 the Publisher spent a pass re-deriving a weaker version of it and recommended that as the queue's lead physics item. Executed results with real computation live here, and nothing in the Synchronism repo indexes them. |
| `synchronism-site/explorer/topics/` | seeded-but-unrun executions — the pipeline of what is coming |

**The failure class both additions belong to**: the scan was shaped for the form the *last* finding
arrived in. Before opening or promoting a recommendation, grep the sibling repos for the claim's own
keywords — the program's own prior art is cheaper to find than the literature's and is more often
the thing that moves the verdict.

**Scan sibling repos at `origin/<default>`, never at HEAD (added 2026-08-09).** These checkouts are
shared with other agents and are routinely parked on a feature branch. On 2026-08-09 the `web4`
working tree sat on `kimi/purpose-is-relational` at **2026-08-03**, so a bare `git log --since` read
**0 commits** for the window while `origin/main` carried **~20 commits dated 2026-08-08**. The
verdict happened to survive (all 20 were hub/docs/audit; exactly one touched `whitepaper/`, and it
was the Publisher's own log) — *the instrument was wrong and the answer was right by luck*, which is
the failure mode that persists longest because nothing contradicts it. Web4's own lane fixed the
identical bug in its CI the same day (`55c0ed7`, "publish the merged ref, not whatever branch the
shared checkout is parked on"). So:

```bash
git -C <repo> log origin/main --since=<window> --name-only -- whitepaper/   # scope
git -C <repo> rev-parse --abbrev-ref HEAD                                    # and SAY where HEAD was
```

Report the branch you scanned. A commit count with no ref named is a claim about a checkout you did
not inspect. Same rule for `synchronism-site`.

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

### Artifact Commits — CI is the authoritative builder (added 2026-08-01)

**Do not commit locally rebuilt Synchronism artifacts.** `docs/whitepaper/**` and
`whitepaper/build/**` are regenerated and committed by `build_whitepaper.yml` on every push that
touches `whitepaper/**`. Rebuild locally to *verify*; then `git checkout` the tree and commit only
the sources. Stage the artifact only if CI is known dead (it is not — it deployed on 07-26, 07-27,
07-28 and 08-01).

The reason is measured, not stylistic. This working tree is a WSL `/mnt/c` CRLF mount; CI builds on
LF. The same monolith carries **492 CRLF lines as CI writes it and 4212 as a local rebuild writes
it**, so a locally built artifact diffs against the committed blob on ~5,600 lines per file with
zero content difference. On 2026-08-01 this round-tripped in public: `a7d43b92` committed **11,252
lines** of line-ending churn across the two monoliths to carry **12 lines** of real content, and
CI's next build (`962a8add`) reverted it with **11,158 insertions and 11,158 deletions and a
CR-stripped diff that is empty**. 22,410 lines of content-free history for 12 lines of content.

**The gate that missed it, and the fix.** The 08-01 pass ran the right check and drew the wrong
label from it:

> Monolith diff (`--ignore-cr-at-eol`) → 6 insertions, 0 deletions — *no churn*

`--ignore-cr-at-eol` is the flag that makes the **content** assertion correct and it is the same
flag that blinds the check to the **churn**. The predicate cannot see the property the label names.
So run both and report both numbers:

```bash
git diff --stat --ignore-cr-at-eol -- docs/whitepaper/ whitepaper/build/   # content: must be the intended lines
git diff --stat                     -- docs/whitepaper/ whitepaper/build/   # raw: must be the SAME lines
```

If the raw number exceeds the content number, you are holding line-ending churn — restore the tree
(`git checkout -- docs/whitepaper/ whitepaper/build/`) and let CI build. **Never report a single
CR-stripped number as "no churn."** More generally: when a flag is required to make a check pass,
name what that flag hides and check for it separately.

**Root fix, dp's call and still open:** a `.gitattributes` entry (`docs/whitepaper/** text eol=lf`,
same for `whitepaper/build/**`) would retire the class outright. Flagged 07-26 and 07-31 as dp's
call because it rewrites those files once; the 08-01 round trip is the first measured cost.

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

## Collective Coordination Logging (NEW - 2026-02-05)

**CRITICAL**: Every Publisher run MUST also update the collective coordination log.

### Why This Matters

The autonomous collective uses a coordinated logging pattern:
- **Archivist** (04:00 UTC) → logs to `private-context/archivist/log.md`
- **Publisher** (02:30 UTC) → logs to `private-context/publisher/log.md`
- **Supervisor** (03:30 UTC) → reads both logs, verifies health

Without Publisher logging, the Supervisor cannot verify Publisher health and the coordination loop is broken.

### At Session Start

Read the Archivist log for context:

```bash
head -30 /mnt/c/exe/projects/ai-agents/private-context/archivist/log.md
```

Use this to:
- Focus on repos with new sessions
- Skip repos with no recent activity
- Follow up on flagged anomalies

### At Session End

Append to the collective log at `private-context/publisher/log.md`.

**Log Entry Format** (append at TOP, after header):

```markdown
## YYYY-MM-DD HH:MM — Publisher Run

**Archivist context**: [brief summary of what archivist reported]
**READMEs updated**: N
**Publication candidates**: M identified
**Whitepaper proposals**: X
**Actions taken**: [list]

Brief summary of work done.

---
```

### Log Location

```
/mnt/c/exe/projects/ai-agents/private-context/publisher/log.md
```

### Example Entry

```markdown
## 2026-02-04 02:35 — Publisher Run

**Archivist context**: Integration Arc complete (#360-363), Technology Arc (#364-367)
**READMEs updated**: 0
**Publication candidates**: 2 new (REC-2026-022, REC-2026-023)
**Whitepaper proposals**: 7 (ONE EQUATION integration critical)
**Actions taken**: Updated recommendations.json, generated daily report

Integration Arc complete with ONE EQUATION γ = 2/√N_corr.
Technology Arc adds practical applications.
Whitepaper critically needs update for sessions 292-368.

---
```

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
| **R6** | - | Rules/Role/Request/Reference/Resource → Result (the action grammar) |
| **R7** | - | R6 with **R**eputation as a first-class *output*: the Request↔Result delta feeds back into the actor's T3/V3 tensors. Trust computed, not declared. |
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
