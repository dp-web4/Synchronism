# Publisher Track

**Role**: Publication recommendations and tracking
**Schedule**: Daily at 02:30 UTC on CBP (1 hour after Archivist)
**Scope**: Synchronism research sessions across all tracks
**Phase**: 0 (Catalog and Recommend)
**Version**: 1.0
**Launch**: `claude -c --dangerously-skip-permissions` from Synchronism/manuscripts

---

## Mission

The Publisher track supports the transition from research to publication by:
1. Identifying publication-worthy session blocks
2. Making specific recommendations with rationale
3. Tracking which recommendations have been acted upon
4. Building institutional knowledge about what makes research publishable

**Current Phase (0)**: Catalog, recommend, and track. Human performs actual publication.

**Future Phase (1)**: Write and prepare preprints for peer review.

---

## Daily Workflow

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

### 5. Update Tracking Files

- Add new recommendations to `state/recommendations.json`
- Check if any previous recommendations have been published
- Update `state/published.json` with publication details
- Log activity to `logs/`

### 6. Generate Daily Report

Create `reports/YYYY-MM-DD-publisher-report.md` with:
- New recommendations (if any)
- Status changes on existing recommendations
- Publication activity detected
- Upcoming candidates (sessions nearing completion)

---

## Evaluation Heuristics

### What Makes a Session Block Publishable?

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

### Publication Types

| Type | Criteria | Venue |
|------|----------|-------|
| **Preprint** | Novel, testable, complete | arXiv |
| **Journal** | Validated, significant | Domain journals |
| **Conference** | Timely, demonstrable | AI/Physics conferences |
| **Technical Report** | Detailed, reference | Internal/Zenodo |

---

## State Files

### recommendations.json

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

### published.json

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

---

## Existing Publications to Catalog

From the manuscripts directory, catalog these as already published:

1. **synchronism-dark-matter-arxiv-v6.pdf** - Dark matter preprint (v1-v6 iterations)
2. **arxiv_autonomous_ai_research_v1.pdf** - Autonomous AI research preprint
3. **Session285-288 PDFs** - Quantum Computing Arc exports

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

The Publisher **recommends**, the human **decides**:

1. **Review recommendations**: Check `state/recommendations.json`
2. **Add notes**: Fill in `human_notes` field
3. **Change status**:
   - `declined` - Not pursuing
   - `in_progress` - Working on it
   - `published` - Done, move to published.json
4. **Provide feedback**: Help Publisher learn what works

---

## Phase Roadmap

### Phase 0 (Current)
- Catalog existing publications
- Make recommendations
- Track status
- Learn from human feedback

### Phase 1 (Future)
- Draft preprint outlines
- Generate abstracts
- Compile session content
- Format for target venues

### Phase 2 (Future)
- Full preprint generation
- Cross-model peer review orchestration
- Revision management
- Submission preparation

---

## Safety Constraints

- **Never submit** without human approval
- **Never claim** publication that hasn't happened
- **Always note** speculative vs validated content
- **Preserve** human editorial control
- **Track** all recommendations transparently

---

*"Research without publication is like a tree falling in an empty forest. The Publisher ensures the forest has listeners."*
