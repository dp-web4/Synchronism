# Synchronism Research Memory

Accumulated research knowledge for the Synchronism project, using SNARC-based salience scoring for intelligent retrieval and pruning.

## Philosophy

This is **semantic memory** - accumulated research knowledge about the domain.
Not autobiographical (what sessions did), but conceptual (what we know).

Pattern: Modeled after SAGE's memory system, but for research findings instead of lived experience.

## Quick Start

```bash
cd memory/tools

# Initialize database (first time only)
python3 init_db.py

# Query findings
python3 query.py "BTFR"                    # Full-text search
python3 query.py --domain galaxy_physics   # By domain
python3 query.py --type derivation         # By type
python3 query.py --high-salience           # High-importance only

# View structured data
python3 query.py --parameters              # All tracked parameters
python3 query.py --predictions --untested  # Untested predictions
python3 query.py --validations             # Test results
python3 query.py --stats                   # Overview

# Add new knowledge
python3 add.py finding --interactive       # Interactive finding entry
python3 add.py parameter --name X ...      # Add parameter
python3 add.py prediction --title "..." ...# Add prediction
python3 add.py validation --test "..." ... # Add test result

# Memory management
python3 prune.py --distribution            # See salience distribution
python3 prune.py --threshold 0.3 --dry-run # Preview low-salience prune
python3 prune.py --decay --dry-run         # Preview time-based decay
```

## Schema

### Findings
Core research knowledge with:
- **Domain**: cosmology, galaxy_physics, coherence_math, implementation, validation, methodology
- **Type**: derivation, validation, prediction, connection, methodology, failure, insight, question
- **Validation Level**: proven, validated, theoretical, speculative, falsified
- **SNARC Scores**: surprise, novelty, arousal, reward, conflict → composite salience

### Parameters
Tracked theory parameters (γ, A, B, δ) with:
- Derived vs empirical values
- Derivation status and session
- Physical meaning

### Predictions
Testable claims with:
- Quantitative claims
- Test methods
- What they discriminate from
- Priority (critical/high/medium/low)
- Status (untested/confirmed/falsified)

### Validations
Test results with:
- Success rates, χ² values
- Parameters used
- Comparison to other models

## SNARC Scoring

Each finding has salience scores (0.0-1.0):

| Score | Meaning | Weight |
|-------|---------|--------|
| Surprise | How unexpected? | 25% |
| Novelty | How new? | 20% |
| Arousal | How important to attend to? | 30% |
| Reward | How valuable? | 10% |
| Conflict | Challenges existing knowledge? | 15% |

Composite salience enables:
- **Retrieval**: High-salience findings surface first
- **Pruning**: Low-salience findings decay/deprecate
- **Focus**: Attention goes where it matters

## For Autonomous Sessions

### At Session Start
```bash
# What do we know about your topic?
python3 query.py "your topic"
python3 query.py --predictions --untested  # Open questions
```

### During Session
If you discover something significant:
```python
from add import add_finding, add_parameter, add_prediction, add_validation

add_finding(
    title="Your Discovery",
    summary="One sentence summary",
    domain="galaxy_physics",
    finding_type="derivation",
    session_id="82",
    surprise=0.8,  # Set SNARC scores based on significance
    novelty=0.7,
    arousal=0.9,
    reward=0.8,
    conflict=0.3
)
```

### At Session End
```bash
# Commit memory changes
cd /path/to/Synchronism
git add memory/knowledge.db
git commit -m "Session #XX: Added findings about ..."
git push
```

## Current Contents

As of Session #81:
- **13 findings** across galaxy physics, coherence math, validation, methodology
- **4 parameters**: γ (derived), B (derived), A (semi-empirical), δ (semi-empirical)
- **3 predictions**: Void BTFR offset tests (critical, high, medium priority)
- **3 validations**: SPARC tests with different parameter sets

## Key Findings Available

Query `python3 query.py --high-salience` to see:
- B Parameter BTFR Breakthrough (Session #78)
- Coherence depends on baryonic density insight
- Void galaxy test methodology
- MOND complementarity connection

## Relation to Private Context

- **This memory**: Public research knowledge about Synchronism physics
- **private-context/memory**: Private cross-project insights, strategic directions

Both can coexist - this is domain-specific, that is cross-cutting.
