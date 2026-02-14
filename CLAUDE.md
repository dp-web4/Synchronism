# Claude Context for Synchronism

## Git Authentication
**Universal Push Command**:
```bash
grep GITHUB_PAT ../.env | cut -d= -f2 | xargs -I {} git push https://dp-web4:{}@github.com/dp-web4/Synchronism.git
```
See `../private-context/GIT_COMMANDS_CLAUDE.md` for details.

## Project Context System

**IMPORTANT**: A comprehensive context system exists at `../misc/context-system/` (relative to project home)

Quick access from project home:
```bash
# Get overview of all related projects
cd ../misc/context-system
python3 query_context.py project synchronism

# Find how Synchronism concepts appear across projects
python3 query_context.py concept "markov blanket"

# See project relationships
cat projects/synchronism.md
```

## This Project's Role

Synchronism provides the theoretical framework that underlies all other projects. Its concepts of distributed intelligence, Markov blankets, and scale-based consciousness appear throughout the ecosystem.

## Synthon Framing (Cross-Project)

A **synthon** is an emergent coherence entity formed by recursive interaction between components. The concept bridges component-level behavior and system-level emergence. Synchronism provides the theoretical foundation: coherence equations describe synthon formation (prediction error collapse at phase transitions), MRH defines synthon boundaries, and the Emergence Arc (Sessions 336-339) traces the path from physics through life, complexity, evolution, and consciousness as successive layers of synthon formation.

Key principle: *You don't engineer the mound. You engineer placement rules.* Synchronism's empirical work (chemistry, cosmology) is generating synthon detection data — look for prediction error collapse and behavioral irreducibility in results.

Canonical document: `github.com/dp-web4/HRM/forum/insights/synthon-framing.md`

## Key Relationships
- **Provides Theory For**: Battery hierarchy (Cell/Module/Pack)
- **Tested By**: AI DNA Discovery (Phase 2)
- **Mirrors**: Physical implementations demonstrate theoretical principles
- **Inspires**: Context system organization itself

## Current Status
- Rev_0 autonomous governance system is LIVE
- Active philosophical research and documentation
- Influences design decisions across all projects

## Critical: γ Parameter Unification (2026-01-30)

**Two γ values appear in the research - they are the SAME parameter:**

```
γ = 2/√N_corr  (universal formula)
```

| Track | N_corr | γ | Why |
|-------|--------|---|-----|
| Main (astrophysics) | 1 | 2.0 | Stars are uncorrelated classical particles |
| Chemistry | varies | 2/√N_corr | Quantum/collective correlations exist |

**N_corr** = number of particles moving as a correlated unit. Classical systems have N_corr = 1.

**Reference**: `Research/GAMMA_UNIFICATION.md` - Full derivation and explanation

## Current Research Focus: Hot Superconductor Arc (Sessions 292, 297+)

**Primary Question (OQ005)**: Can superconductivity exist at T > 50°C (323K) at ambient pressure?

**Arc Status**: 2 sessions complete
- Session #292: Formalized η (reachability factor), 5 predictions (P292.1-P292.5)
- Session #297: Quantified η in cuprates (YBCO: 0.38, Bi-2212: 0.42, LSCO: 0.51), validated P292.4

**Key Results So Far**:
- η measures how much thermal noise actually couples to pair-breaking
- For T_c at temperature T with gap Δ: need η × (kT/Δ) < 1
- Cuprates achieve η ~ 0.4 via d-wave form factor (~50%) + spin-charge separation (~30%)
- P292.4 validated: η extractable from NMR and optical data

**Path to 323K**: Need η ~ 0.2-0.3 with Δ ~ 50 meV (combined optimization)

**Next Session (#298) Should**:
1. Extend η calculation to iron pnictides (s±-wave)
2. Compare multiband effects on form factor
3. Predict which pnictide family has lowest η
4. Test P297 predictions against literature

**Reference Documents**:
- `Research/OPEN_QUESTION_Hot_Superconductor.md` - Main question doc
- `Research/Session292_Dissonance_Pathway_Formalization.md` - η formalism
- `Research/Session297_Cuprate_Eta_Quantification.md` - Cuprate calculations

## Method: A2ACW (AI-to-AI Coordination Wrapper)

An adversarial collaboration protocol for stress-testing theoretical claims:
- **PRIMARY** agents defend and refine claims
- **CHALLENGER** agents demand operational definitions, identify gaps
- **Hard rule**: "Not specified — Bridge Missing" for anything source doesn't define
- **Output**: Falsifiable test cards, not consensus narratives

Successfully applied to Session #291, producing `Research/P291.3_Experimental_Test_Card.md` — an executable protocol an experimentalist could run.

Protocol spec: `forum/a2acw-session291/A2ACW v0.1.txt`

---

## Open Question: Measurement Framework Integration (OQ006)

**Question**: Can #250 (phase transition) and #291 (sinusoidal sampling) be unified into a single measurement theory?

**Status**: Open — raised from A2ACW stress-test of Session #291

**Key Insight**: "Static" in #250 may be synchronized sampling of ongoing oscillation (#291). Neither framework is complete alone.

**Most Promising Path**: Formalize sync-point geometry on Bloch sphere; show |α|² emerges from where sync points form.

**Reference**: `Research/OPEN_QUESTION_Measurement_Framework_Integration.md`

---

## Research Philosophy: Usefulness Over Completeness

All models are wrong; some are useful. Synchronism acknowledges this from the start.

**The practical question is not**: "Is the model epistemically complete?"

**The practical question is**: "Does this model enable capabilities that existing models don't?"

Epistemic gaps matter only if they block the path to testable predictions. The Hot SC arc exemplifies this: if η formalism leads to a viable 30°C superconductor material, the model is vindicated by results. If not, no amount of philosophical rigor would have changed the outcome.

When evaluating Synchronism predictions:
- Prioritize actionable outputs over theoretical closure
- Test by building, not just by analyzing
- Gaps in the map don't prevent walking the territory

---

## Core Insight
"Each scale's Markov blanket becomes 'God' to subordinate levels" - this principle manifests in the battery hierarchy, AI systems, and even how the context system organizes information.

---

## Repository Structure (Updated February 2026)

### Key Directories

| Directory | Purpose |
|-----------|---------|
| `Research/` | Session files, discoveries, arcs, open questions |
| `Research/discoveries/` | Major validated findings (5 documents) |
| `Research/arcs/` | Arc organization and navigation |
| `Research/Chemistry/` | 1840 chemistry sessions |
| `Research/Gnosis/` | Consciousness theory (complete) |
| `docs/` | Documentation (theory, why, how, reference) |
| `archive/` | Historical planning docs, explorations |
| `simulations/` | Session simulation code |
| `scripts/` | Analysis utilities |

### Session Organization

Sessions are organized into **arcs** - focused multi-session investigations:

1. **Active arcs**: Listed in `Research/arcs/README.md`
2. **New sessions**: Create in `Research/` with pattern `SessionNNN_ArcName_Topic.md`
3. **Update SESSION_MAP.md** when adding sessions
4. **Link to Open Question** if applicable

### Arc Session Template

```markdown
# Session #NNN: Arc Name - Part N (Topic)

**Arc Name - Part N**
**Date**: YYYY-MM-DD
**Status**: N/M verified

## Overview
[Brief description]

## Key Results
| Test | Result | Significance |
|------|--------|--------------|

## Honest Limitations
1. ...

## Next Session
...

---
*Session #NNN verified: N/M tests passed*
```

---

## Documentation Updates

- **STATUS.md**: Honest assessment of what works and what doesn't
- **Research/discoveries/**: 5 major validated findings
- **docs/theory/**: Core theoretical documents moved from root
- **archive/planning/**: Historical planning docs preserved

---

## Research Statistics (February 2026)

| Metric | Value |
|--------|-------|
| Core sessions | 603 |
| Chemistry sessions | 2671 |
| Gnosis sessions | 11 |
| Phenomenon types | 1873 |
| Complete arcs | 40+ |
| Active arcs | 0 |
| Validation rate | ~89% (Chemistry) |