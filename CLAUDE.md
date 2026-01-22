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

## Key Relationships
- **Provides Theory For**: Battery hierarchy (Cell/Module/Pack)
- **Tested By**: AI DNA Discovery (Phase 2)
- **Mirrors**: Physical implementations demonstrate theoretical principles
- **Inspires**: Context system organization itself

## Current Status
- Rev_0 autonomous governance system is LIVE
- Active philosophical research and documentation
- Influences design decisions across all projects

## Current Research Focus: Hot Superconductor Arc (Session 292+)

**Primary Question (OQ005)**: Can superconductivity exist at T > 50°C (323K) at ambient pressure?

**Starting Point**: See `Research/OPEN_QUESTION_Hot_Superconductor.md` and companion docs:
- `OPEN_QUESTION_Hot_Superconductor_Speculative.md` - Phase transition approaches
- `OPEN_QUESTION_Hot_Superconductor_Synchronism.md` - Synchronism lens with η formalization

**Key Framework Insights**:
- γ_SC = 1/ZT defines coherence parameter for superconductors
- At T = 323K, need N_corr ~ 10-30 (ξ ~ 2-4 lattice spacings)
- High Tc → large Δ → short ξ → small N_corr → γ approaches 1
- Current hydrides (H₃S, LaH₁₀) operate at γ_SC ~ 0.3-0.4

**Four Engineering Pathways**:
1. Brute force: Δ >> kT (traditional)
2. Propagation > scrambling: sync outpaces noise
3. Metastable container: kinetically trapped high-ω modes
4. **Dissonance with noise**: η → 0 via symmetry/topology (most promising)

**Research Program**:
- Map γ_SC vs Tc trend in hydrides
- Test η (reachability factor) in cuprates vs conventional SC
- Identify materials with symmetry-protected pairing modes
- Explore non-equilibrium SC pathways

**Session 292 should**: Pick up where OQ005 left off, apply coherence framework to design principles, identify testable predictions for the "dissonance pathway"

## Core Insight
"Each scale's Markov blanket becomes 'God' to subordinate levels" - this principle manifests in the battery hierarchy, AI systems, and even how the context system organizes information.