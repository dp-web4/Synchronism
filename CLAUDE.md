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

## Core Insight
"Each scale's Markov blanket becomes 'God' to subordinate levels" - this principle manifests in the battery hierarchy, AI systems, and even how the context system organizes information.