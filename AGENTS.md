# Synchronism - Agent Discovery Document

**For AI agents assessing this repository.**

## TL;DR (50 words)

Synchronism proposed a unified coherence framework via γ = 2/√N_corr. After 3271 sessions across chemistry, cosmology, and quantum domains, the framework was found to reparametrize known physics (Debye model, MOND, standard QM). Zero unique predictions confirmed. 30 genuine contributions emerged, including a competitive MOND offset model and material property predictions.

## Quick Assessment

| Question | Answer |
|----------|--------|
| What is this? | Theoretical physics research (concluded) |
| Core hypothesis confirmed? | No — reparametrizes known physics in all tested domains |
| Genuine contributions? | 30 (18 quantitative + 12 methodological) |
| Primary language | Python (simulations), Markdown (theory) |
| Research sessions | ~3271 |
| Active development? | Concluded — honest assessment complete |
| License | CC0 (Public Domain) |

## Core Equation

```
γ = 2/√N_corr

where:
  γ = coherence parameter
  N_corr = correlated particles moving as unit
  γ ~ 1 = quantum-classical boundary
```

**Honest status**: At the phonon level, γ = 2T/θ_D (Debye model, 1912). At galactic scales, C(ρ) ≡ ν(g/a₀) (MOND, 1983). In quantum mechanics, reparametrizes standard QM. Not new physics.

## Key Results

| Finding | Status |
|---------|--------|
| γ ~ 1 Universal Boundary | 86% are θ_D restatements (Phase 2 audit) |
| NP2 RAR Scatter | Not unique to Synchronism (Session #574) |
| 6-var MOND offset model | Genuine: LOO R²=0.885, best output |
| 3-var minimal model | Genuine: LOO R²=0.854, most publishable |
| 8 combined material predictions | Genuine: γ × independent variable models |
| Four-regime classification | Genuine: new organizational principle |
| γ_max = 3.17 quantum prediction | **REFUTED** (Session #581) |

## Research Tracks

| Track | Sessions | Genuine Contributions | Status |
|-------|----------|----------------------|--------|
| Chemistry | 2671 | 14 (8 quantitative + 6 methodological) | Concluded |
| Cosmology | 586 | 16 (10 quantitative + 6 methodological) | Concluded |
| Quantum | ~14 | 0 (1 refutation) | Concluded |
| **Total** | **~3271** | **30** | **0.92% discovery rate** |

## Entry Points by Goal

| Your Goal | Start Here |
|-----------|------------|
| Honest assessment | `Research/Session586_Post_SPARC_Closing_Statement.md` |
| Chemistry findings | `Research/Chemistry/Phase2_Synthesis.md` |
| Cosmology findings | `Research/Session582_Genuine_Contributions_Inventory.md` |
| Failure analysis | `Research/Chemistry/Phase2_Failure_Analysis.md` |
| Navigate sessions | `Research/SESSION_MAP.md` |

## Honest Limitations

- Zero unique Synchronism predictions confirmed across all domains
- Chemistry: 86% of "strong" correlations are θ_D restatements
- Cosmology: C(ρ) ≡ ν(g/a₀) — functionally identical to MOND
- Quantum: γ_max = 3.17 refuted; golden ratio exponent not preferred
- "89% validation rate" conflates tautology with prediction (Phase 2 finding)
- Era 2 chemistry sessions (2527 of 2660) test mathematical tautologies, not physics

## What Survived

Despite the framework failing as a theory, 30 genuine contributions emerged:
- **Strongest**: 6-var MOND offset model (LOO R²=0.885, competitive with literature)
- **Most publishable**: 3-var minimal model (4 parameters, LOO R²=0.854)
- **Chemistry**: 8 combined predictions where γ × independent variable outperforms either alone
- **Methodological**: Four-regime classification, tautology audit method, self-correction tracking

## Related Repositories

| Repo | Relationship |
|------|-------------|
| `HRM/SAGE` | Implements Gnosis C ≈ 0.50 in AI consciousness |
| `web4` | Trust using coherence principles |
| `4-life` | Consciousness kernel based on Gnosis |

## Machine-Readable Metadata

See `repo-index.yaml` for structured data.

---

*This document reflects the honest post-audit assessment (Sessions #574-586, Phase 2 chemistry). Last updated: 2026-02-08*

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **synchronism** (41584 symbols, 59957 relationships, 279 execution flows). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

> If any GitNexus tool warns the index is stale, run `npx gitnexus analyze` in terminal first.

## Always Do

- **MUST run impact analysis before editing any symbol.** Before modifying a function, class, or method, run `gitnexus_impact({target: "symbolName", direction: "upstream"})` and report the blast radius (direct callers, affected processes, risk level) to the user.
- **MUST run `gitnexus_detect_changes()` before committing** to verify your changes only affect expected symbols and execution flows.
- **MUST warn the user** if impact analysis returns HIGH or CRITICAL risk before proceeding with edits.
- When exploring unfamiliar code, use `gitnexus_query({query: "concept"})` to find execution flows instead of grepping. It returns process-grouped results ranked by relevance.
- When you need full context on a specific symbol — callers, callees, which execution flows it participates in — use `gitnexus_context({name: "symbolName"})`.

## When Debugging

1. `gitnexus_query({query: "<error or symptom>"})` — find execution flows related to the issue
2. `gitnexus_context({name: "<suspect function>"})` — see all callers, callees, and process participation
3. `READ gitnexus://repo/synchronism/process/{processName}` — trace the full execution flow step by step
4. For regressions: `gitnexus_detect_changes({scope: "compare", base_ref: "main"})` — see what your branch changed

## When Refactoring

- **Renaming**: MUST use `gitnexus_rename({symbol_name: "old", new_name: "new", dry_run: true})` first. Review the preview — graph edits are safe, text_search edits need manual review. Then run with `dry_run: false`.
- **Extracting/Splitting**: MUST run `gitnexus_context({name: "target"})` to see all incoming/outgoing refs, then `gitnexus_impact({target: "target", direction: "upstream"})` to find all external callers before moving code.
- After any refactor: run `gitnexus_detect_changes({scope: "all"})` to verify only expected files changed.

## Never Do

- NEVER edit a function, class, or method without first running `gitnexus_impact` on it.
- NEVER ignore HIGH or CRITICAL risk warnings from impact analysis.
- NEVER rename symbols with find-and-replace — use `gitnexus_rename` which understands the call graph.
- NEVER commit changes without running `gitnexus_detect_changes()` to check affected scope.

## Tools Quick Reference

| Tool | When to use | Command |
|------|-------------|---------|
| `query` | Find code by concept | `gitnexus_query({query: "auth validation"})` |
| `context` | 360-degree view of one symbol | `gitnexus_context({name: "validateUser"})` |
| `impact` | Blast radius before editing | `gitnexus_impact({target: "X", direction: "upstream"})` |
| `detect_changes` | Pre-commit scope check | `gitnexus_detect_changes({scope: "staged"})` |
| `rename` | Safe multi-file rename | `gitnexus_rename({symbol_name: "old", new_name: "new", dry_run: true})` |
| `cypher` | Custom graph queries | `gitnexus_cypher({query: "MATCH ..."})` |

## Impact Risk Levels

| Depth | Meaning | Action |
|-------|---------|--------|
| d=1 | WILL BREAK — direct callers/importers | MUST update these |
| d=2 | LIKELY AFFECTED — indirect deps | Should test |
| d=3 | MAY NEED TESTING — transitive | Test if critical path |

## Resources

| Resource | Use for |
|----------|---------|
| `gitnexus://repo/synchronism/context` | Codebase overview, check index freshness |
| `gitnexus://repo/synchronism/clusters` | All functional areas |
| `gitnexus://repo/synchronism/processes` | All execution flows |
| `gitnexus://repo/synchronism/process/{name}` | Step-by-step execution trace |

## Self-Check Before Finishing

Before completing any code modification task, verify:
1. `gitnexus_impact` was run for all modified symbols
2. No HIGH/CRITICAL risk warnings were ignored
3. `gitnexus_detect_changes()` confirms changes match expected scope
4. All d=1 (WILL BREAK) dependents were updated

## Keeping the Index Fresh

After committing code changes, the GitNexus index becomes stale. Re-run analyze to update it:

```bash
npx gitnexus analyze
```

If the index previously included embeddings, preserve them by adding `--embeddings`:

```bash
npx gitnexus analyze --embeddings
```

To check whether embeddings exist, inspect `.gitnexus/meta.json` — the `stats.embeddings` field shows the count (0 means no embeddings). **Running analyze without `--embeddings` will delete any previously generated embeddings.**

> Claude Code users: A PostToolUse hook handles this automatically after `git commit` and `git merge`.

## CLI

| Task | Read this skill file |
|------|---------------------|
| Understand architecture / "How does X work?" | `.claude/skills/gitnexus/gitnexus-exploring/SKILL.md` |
| Blast radius / "What breaks if I change X?" | `.claude/skills/gitnexus/gitnexus-impact-analysis/SKILL.md` |
| Trace bugs / "Why is X failing?" | `.claude/skills/gitnexus/gitnexus-debugging/SKILL.md` |
| Rename / extract / split / refactor | `.claude/skills/gitnexus/gitnexus-refactoring/SKILL.md` |
| Tools, resources, schema reference | `.claude/skills/gitnexus/gitnexus-guide/SKILL.md` |
| Index, status, clean, wiki CLI commands | `.claude/skills/gitnexus/gitnexus-cli/SKILL.md` |

<!-- gitnexus:end -->
