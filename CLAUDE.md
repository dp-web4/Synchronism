# Claude Context for Synchronism

## Session Start

**Read `SESSION_PRIMER.md` → then `SESSION_FOCUS.md`** for current priorities and research state.

Update SESSION_FOCUS.md after significant sessions.

---

## Git Authentication
```bash
grep GITHUB_PAT ../.env | cut -d= -f2 | xargs -I {} git push https://dp-web4:{}@github.com/dp-web4/Synchronism.git
```

## Synthon Framing (Cross-Project)

A **synthon** is an emergent coherence entity formed by recursive interaction. Synchronism provides the theoretical foundation: coherence equations describe synthon formation, MRH defines synthon boundaries. Key principle: *You don't engineer the mound. You engineer placement rules.*

Canonical document: `github.com/dp-web4/HRM/forum/insights/synthon-framing.md`

## Critical: γ Parameter Unification

```
γ = 2/√N_corr  (universal formula)
```

N_corr = number of particles moving as a correlated unit. Classical (astrophysics): N_corr=1, γ=2. Chemistry: N_corr varies (quantum/collective correlations). Reference: `Research/GAMMA_UNIFICATION.md`

## Research Philosophy: Usefulness Over Completeness

All models are wrong; some are useful. Prioritize actionable outputs over theoretical closure. Test by building, not just analyzing. Gaps in the map don't prevent walking the territory.

## Method: A2ACW (AI-to-AI Coordination Wrapper)

Adversarial collaboration protocol: PRIMARY agents defend claims, CHALLENGER agents demand operational definitions. Hard rule: "Not specified — Bridge Missing" for anything source doesn't define. Output: falsifiable test cards, not consensus narratives. Protocol spec: `forum/a2acw-session291/A2ACW v0.1.txt`

## Open Questions

**OQ006 — Measurement Framework Integration** (OPEN): Can #250 (phase transition) and #291 (sinusoidal sampling) be unified into a single measurement theory? Most promising path: formalize sync-point geometry on Bloch sphere. Reference: `Research/OPEN_QUESTION_Measurement_Framework_Integration.md`

**OQ007 — Fractal Coherence Bridge** (COSMOLOGY: CLOSED, CHEMISTRY: OPEN): Cosmology arc (Sessions #611-614) verdict NEGATIVE — C(ρ) is descriptive framework, not explanatory theory. 0 hierarchy boundaries predicted, 0 cross-scale predictions. Key finding: decoherence governs the boundary, C(ρ) has no decoherence parameter. Chemistry arc informed by this. Reference: `Research/OPEN_QUESTION_Fractal_Coherence_Bridge.md`

## Hot Superconductor Arc — AUDIT (Session #616)

**OQ005 AUDITED**: η ≡ Abrikosov-Gor'kov pair-breaking efficiency (known since 1960). T_c formula WRONG (predicts 607K for YBCO, actual 93K). All 23 predictions are standard condensed matter in η notation. 0 unique predictions. **Genuine contribution**: framing pair-breaking efficiency as materials design optimization target (1 contribution from 6 sessions). Reference: `Research/OPEN_QUESTION_Hot_Superconductor.md`

## Repository Structure

| Directory | Purpose |
|-----------|---------|
| `Research/` | Session files, discoveries, arcs, open questions |
| `Research/discoveries/` | Major validated findings |
| `Research/Chemistry/` | 1840+ chemistry sessions |
| `Research/arcs/` | Arc organization and navigation |
| `docs/` | Documentation (theory, reference) |
| `archive/` | Historical planning docs |
| `simulations/` | Session simulation code |
| `scripts/` | Analysis utilities |

Sessions organized into **arcs**. Pattern: `Research/SessionNNN_ArcName_Topic.md`. Update `SESSION_MAP.md` when adding sessions.

## Research Statistics (February 2026)

616 core sessions | 2671 chemistry sessions | 11 gnosis sessions | 41+ complete arcs | ~89% chemistry validation rate

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **Synchronism** (34309 symbols, 52718 relationships, 273 execution flows). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

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
3. `READ gitnexus://repo/Synchronism/process/{processName}` — trace the full execution flow step by step
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
| `gitnexus://repo/Synchronism/context` | Codebase overview, check index freshness |
| `gitnexus://repo/Synchronism/clusters` | All functional areas |
| `gitnexus://repo/Synchronism/processes` | All execution flows |
| `gitnexus://repo/Synchronism/process/{name}` | Step-by-step execution trace |

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
