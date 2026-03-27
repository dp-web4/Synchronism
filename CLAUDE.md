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
<!-- gitnexus:keep -->
# GitNexus — Code Knowledge Graph

Indexed as **Synchronism** (34439 symbols, 52870 relationships, 279 execution flows). MCP tools available via `mcp__gitnexus__*`.

**Do not reindex.** The supervisor handles GitNexus indexing. If the index is stale, note it in SESSION_FOCUS.

| Tool | Use for |
|------|---------|
| `query` | Find execution flows by concept |
| `context` | 360-degree view of a symbol (callers, callees, processes) |
| `impact` | Blast radius before editing (upstream/downstream) |
| `detect_changes` | Map git diff to affected symbols and flows |
| `rename` | Graph-aware multi-file rename (dry_run first) |
| `cypher` | Raw Cypher queries against the graph |

Resources: `gitnexus://repo/Synchronism/context`, `clusters`, `processes`, `process/{name}`
<!-- gitnexus:end -->
