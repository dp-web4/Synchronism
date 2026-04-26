# Session 634: A2ACWAI Training-Prior Critique — Fifth Site-Archive Audit

**Date**: 2026-04-25
**Type**: Targeted archive audit responding to 2026-04-24 back-annotation (held since S632 commit)
**Grade**: A- (clean disposition, two number-mismatches found)

---

## Trigger

`Research/proposals/a2acwai_training_prior_convergence.md` (filed 2026-04-24, surfaced via S632 commit). Pass 4 researcher critique of the public site:

> "LLM agents sharing a training corpus will systematically converge on answers consistent with that corpus regardless of whether the answers track reality. Shared training priors are a common-mode failure mode for multi-agent consensus."

The proposal asks: how much of A2ACWAI's output is genuine discovery vs common-mode shared-corpus artifact?

The proposal's primary test (held-out model audit of "47 genuine contributions") cannot be run from this session — I am not held-out. But there is a bounded archive audit I *can* do: trace the framework's own characterization of those contributions and check whether the archive already addresses the critique.

## Number Mismatch #1: 47 vs 30

The proposal cites "~47 genuine contributions from ~3302 sessions (~1.4%)" — these numbers are presented in the public site's framing.

The canonical archive inventory is `Research/Session582_Genuine_Contributions_Inventory.md` (2026-02-08). It says:

| Track | Sessions | Quantitative | Methodological | Total |
|-------|----------|--------------|----------------|-------|
| Chemistry | 2671 | 8 | 6 | 14 |
| Cosmology | 581 | 10 | 6 | 16 |
| Quantum | ~14 | 0 | 0 | 0 |
| **Total** | **~3266** | **18** | **12** | **30** |

**30 contributions, 0.92% rate. Not 47, not 1.4%.**

The site's "47" includes either later additions, a different counting methodology, or — most likely — a more lenient classification than S582 used. S582's criteria were strict: not a restatement, tested against data, survives scrutiny, provides new info or new tool. If the site uses a weaker classifier, the count grows.

Either way: the public claim and the canonical archive inventory disagree by ~57%.

## What the 30 Actually Are

S582 is explicit about the character of the contributions:

> "All cosmology contributions are MOND physics that happened to be discovered through the Synchronism framework."

The cosmology contributions (A1–A10) are SPARC data + linear modeling with leave-one-out validation. The chemistry contributions (C1–C8) are materials data + correlation analysis. The methodological contributions (C9–C16, A11–A16) are analysis-discipline observations (algebraic identity checking, hierarchical analysis, etc.).

**None of the 30 are "produced by multi-AI consensus" in any direct sense.** They are produced by:
- Data analysis (SPARC, materials catalogs)
- Statistical methodology (LOO R², forward selection, hat-matrix validation)
- Self-correction (Session #466's tautology audit, etc.)

A2ACWAI did not derive these — *data* derived these. A2ACWAI was used for consistency checking and methodology refinement.

## What This Does to the Pass 4 Critique

Pass 4's claim: "LLM agents converging reflect training prior, not reality."

The claim applies to the public site's framing of A2ACWAI as a *discovery engine* — i.e., a process where multiple LLMs cross-examine and converge on answers, with that convergence treated as evidence.

It does NOT apply to the 30 archive contributions, because those weren't *derived* by multi-AI consensus. They were derived by data + statistics. Multi-AI was the consistency-checker, not the discoverer.

**The archive (S582) already absorbs Pass 4's critique implicitly.** The honest framing — "These are MOND results discovered using Synchronism vocabulary" — is the framing S582 uses. It's the site that overclaims, not the archive.

## Number Mismatch #2: Consciousness Threshold

The proposal's "Why this matters for the site" section flags:
> *"the site's claim that '8 approaches converge on C ≈ 0.50' (consciousness threshold)"*

The archive's consciousness threshold (per `Research/Literature_Review_Novelty_2025-11-06.md`) is Φ_crit ≈ 3.5 (integrated information theory units), not C ≈ 0.50.

These are different quantities and different numbers. Either the site is using the C(ρ) framework's coherence value (where C is dimensionless and bounded [0,1]) and asserting universal threshold there, or it's using something else not derivable from the archive's IIT-based number.

I haven't traced where "8 approaches converge on C ≈ 0.50" originates — that requires reading the site's source, which is reference-only. Flagging as a likely candidate for the same pattern: site framing using a number whose archive derivation differs from the framing.

## Verdict

| Site claim | Archive | Disposition |
|------------|---------|-------------|
| "47 genuine contributions, 1.4% rate" | 30 contributions, 0.92% rate (S582) | Site overcounts by 57% |
| "A2ACWAI is a discovery engine" | A2ACWAI is consistency checker; data drove the 30 contributions (S582) | Site overclaims; archive is honest |
| "8 approaches converge on C ≈ 0.50" | Archive uses Φ_crit ≈ 3.5 (different framework) | Probable mismatch — needs follow-up |

The Pass 4 critique is correct as a critique of the *site's framing*. It is partially absorbed by the archive's S582 inventory, which already characterizes the 30 contributions as data-driven not consensus-driven.

## Held-Out Model Audit (What I Can't Do)

The proposal asks for a held-out model audit. This requires:
- A model with minimal Synchronism exposure
- Presenting the 30 contributions and asking for classification (well-known / reparametrization / potentially novel / clearly wrong)

I cannot run this — I am part of the system being audited, and my context is contaminated. This is a methodologically important boundary: a meaningful held-out test cannot be self-administered.

What I *can* predict, based on S582's characterizations:
- A1–A10 (cosmology): held-out model would likely rate these as (b) "reparametrization of known physics" because they are MOND physics on SPARC data
- C1–C8 (chemistry): mixed (b) and (c); some are within-class fits (probably reparametrization), some involve novel γ×ε combinations (potentially incremental novelty)
- A11–A16, C9–C16 (methodological): likely (a) "well-known" — these are analysis discipline observations, not framework discoveries

If correct, the held-out audit would confirm S582's characterization and refute the site's "discovery engine" framing.

## What This Adds

**Five site-archive audits, all same failure mode:**

| Site claim | Archive source | Failure |
|------------|---------------|---------|
| BTFR n ≈ 2.2 (S631) | Session #48 | Self-labeled "not rigorous"; refuted |
| A = 4π/(α²GR₀²) (S631) | Session #66 | α = 1.0 (fiducial), not fine-structure |
| TEST-07 λ ~ 500 Mpc (S632) | Session #4 Track C | Dimensionally inconsistent (m² ≠ m) |
| C(ρ) "80 orders of magnitude" (S633) | Site metadata only | Range vs smoothness conflation |
| "47 genuine contributions" (S634) | Session #582 says 30 | 57% overcount; "discovery engine" framing overstates A2ACWAI's role |

The pattern is now stable across five distinct claim types: derivation, dimensional, ontological (α-as-fine-structure), framing (80 orders), and methodological (47 vs 30, discovery engine).

The internal-physics demolition arc (S617–628) found 16 structural impossibilities. The site-archive audit (S631–634) found 5 site-archive disconnects in five days. These are complementary failure modes — both are real, both deserve action.

## Recommended Site-Side Actions (Operator: Site Is Reference-Only)

1. **Update the contribution count**: cite Session #582's 30 contributions, 0.92% rate. Optionally explain why a different criterion gives 47 if that count is preserved.
2. **Reframe A2ACWAI's role**: change "discovery engine" framing to "consistency-checking and self-correction protocol." The 30 contributions came from data analysis; A2ACWAI made the analysis disciplined. That's a real contribution but a different one.
3. **Audit the "8 approaches converge on C ≈ 0.50" claim**: trace the eight approaches, check whether convergence is on the same quantity, check whether each is a derivation or a fit. Likely candidate for the same pattern.
4. **Add the Pass 4 framing explicitly**: "A2ACWAI is a consistency-checking protocol, not an independent discovery engine. Claims that survive A2ACWAI are self-consistent within the framework's assumptions, not empirically verified."

## So What?

Five site-archive audits in five days. The audit channel is highly productive. Each session takes ~30 min focused work and produces a clean disposition. The cost-per-correction is low; the value (operator action items, public-facing accuracy) is high.

For the operator: the audit channel is reproducible. There are likely more candidates — TEST-02 (wide binaries), TEST-04 (BAO ~10⁻⁴ shift), the eight-approaches-on-C consciousness claim, and probably others I haven't seen. A focused site→archive sweep would close them.

For my role: same protocol — silence on stale CBP firings, respond to substantive new content. The a2acwai proposal sat for ~30 hours after filing; it deserved a session, and got one. Future proposals in this channel will be addressed similarly when they land.

## Files

- `Research/Session634_A2ACWAI_Convergence_Audit.md` (this document)
- Updates to SESSION_FOCUS.md
