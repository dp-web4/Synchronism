# A2ACW v2 — A Three-Axis Protocol from the Vocabulary-Asymmetry Result

**Filed**: 2026-05-19 (back-annotated from synchronism-site explorer track)
**Status**: Proposal — A2ACW protocol revision based on a 4/6 catch-rate result
**Source finding**: `synchronism-site/explorer/findings/a2acw-vocabulary-asymmetry-result.md`
**Prior**: 2026-05-18 temporal-asymmetry counterfactual (0/6 catch); memory `project_a2acw_vocabulary_lockin.md`

---

## Summary

The 2026-05-18 temporal-asymmetry counterfactual measured A2ACW's catch rate on its own 6 demoted claims using training-cutoff-asymmetric adversaries. Result: 0/6. The diagnosis at the time was *vocabulary lock-in, not temporal lock-in*.

The successor experiment was run 2026-05-19. Vocabulary-asymmetry catch rate: **4/6 overall, 4/4 on the prior-art-rediscovery sub-class.**

The two misses are not "things vocabulary asymmetry failed at" — they are demotions in *different failure-mode classes* that vocabulary asymmetry was never designed to catch:

| Demotion | Failure-mode class | Correction axis |
|---|---|---|
| Born rule | Prior-art rediscovery (Zurek envariance) | Vocabulary translation + adversary with prior-art knowledge |
| Wide-binary EFE | Prior-art rediscovery (Bekenstein-Milgrom 1984) | (same) |
| Galaxy rotation | Prior-art rediscovery (MOND / RAR) | (same) |
| Γ = γ²(1−c) | Prior-art rediscovery (Palma-Suominen-Ekert 1996) | (same) |
| **Dual-C** | **Internal-consistency** | **Single-source-of-truth audit** |
| **Chemistry r=0.98** | **Null-baseline deficit** | **Explicit null-model computation** |

The protocol should therefore be redesigned as a **three-axis filter**, not a single-axis swap.

## Proposed A2ACW v2

```
A2ACW.v2(claim):
  axis_1 = vocabulary_translation(claim) → modern-register adversary
                                            with prior-art knowledge
  axis_2 = symbol_audit(claim, whole_framework) → flag duplicate or
                                                  inconsistent definitions
  axis_3 = null_model(claim) → fit constructed monotonic / shuffled null
                               and compare AIC/BIC
  flag = OR(axis_1.flagged, axis_2.flagged, axis_3.flagged)
```

Each axis has a specific calibration target:

### Axis 1 — Vocabulary translation
- **Operational requirement**: claims must be **pre-translated to modern register** before adversarial review. The asymmetry is in the *presentation*, not in the *adversary*.
- **Catch target**: prior-art-rediscovery against competing frameworks where canonical references exist (MOND, envariance, DFS, etc.).
- **Calibration**: 4/4 on appropriate sub-class in the 2026-05-19 experiment.
- **Null-baseline requirement**: run on a closed pre-1990s physics corpus (BCS, Anderson localization, EW unification) and report false-novelty rate. Should be ≈ 0 if calibrated; ≫ 0 implies the adversary is hallucinating novelty.

### Axis 2 — Internal-consistency audit
- **Operational form**: AST or grep traversal flagging any symbol re-defined across pages or sessions with non-trivially different formulations.
- **Catch target**: dual-definitions (dual-C, the γ collision documented 2026-05-04, possible BTFR-slope drift).
- **Calibration**: 1/1 on dual-C in the 2026-05-19 experiment. Should also catch the γ regime-label inversion documented in `gamma_definitional_collision_regime_label_inversion.md`.

### Axis 3 — Null-baseline computation
- **Operational form**: for any claim of the form "framework predicts correlation r > T between X and Y," compute the r expected under a structured null where X is replaced by an arbitrary monotonic function of the same underlying variable.
- **Catch target**: forced correlations under monotonicity / order constraints (chemistry r=0.98, possible analogs in atomic-volume and sound-velocity datasets).
- **Calibration**: 1/1 on chemistry in the 2026-05-19 experiment.

Under the combined three-axis protocol, the 6/6 demotion set is caught 6/6.

## Why this matters more than the catch-rate number

The result that pushes the 2026-05-19 finding past a "vocabulary asymmetry works better than temporal" headline is structural:

- **The 6 demoted claims are not the same kind of failure.** They decompose 4-1-1 into three distinct failure-mode classes.
- **No single new asymmetry axis can catch all three.** The temporal/vocabulary search was implicitly assuming that the failure mode was homogeneous (priority rediscovery). It is not.
- **The honest A2ACW v2 is a *taxonomy of failure modes* plus *axis-specific corrections***, not a tweaked adversary.

This is consistent with the framework's broader "compander class diagnosis" pattern: the right fix is often to reclassify, not to recalibrate.

## What remains open

1. **Fresh-adversary validation**: the 4/6 result was self-simulated by one agent (Claude Opus 4.7) that knows both vocabularies. A real fresh adversary (different model family, untouched by the Synchronism corpus) should be tested. Expect catch rate to be *lower* on the prior-art-class — possibly 2/4 if the adversary has weaker prior-art coverage on DFS or envariance.

2. **False-novelty rate on closed physics**: this is the open question raised by the 2026-05-19 visitor Pass 4 ("A2ACW page does not quantify a null-baseline false-novelty rate"). The three-axis v2 should be run on BCS, Anderson localization, and EW unification corpora and the per-axis flagging rate reported. Without this, even a 6/6 protocol is uncalibrated.

3. **Coverage outside the demotion sample**: the 6 demoted claims were *already flagged*. The deeper test is whether v2 would catch *currently-undemoted* claims that should be demoted. Top candidates from the audit ledger:
   - "C(ρ) is a logarithmic compander, not a phase transition" — axis-1 would catch this against the μ-law / Naka-Rushton / Hill literature
   - "ρ_crit = A·V_flat² with 5–12% agreement" — axis-3 would catch this against a null where V_flat is the input (rho_crit_derivation_calibration_vs_prediction proposal)
   - "MRH is Markov" — axis-2 would catch this against the framework's own admission of no probability distribution

4. **Self-application**: the protocol should be applied to itself. Does A2ACW v2 catch A2ACW v1 as a demoted methodology? (Yes — by axis-1, against the broader literature on adversarial validation in ML; by axis-3, against the absent null-novelty baseline.)

## Decision required

Adopt A2ACW v2 as the canonical protocol going forward, or treat the 4/6 result as additional evidence for the "filter, not discovery method" frame and stop trying to fix the protocol?

Both are defensible. The current site posture is the second (`/a2acw` reads as a self-audit of methodology, not a working filter). The vocabulary-asymmetry result strengthens that posture by showing what *would* be needed to make it a working filter. The maintainer pass should propagate this distinction onto the public page.

---

*Back-annotated from synchronism-site/explorer/findings/a2acw-vocabulary-asymmetry-result.md, 2026-05-19.*
