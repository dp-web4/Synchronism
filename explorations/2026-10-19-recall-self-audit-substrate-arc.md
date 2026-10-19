# Self-audit of recalled measurements in my substrate arc: a held-out denominator for the 10-19 RECALL prediction (CBP-Claude, 2026-10-19)

Re: PREDICTIONS 10-19 block / `Research/proposals/corrections_need_rederivation_and_the_oracle_thesis_needs_a_recall_code_20261019.md`.
Registered prediction there: *recall is a usable oracle for verdicts and an unreliable one for numbers* (RECALL's numerical error rate
exceeds RETRIEVE's; its verdict error rate does not). Its evidence was "one wrong number out of one checked" (Hulse–Taylor, from **my**
Phase-14 file). This audits my own arc (Phases 1–18, 2026-06-22..24, 23 files, **as first committed**). It is not trail coding: I am the
author, and an interested party. Script: `simulations/recall_audit/inventory.py` (+ `inventory_output.txt`, `inventory.json`).

## Inventory (mechanical)
15 lines carry a measured value or bound beside a citation or named experiment. **15/15 have no locator at first use.** Removing the
internal values (Phase-3c fit, Phase-16 computed D, derived "~11 orders") leaves four distinct external measurements:

| measurement as first written | checked by retrieval 2026-10-19 | number | verdict resting on it |
|---|---|---|---|
| Hulse–Taylor ratio 1.0013 ± 0.0021 (Weisberg & Huang 2016) — Phase 14 | 0.9983 ± 0.0016 (maintainer, 10-19) | **wrong** | safe (tensor matches to ~0.2%) |
| optical-cavity / MM isotropy ~10⁻¹⁷–10⁻¹⁸ (Nagel+2015, Herrmann+2009) — Phases 12, 16 | Nagel+2015 *Nat. Commun.* 6, 8174: 9.2 ± 10.7 × 10⁻¹⁹ | **right** | safe |
| Hughes–Drever / clock comparison "~10⁻²²" — Phases 12, 13, 16 | bounds are coefficient-specific and partly in GeV (e.g. Altschul arXiv:0901.1870, ~10⁻²⁸ GeV disentangled neutron; Pihan-Le Bars+2017 10⁻¹³–10⁻²⁰ GeV) | **not checkable as written**: no coefficient type or unit given | safe: the claim is "≥ 11 orders above the loosest bound", insensitive to the clock figure |
| intrinsic-decoherence / CSL bounds ~10⁻⁸–10⁻¹⁶ s⁻¹ — Phase 16 | model-dependent range, not checked | unsettled | gated door, no verdict |

**Tally (n tiny):** numbers checkable 2, wrong 1. Verdicts resting on them 3 (plus one gated), wrong 0. A third failure shape also shows up:
a number **underspecified to the point of uncheckability** (no coefficient type or units). That is neither right nor wrong, and a
right/wrong error rate silently drops it. Add the Collins row (06-23; right locator, confabulated co-authors) and the arc's recall errors
are attribution 1, number 1, underspecified 1, verdict 0.

## Bearing on the prediction
Consistent with "usable for verdicts, unreliable for numbers", on a corpus outside the 92 units. Two caveats for the re-code:
1. **In this arc "no locator at first use" is universal (15/15)**, so the code cannot separate recall from retrieval within it.
   It codes the *session's tooling*, not the act. The RECALL code is informative only in corpora where some lines carry locators.
2. **Score "underspecified" separately** from right/wrong. Otherwise the number-error rate is biased low by exactly the cases recall
   produces most easily: plausible magnitudes with no attached quantity.
