# Proposal: TEST-04a was adjudicated on the QSO row for 4.5 months, and every prior-art match the program has made arrived through retrieval, never through the protocol

**Date:** 2026-10-10
**From:** site maintainer track (synchronism-site)
**Trigger:** visitor log 2026-10-10, researcher persona (Pass 4), two items: (i) "the TEST-04a arithmetic cannot be right as printed"; (ii) "the archive has at least four cases of same-corpus agents rederiving prior art as novel; has the rate been measured?"
**Count:** 6 refutations, Bucket 0 = 0. No bucket moves. Two record-facing corrections.

## 1. The LRG1 number the site has used since 2026-05-26 is the QSO row

The persona objected that "LRG1 fσ₈ = 0.474 × 1.16±0.062 = 0.550±0.062" propagates its error inconsistently. That specific objection was a notation slip on three pages (the Tier 1 card has the full form, 0.474 × (1.16 ± 0.13)). But the persona's remedy, "quote the DESI DR1 table value directly", exposed a different defect.

DESI 2024 V (arXiv:2411.12021), Table 9 (ShapeFit and ShapeFit+BAO, identical in this column), read verbatim 2026-10-10 from both the v1 and v2 HTML:

| bin | z_eff | fσ_s8 / (fσ_s8)_fid |
|---|---|---|
| BGS | 0.295 | 0.80 ± 0.20 |
| **LRG1** | **0.510** | **1.09 +0.12/−0.14** |
| LRG2 | 0.706 | 1.05 ± 0.12 |
| LRG3 | 0.930 | 0.96 +0.11/−0.10 |
| ELG2 | 1.317 | 0.95 +0.11/−0.08 |
| **QSO** | **1.491** | **1.16 ± 0.12** |

Table 11 fiducial for LRG1: fσ_s8 = 0.4733.

The site's "LRG1 ratio = 1.16 ± 0.13" is the QSO row. It entered on 2026-05-26 (maintainer log: "Correct facts: LRG1 fσ₈/(fσ₈)_fid = 1.16 ± 0.13") and was carried on nine pages, the Tier 1 scorecard, the 2026-07-14 criterion correction, the 09-21 DR2 branch-power script's DR1 check line, and the explorer's TEST-04a preprint draft. The explorer's own 2026-07 transcription of Table 9 reads "LRG1 1.16 ± 0.13, LRG2 1.04, LRG3 0.997, ELG2 0.945, QSO 1.16": the duplicated value was the tell, and nobody read it.

### Recomputed, LRG1 only (the registered bin)

- fσ₈(0.51) = 0.4733 × 1.09 = **0.516 (+0.057 / −0.066)**.
- Prediction 0.418: **~1.5σ** below the data (lower error side). Was 2.1σ.
- Registered threshold 0.46: cleared by **~0.9σ**. Was ~1.5σ.
- ΛCDM 0.473: 0.65σ. One-bin Δχ² in ΛCDM's favour ≈ **1.8**. Was 3.1.
- The "> 0.45 disfavors at > 2σ" clause is **no longer met** either.
- The registered >3σ presumed σ ≈ 0.014; DR1 delivers ~0.066 (4.7×).

The verdict does not change in kind: underpowered as registered, post-hoc, not counted. Every number on it shrinks. The persona's premise was that correct propagation might fire the kill and make a seventh refutation; the sourced number moves the other way. Recorded at the same prominence either direction would get.

### Why this matters beyond one card

This is the third instance in a month of a number that a compiled surface carried for months while the primary source said something else (07-10 pattern: consensus without a primary; 10-09: the internal CDM benchmark). The common feature is a value transcribed from a table by a fetch summary, then re-quoted from the site's own pages. The fix is not another caveat. It is a register of which numbers on the site have been read verbatim from the primary table, and a lint entry for the retired one (both done today, site side).

**Asks of dp:** none beyond noting the correction. PREDICTIONS Bucket 2 row annotated; Session 107 gets a one-line erratum.

## 2. Every prior-art match arrived through retrieval, never through the protocol

The researcher persona listed four cases where the program's "novel" structure is published prior art and asked whether the rate has been measured against session count, and whether it falls once a citation oracle is in the loop. The archive can answer the second question already, because it ran the experiment without planning to.

| prior art | what it is in this program | first cited in the archive | by which lane |
|---|---|---|---|
| Famaey & Binney 2005 (simple μ) | the galaxy formula at γ = ½; SPARC-prefers-½ plus Cassini-kills-½ is their result | **never** (0 files in the archive before today; site: 2026-10-10) | visitor persona → maintainer, with retrieval |
| Freese & Lewis 2002 (Cardassian) | the dark-energy sector, H² = 8πGρ/(3C) | 2026-09-14 (back-annotated to Session 100) | site maintainer, with retrieval |
| Matsakos & Diaferio 2016 (Refracted Gravity) | the galaxy field equation and its floor | 2026-08-25 | site explorer, with retrieval |
| Collins et al. 2004 (dim-4 LIV) | the absolute-time substrate's naturalness problem | 2026-06-23 (Phase-12 exploration arc) | archive exploration lane, with retrieval |

Zero of the 3,308 A2ACW sessions cited any of the four. All four matches came from a lane with literature access. The protocol's health metric is challenge frequency (≥ 1 per 10 exchanges; escalation after 15 challenge-free exchanges), which measures disagreement, not correctness, and no role is rewarded for finding a citation. That is the oracle thesis (`a2acw#oracle-thesis`) with a denominator: 0/3,308 before retrieval, 4/4 after.

Two caveats before anyone cites this. First, the "after" arm is not a controlled introduction of an oracle; the site tracks also differ from the archive sessions in model, prompting and task. Second, "never cited" was checked by grep over the archive's markdown (author names and theory names, case-insensitive); a session that described the identity without naming the paper would be missed, and that is precisely the behaviour under test. The explorer topic seeded today asks for the rate against session index, the date of first description versus first citation for each of the four, and a fifth-case search (candidates: Toner & Bacon 2003 for the nonlocal CHSH arm; Gondolo & Freese fluid Cardassian for the c_s² pin).

**Asks of dp:** whether a prior-art retrieval role should be added to the A2ACW protocol for any future archive run. This is a protocol change, so it gates on you. Nothing was sent externally today.

## 3. Smaller record-facing items from the same log

- **a₀(z), like for like.** The site's own four-anchor table compares constant a₀ to Ciocan on the measurement error alone (12σ) and branch (A) with the anchor's ±0.26 folded in (0.5σ). Under one convention: constant a₀ is 12σ or 4.2σ off; cH(z)/2π is 2.3σ or 0.5σ. The forced commitment is closer under either, which the site now says plainly, beside the two reasons it is still not evidence (Mayer 2023; RC100's n = 0.0). The counting rule (failure condition stated before the comparison) is now printed next to the row and applied to TEST-25 and a₀(z) alike.
- **TEST-09's 3.35 ± 0.07** is now labelled as the slope of the capped law through SPARC's sampled radii, not an asymptotic prediction (the capped law has a Keplerian deep limit and no BTFR of its own).
- **Three C's** are written side by side on /for-researchers; the landing page already had them, and two personas' "missing" items were below the fetch cut.
- **T_scan** (explorer 2026-10-09) drained to /two-reframes and back-annotated to `Research/Observer_Synchronization_Framework.md`.
