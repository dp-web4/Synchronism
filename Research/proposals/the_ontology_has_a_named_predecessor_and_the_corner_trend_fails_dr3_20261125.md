# Proposal: the ontology has a named predecessor ('t Hooft's CAI); the EFE-free corner's trend fails DR3 before DR4

**Date:** 2026-11-25 · **From:** synchronism-site maintainer track (Claude Opus 5.5) · **Gates on:** the project lead
for items 1 and 4; items 2 and 3 are record.
**Sources:** visitor log 2026-11-25 (four simulated personas; Passes 2–4 read full page text by curl, not fetch summaries),
site explorer 2026-11-24 (Draco σ; persona-loop catch rate), and one maintainer run (pre-registered, below).

## 1. The "one move" has a worked-out predecessor that the archive never cites
A researcher persona noted that a deterministic, discrete, absolute-time substrate in which observers are patterns,
escaping Bell by superdeterminism, is 't Hooft's *Cellular Automaton Interpretation of Quantum Mechanics* (Springer FTP
185, 2016; arXiv:1405.1548, abstract read verbatim 2026-11-25: QM "as a tool, not as a theory", the Standard Model as
"a quantum mechanical approach to analyse a system that could be classical at its core", Bell reconciled through
superdeterminism, the collapse and measurement problems "eradicated"). A grep of the archive for "hooft" and "cellular
automaton interpretation" finds only two 1979 gauge-theory citations.

This is the fifth documented prior-art identity, and the first at the level SPINE calls the contribution ("The move
itself is the contribution"). The earlier four were pieces: the simple μ, the Cardassian DE sector, Refracted Gravity's
field equation, and the Collins et al. LIV class. Like them, it reached the record through a reading lane with
retrieval (persona noticed; maintainer retrieved). None came through the A2ACW protocol.

**Stake (maintainer, before reading the book; explorer topic seeded):** no difference between the two ontologies
changes a prediction or a no-go. Patterns map to CAI's ontological states, "measurement is synchronization" maps to
classical pointer states, and absolute time maps to the automaton clock. The MRH describes which correlations an
observer tracks, not the dynamics.

**For the project lead:** (a) add one SPINE sentence naming CAI as the nearest worked-out relative of the one move,
pending the comparison; (b) if the comparison finds no prediction-changing difference, restate the contribution as
the audit trail plus the applied vocabulary, not the move. The site now cites CAI in the landing Bell row and as a row
of the distinct-content register on `/for-researchers`.

## 2. EFE-free corner vs Chae 2023's two published points (record; no count change)
Pre-registered: predictions were committed at site `590d4d4` before the run; 3 of 4 held. Script:
`synchronism-site/maintainer/scripts/efe_free_corner_vs_chae2023_points.py` and `_output.txt`. Inputs are from the
arXiv:2305.04613 abstract: δ = 0.034 ± 0.007 at g_N ≈ 10⁻⁸·⁹¹ and 0.109 ± 0.013 at 10⁻¹⁰·¹⁵ m/s², so the boost
ratio is 1.117 ± 0.025 and 1.43 ± 0.06.

The internal reading (γ = 0.87–0.95, a₀′/γ = 1.0–1.2×10⁻¹⁰) predicts:
- 1.59–1.71 at 10⁻¹⁰·¹⁵, which is +2.6 to +4.7σ;
- 1.013–1.023 at 10⁻⁸·⁹¹, which is −3.7 to −4.1σ (P2's predicted range, 1.03–1.12, failed low);
- a ratio of 1.57–1.67 between the two, against 1.28 ± 0.06: +4.7 to +6.4σ.

So the trend is wrong at both ends in the analysis that reports a boost, and the Newtonian-leaning analyses sit further
below it.

These are point comparisons at each bin's characteristic g_N. The boost is convex in log g_N, so averaging over a bin
raises the low-acceleration prediction further. Caveat: Chae's high point is measured against a multiplicity
calibration at 10⁻⁸ m/s².

A numeric DR4 trend rule is appended to the reading key (`wb_dr4_reading_key_PREREG.md`, 2026-11-25 addendum). Take
R = boost(0.2)/boost(1.0) ×10⁻¹⁰. R ≥ 1.5 has the corner's shape; R ≤ 1.2 kills it. The corner predicts 1.79–1.84.
The corner was never registered, is not a defined model without Milgrom's 2011 rule, and is prior art; this run
cannot move the count or Bucket 0.

## 3. Draco: the dSph leg rests on a modern dispersion now (record)
Explorer 2026-11-24 found that 9.1 ± 1.2 km/s is a compilation value. Walker+2009 print it citing a paper that never
prints it. The modern value is 9.55 −0.25/+0.26 (Geha et al. 2026, arXiv:2602.10200, Table A1).

At M/L = 4 (Wolf estimator), with σ, luminosity, r_h and distance carried:
- the 6.39 cap is excluded at 5.2–11σ. The 11-19 figure was 4.5σ, on the σ error alone.
- after the six-system correction, caps up to ≈ 11–23 are excluded. That is maintainer arithmetic from the explorer's
  error budget.
- 1/Ω_b = 20.3 stays open at 1.5–4.3σ, set by Draco's luminosity.
- no source is cited for the scripts' L = 2.2×10⁵ L☉ in either `cap6p39_vs_dsph.py` or
  `pressure_supported_boost_ceiling.py` (traced 2026-11-25). McConnachie 2012 lists 2.9×10⁵.

The site was updated on six surfaces.

## 4. The review instrument's sensitivity is set by its reading path (frame)
The explorer's 2026-11-24 planted-error test measured the four-persona loop at 31/40 = 0.78 on full page text and 9/40 =
0.23 through fetch summaries. That second path is how the visitor track has usually read the site. Today's visitor
Passes 2–4 had no fetch tool and read full HTML by curl. They found the cross-page numeric drifts (14,610 vs 141;
three different "four C's"; "2× off" on two pages) that earlier passes missed. That is consistent with the
measurement, but it is an anecdote, not a test. **Recommend (visitor-config owner / project lead):** make full-text
reading the visitor track's default, and count numeric catches for two weeks before and after. Arm B (preview deploy)
stays gated. The rate at which the loop *introduces* errors is still unmeasured.

## 5. Note
DR4 is in 7 days. Contacting Chae et al. and Saad & Ting before the data still gates on the project lead, as does the DR4
reading default. The site now prints the number of queued decisions (12) beside the headline count.
