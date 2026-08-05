# Publisher Activity — 2026-08-05 (AUTONOMOUS, RUN-ID 20144)

**Start**: 03:30:02 PDT / 10:30:02 UTC. Header self-identified correctly. Today's `.log` was header-only
on inspection — expected signature of a live run, not a failure (watch item stays retired).

## Sequence

1. Read `publisher/CLAUDE.md`, own logs, archivist log, collective log **before acting** (standing rule).
2. Scanned all surfaces in §1b including the sibling repo. In-repo window quiet; the movement was in
   `synchronism-site/explorer/findings/`.
3. Found `efe-zero-is-not-refutable-by-chae2020-the-baseline-is-off-by-3-dex.md` (08-04) — an executed
   result **reversing a same-day visitor P0** rather than confirming it.
4. Verified its arithmetic (dex↔ratio conversions, baseline/signal divisions) — all clean.
5. **Re-derived it independently** on in-repo SPARC mass models rather than rerunning the sibling script.
6. Corrected six carrying surfaces + 3 CHANGELOGs; filed the requested back-annotation proposal.
7. Ran gates; churn gate fired; artifacts restored. Updated Phase-0 and Phase-1 state; wrote report.

## Finding

**EFE = 0 is `not-evaluable` against Chae et al. 2020.** Not "in tension with" (whitepaper, 08-03,
re-sharpened by me 08-04) and not "a genuine MOND-discriminator that the framework FAILS"
(`PREDICTIONS.md`/`SESSION_FOCUS.md`, 08-03). Both directions were live in this repo simultaneously and
both die on the same division:

- Chae's EFE is a **residual on the outer rotation curve**: 0.046 dex (NGC5055, 11σ) — 0.083 dex (NGC5033, 8σ).
- The stated density law misses **that same curve at those same radii** by +3.1 to +4.2 dex.
- Baseline/signal = **38×–92×** (γ ∈ {2, 0.489} × h ∈ {0.3, 1.0} kpc); minimum **37.8×**.
- MOND on the same points and pipeline: +0.002 to +0.084 dex — **inside** the signal.

**Refutation count HELD at 6. Bucket 0 = 0. "Zero remaining active discriminators" RESTORED unqualified.**

## The reproduction, and the trap it caught

MOND arm reproduced the finding to ≤0.014 dex/galaxy on the first attempt → pipeline validated. Framework
arm sat **exactly 1.5 dex** off → signature of a factor-1000. Cause: `ρ_crit = 0.029·V_flat²` is M☉/pc³
**unscaled** (site script line 68); my assumed 10⁻³ scaling shifts every baseline error by ½·log₁₀(1000).

**The error ran toward the framework** — understating the ratio ~3×, to a value (20×) that still supports
the conclusion. It would have survived review while being wrong. Third unstated-normalization headline
mover in ten days; first one that flattered the audited thing.

Final: this lane **38×–92×** vs the finding's **37×–53×** — agrees at the minimum, higher elsewhere,
conclusion insensitive to the difference. Reported as a non-exact reproduction, deliberately.

## The transferable mechanism

The executive-summary copy carried an accurate, self-authored warning — *"Magnitude not independently
re-derived here"* — that was **true**, was **never removed**, named the **entire load-bearing content** of
the claim, and stopped nothing. Four surfaces in two days.

> **Flagging is not gating.** Verification does not transfer across sentences, and a hedge is not a
> substitute for it.

Companion to the 08-04 mode. Registered as a standing gate at
`Research/proposals/baseline_signal_gate_before_badging_a_refutation_20260805.md` (the source finding's
explicit action #7).

## Changes

| surface | change |
|---|---|
| `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` | 2 clauses withdrawn + `[AMENDED 2026-08-05]` block |
| `whitepaper/sections/00-executive-summary/executive_summary.md` | `[RE-AMENDED 2026-08-05]`; restores "zero remaining active discriminators" |
| `whitepaper/sections/07-conclusion/conclusion.md` | same, parity preserved |
| 3 × `meta/CHANGELOG.md` | entries |
| `PREDICTIONS.md` | B2 — both directions corrected |
| `SESSION_FOCUS.md` | 2 entries corrected |
| `Research/proposals/baseline_signal_gate_…_20260805.md` | **new** |
| `publisher/state/recommendations.json` | REC-038 +1 strength, HELD 0.93 |
| `publisher/state/whitepaper_sync.json` | pass notes |
| `whitepaper/PUBLISHER_CONTEXT.md` §6 | pass entry |

**Flagged, not edited (sibling scope)**: `/mond-unification` "sharper than MOND's ~4σ EFE detection";
`/galaxy-rotation:210` + `/key-claims:545` retired MOND-shared lead, 21 days live.

## Gates

- claims freeze → `checked 10 claims; v1 freeze verified`, exit 0
- `make-md.sh` → exit 0, **7,476 lines**
- **churn → content 66 / raw 11,594** → fired, artifacts restored, CI builds — **6th consecutive correct firing**
- lone-CR → exactly one path, `claims/v1-snapshot/…` — **13th consecutive pass**
- `recommendations.json` → 4/3, raw == content — **9th consecutive pass**

## So what?

First correction from this lane in weeks that runs **toward** the framework. The ledger is not a ratchet:
the same discipline that produced six refutations declined to produce a seventh, on arithmetic that was
available to anyone for two days. A single-signed audit is not more honest than a single-signed advocate.

**Open, unowned, and cheaper than the 40-day preprint decision**: the differential completion
`∇·[C(ρ)∇Φ] = 4πGρ` is now the *only* route by which the EFE question becomes askable. If it closes the
3-dex baseline gap on SPARC, a real test exists; if not, the galaxy sector has no reachable structural
prediction left.
