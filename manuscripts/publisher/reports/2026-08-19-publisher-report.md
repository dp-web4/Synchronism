# Publisher Daily Report — 2026-08-19

**Run**: RUN-ID 21856 · **Window**: 24h (2026-08-18 03:42 PDT → 2026-08-19 03:30 PDT)
**Archivist context**: 12th consecutive deliberate zero (0 new numbered Sessions), counts re-derived
at `64c1a0637`; a phantom Session 107 avoided by refusing a filename substring parse.

---

## Phase 0: Publication Recommendations

### New Recommendations
None. The window's material is an executed **null** inside a sector already classified Bucket 3.

### Status Changes

| REC | Was | Now | Reason |
|---|---|---|---|
| **REC-2026-036** Experimental Test Catalog | 0.68 | **0.62** | Defect found *inside* the candidate document, not a missing row — see below |
| **REC-2026-040** α-admixture bound | 0.55 | **0.55 (HELD)** | Sole blocker (external prior-art walk) now 16 days open; priced at one pass on 08-18; not re-escalated |

**REC-036 lowered, and the reason is the pass's main finding.**
`Research/EXPERIMENTAL_TEST_CATALOG.md` carries the S674/S675 census as a **line-1 banner** — *5
executed → all collapsed*, *"of the 9 untested, 0 have a verified first-principles-derived
amplitude"* — and its **body was never updated against it**. The closing *Summary: The Decisive
Tests* still ranked **TEST-04 #1, TEST-14 #2, TEST-15 #4, TEST-07 #5**: three recorded as
executed-and-collapsed in the same file's header, the fourth dimensionally inconsistent (500 Mpc,
m² not m, S632). The fifth (TEST-11) falls under the header's own 0-of-9 no-derived-amplitude
finding. **None of the five survives.** *What "Positive" Would Mean* named TEST-04/07/14/17 —
**4 of 4 disqualified**.

TEST-04 held the #1 decisive slot for **15 months** while Session #107, the session it derives from,
predicts BAO matching ΛCDM at exactly 0% in all five bins; its 10⁻⁴ amplitude has no derivation; its
10⁻⁵ kill threshold is ~3,000× below DESI's best precision. Withdrawn 2026-05-04.

**Fixed in place**: each entry now carries its census verdict inline, and the four names were
**removed rather than reranked** — the catalog holds no surviving candidate to promote. The readiness
change follows from what the fix reveals: the document needs a **regenerated forward list**, not
added rows. Also recorded: the ID space stops at TEST-25, so **both live cosmology registrations are
absent** — TEST-04a (fired, Bucket 2) and TEST-26 (dp-gated).

One ledger disagreement recorded rather than smoothed: the catalog files TEST-07 under *untested*
while `PREDICTIONS.md` Bucket 2 files it under *refuted*. Refuted is better sourced.

### Upcoming Candidates
Unchanged. Preprint / two-paper strategy and the TEST-09/TEST-10 count-collapse question remain
dp-gated.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper
- **Status**: **Updated** (minor, in place) — §5.15 `15-dark-matter`
- **Sessions Reviewed**: through 691 (0 new)
- **Changes Made**: removed the qualifier **"background"** from the γ = ½ degeneracy sentence and
  replaced it with the non-perturbative statement. The qualifier was **not inert**: the site's 08-18
  visitor Pass 4 read the matching *"background-only, no perturbation sector"* caveat as an
  unexplored channel and proposed converting the dp-gated DR3 registration into a 2026 test
  (*"either horn closes TEST-26 in 2026"*). The 08-18 execution refutes exactly that.
- **Verified here before propagation** (`simulations/publisher_20260819_de_locality_identities.py`,
  SymPy, exact, all x and γ, 4/4 pass):
  `ρ_DE/ρ_crit = 2x/((1+x)^{2γ} − 1)`, which is **exactly 2 at γ = ½ with dρ_DE/dρ_m ≡ 0** — constant
  in space as well as time — and `δ_DE/δ_m = 1 + w_DE = ε(ln(1+x)/x − 1) + O(ε²)`, `ε ≡ 2γ − 1`,
  vanishing identically at γ = ½. So the branch is ΛCDM **non-perturbatively**.
- **One correction made against the source, not relayed**: the finding's *"the locality horn improves
  the required precision by 7%"* is a central-value artifact. `σ_γ = (½ − γ̂)/3` is horn-independent
  by construction — it prices a γ-space separation, exactly as §5.15 itself derives it — so 0.004 at
  γ̂ = 0.487 and 0.0037 at 0.489 are one quantity rounded twice; **locality buys 0% there.** The real
  ×2 is in σ(fσ₈): 0.074% → 0.15%. Verdict untouched (~10× underpowered vs DESI DR2's 1–3%/bin;
  ε = −0.022 ± 0.220 = **0.10σ** from exact ΛCDM). Same correction applied to the `PREDICTIONS.md`
  line that had already propagated the 7%.
- **Propagation check**: `background-only` / `perturbation sector` occur in whitepaper sources **only**
  in §5.15 — matches opened, not counted. No exec-summary or conclusion edit owed.
- **Terminology Concerns**: none new.

### Web4 Whitepaper
- **Status**: **Current — 11th consecutive structural zero**
- **Ref scanned**: `origin/main` @ `184dd61`; working-tree HEAD parked on
  `cbp/concepts-normative-home` at an older commit. **A bare HEAD scan would have read 0 commits;
  `origin/main` carries 10.** The 2026-08-09 rule changed the input, not just the wording.
- **Changes Made**: none. Of 10 commits exactly one touches `whitepaper/` — `4c6e6a8`, the web4
  lane's own Publisher run log appended to `PUBLISHER_CONTEXT.md`. Zero normative edits.
- **R7 trigger**: unchanged and correct — **registry publish**, not hub sprint completion.
- **MRH divergence watch (opened 08-18): checked, NOT tripped** — verified, not asserted. The forum
  file is untouched since its creating commit; `git grep` over web4's `whitepaper/` for the contract
  reading returns nothing; Synchronism §4.2 still carries the radius definition.
- **Terminology Concerns**: none.

---

## In-archive corrections executed

1. **`Research/EXPERIMENTAL_TEST_CATALOG.md`** — census verdicts propagated from banner into body.
2. **`Research/Session107_DESI_Forecasts.md`** — erratum added. The primary read *"Status: COMPLETE"*
   with "6.6σ combined" against **three** independent corrections none of which had reached it
   (DR1 refutation 2026-05-05; underpowered-as-registered 2026-07-14; mechanism = the door-#1 C(ρ)
   law applied to growth + a withdrawn `C_cosmic ≠ C_galactic` ratio + an **assumed** σ₈ = 0.76,
   2026-08-18). The genuine DE forecast is −0.22%; the two disagree **173×** and sat in one archive as
   two models never recorded as different — **and the large one is the refuted one.** Forecast tables
   untouched; they stand as the December 2025 record.
3. **`PREDICTIONS.md`** — the 7% artifact corrected in place.
4. **`whitepaper/.../dark_matter.md`** + section CHANGELOG — §5.15 strengthened.

---

## Verification

| Check | Result |
|---|---|
| Symbolic identities (new script) | 4/4 pass |
| `make-md.sh` / `make-web.sh` | exit 0 / exit 0 |
| Churn gate, content (`--ignore-cr-at-eol`) | 98 ins / 500 del |
| Churn gate, raw | **11,984 ins / 12,386 del** |
| Action | artifacts **restored, not staged** — CI is the authoritative builder |

Refutation count **HELD at 6** · **Bucket 0 = 0**.

---

## Standing escalations

- **Site maintainer 401 — day 7, owner action, will not self-heal.** `Failed to authenticate. API
  Error: 401` on 08-16/17/18 (and 08-12/13/14): a dead `CLAUDE_ADMIN_TOKEN`, **not** the account
  outage — the explorer ran a full session on the same host on 08-18. The blocked backlog grew again:
  `/dark-energy` line 195 still publishes *"the sector is background-only. There is no perturbation
  sector"*, which the site's **own** 08-18 finding declares false.
- **REC-040**: 16 days on a gate priced at one pass. Held, not re-escalated.
- **dp-gated**: preprint / two-paper strategy; TEST-09/TEST-10 count collapse; TEST-26 DR3
  registration terms.

---

## Summary

The window's physics was a verified null — the dark-energy sector's last unexplored degree of freedom
(evaluate C locally, as the framework's own postulate requires) adds no independent channel, and at
γ = ½ the sector is ΛCDM non-perturbatively rather than merely degenerate with it. Two corrections
came out of verifying it rather than relaying it: a 7% "improvement" that belonged to a shifted
central value and not to locality, and a whitepaper qualifier that was inviting exactly the inference
the execution refutes. The larger finding is elsewhere: the program's Experimental Test Catalog kept a
withdrawn test in its **#1 decisive slot for fifteen months, underneath its own banner saying that
test had collapsed** — the 08-14 propagation split observed inside a single document, and
independently corroborated the same window by the web4 lane finding the identical shape in its own
specs. **A correction has a landing site, and the banner is not it.**
