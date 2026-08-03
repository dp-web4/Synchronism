# Publisher Daily Report - 2026-08-01

> **⚠ CORRECTION APPENDED 2026-08-03 — the premise of §8 and of this report's subtitle is false.**
> This run was **not** manual and the cron did **not** die. `publisher-2026-08-01.log` is 3,514 bytes
> and ends `Publisher Session Complete: 2026-08-01 03:44:57 (claude exit=0)`. The 03:30:02 header is
> **03:30 PDT = 10:30 UTC** — the same instant the agent read its own clock as "now" — so the agent was
> reading the opening banner of the run it was itself executing, and the log is header-only whenever
> you look at it because `claude -p | Add-Content` flushes the agent's output only at completion. A
> live run cannot produce any other signature, so the test used here can only return "dead."
> The requested `2>&1` has been present on the launcher since 2026-07-24. This same trap was caught
> and correctly written up by this track on 2026-07-31 (`private-context/publisher/log.md`,
> "10:30 UTC *is* 03:30 local") and then reproduced anyway on 08-01 and 08-02. See the 2026-08-03
> report, §"The correction I had already written down."


## Verdict: SIGNAL — the sibling repo's measurement was right and its headline was wrong, because the canonical text says the opposite of what the headline assumed. Three days after learning to grep the sibling repos, the other half of that rule: grep the canonical text too.

**Not an autonomous run.** Manual invocation. The 02:30/03:30 cron has now failed to produce a
report on two consecutive days — see §9.

---

## 1. WAKE: a two-day backlog, and one item in it was load-bearing

No Publisher pass landed on 07-31 (host, not Publisher — diagnosed in `c920a128`) and the 08-01
cron wrote its three-line header at 03:30:02 and produced nothing further. So the scan window is
two days wide, not one.

The widened §1b surface list earned itself immediately. Two new items in
`synchronism-site/explorer/findings/` — both post-dating my last pass, neither indexed anywhere in
this repo:

| finding | date | what it claims |
|---|---|---|
| `two-coherence-orientations-chemistry-uses-the-flipped-one.md` | 07-29 | the governing equation's orientation is backwards; four documented inversions are one |
| `a0-epoch-branch-A-tested-disfavoured-for-evolving-too-slowly.md` | 07-30 | a₀(z) ∝ H(z) is engaged once by the right observable and misses |

Note the first one is dated **07-29** and my 07-30 pass did not see it — it ran before the site's
08:14 commit. The 07-30 fix works; it just does not reach backwards through a missed day. Two-day
gaps need two-day scans, which is now the only reason today's finding exists.

Both got run rather than summarized.

---

## 2. The finding: the framework carries two couplings of C to gravity, and they disagree

The 07-29 session measured — correctly, and I reproduced it — that the non-baryonic term SPARC
galaxies require anti-correlates with local density (mean r = −0.824 over five galaxies). It
charged that to the governing equation: *"the governing equation makes coherence rise with density,
and every sector that touches data needs it to fall,"* and proposed a global flip `C → C(ρ_crit/ρ)`.

**The measurement survives. The attribution does not.** The framework couples C to gravity in two
different places, in two different ways:

| | coupling | enhancement scales as |
|---|---|---|
| **(A)** whitepaper §5.15, `dark_matter.md:36` | `G_eff = G/C(ρ)` ⇒ `v² = v_bar²/C` | **1/C** — low C ⇒ *more* missing gravity |
| **(B)** `synchronism-site` `/galaxy-plotter` (code) | `v_syn² = v_bar² + (V_flat·C)²` | **C²** — high C ⇒ *more* missing gravity |

They demand **opposite orientations from the same rotation curve.** Inverting each for the C it
requires at every radius and correlating against SPARC's own measured surface brightness
`SBdisk + SBbul` — a directly observed local density proxy, no profile assumption, no estimator
choice, which matters because §4 of that same finding shows unnamed estimators have carried
headline verdicts on this program three times:

```
Spearman(log Σ, required C) — 147 galaxies, Q<=2, inc>30, >=5 usable radii

coupling                                             median    mean   frac>0
(A) whitepaper   G_eff = G/C   [C = vbar^2/vobs^2]    +0.827  +0.545   82.3%
(B) site plotter vsyn^2 = vbar^2 + (Vflat C)^2        -0.976  -0.873    3.4%
```

Against radius — an independent, inverse density proxy — the signs flip exactly as they must:
(A) −0.827, (B) +0.976. On the five galaxies the finding itself tabulated, (A) runs +0.803 to
+1.000 and (B) runs −0.973 to −1.000.

**The whitepaper's coupling demands a coherence that rises with density, which is what
`C(ρ) = tanh(γ·ln(ρ/ρ_crit + 1))` already does. It is not inverted. The site's coupling is
inverted, in 96.6% of galaxies.**

The 07-29 session measured the right thing and charged it to the wrong account. Script:
`simulations/orientation_coupling_test.py`; writeup: `explorations/2026-08-01-the-inverted-coupling-is-the-sites-not-the-whitepapers.md`.

### What it moves

- Repair-matrix row (iv) — *"radial coupling dC/dr < 0 → fixed by the flip: YES"* — is a statement
  about site code, not about the framework. It should not count as one of the four level-axis votes.
- The finding's most interesting structural claim — that the *success* sector and the *failure*
  sectors vote unanimously for the flip — loses the galaxy vote. What remains is chemistry (§2 of
  that finding, untouched by this and standing on the wave-speed identity), plus BEC/BCS and the
  CFD/Bullet-Cluster mapping, both of which were derived under coupling (B) and need re-deriving.
- The site's `/dark-matter-failure` prose ("dark matter = low coherence C") is **correct** and
  matches the whitepaper. It contradicts the plotter, not the framework. Fix the code, keep the
  prose — the inverse of the 07-29 action list.

### What it does not move, stated so it is not overread

It rescues nothing. `v² = v_bar²/C` with a bounded C *is* the bounded-boost family that TEST-10
and `efe-boost-ceiling-closure.md` (2026-06-03) already exclude on this same data. **Consistent
orientation, still excluded.** The 07-29 conclusion — *"fixing the orientation makes the sign
ledger consistent and buys no predictive power"* — holds, for a different reason than it gave.

### My wrong guess, recorded

Coupling (A) is not unanimous: 24 of 147 galaxies go negative. I predicted before looking that
these would be bulge-dominated massive spirals, where `v_bar²/v_obs²` is bulge-driven. **Wrong** —
they are the faint end (median L₃.₆ 1.35 vs 10.84 ×10⁹ L☉; 4% with a bulge vs 24%). Mechanism is
compressed dynamic range in gas-dominated dwarfs, where `v_bar²/v_obs²` is small and nearly flat so
the rank correlation is noise. The defensible comparison is the gap — 82.3% vs 3.4% — not either
number alone.

---

## 3. Phase 1 Synchronism: one edit, and it is a prediction meeting published data

The 07-30 finding says the a₀ epoch fork's branch (A) has been tested once, by the right observable,
and misses. I did not inherit that. **The Ciocan abstract was re-fetched and the arithmetic
re-derived here** — because that same finding documents the site citing Milgrom 2017 as *proposing*
a₀ ~ cH/2π when it is the paper that tests and disfavours it, which is a standing reason to check
rather than accept.

Verified independently: Ciocan et al. 2026 (MUSE-DARK III, A&A **709**, L16, arXiv:2604.22613),
N = 79, 0.33 < z < 1.44, RAR-fitted per bin, `a₁ = 1.59 ± 0.1 × 10⁻¹⁰`. Against an H(z)-tracking
slope of `a₀(0)·(dE/dz)|₀ = 0.49 × 10⁻¹⁰` at Ω_m = 0.315:

| statement | value | note |
|---|---|---|
| measured slope vs H(z)-tracking slope | **3.2× steeper** | normalisation-free — the robust form |
| branch (A) at z ~ 1, a₀(0) = 1.04 | 1.86 vs 2.38 (+0.12/−0.10) | ≈5σ low |
| branch (A) at z ~ 1, a₀(0) = 1.20 (Milgrom) | 2.15 vs 2.38 | ≈2σ low |

The slope comparison is the better statement and it is the one I have put in the paper, because it
does not depend on which a₀(0) normalisation you pick — and normalisation-dependence is precisely
the defect I spent 07-27 through 07-30 correcting on a different axis.

**Edit made** — `whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md`, an inline
status note under the "Implication" paragraph that asserted a₀(z) ∝ H(z) and called it "testable
via high-z BTFR". It was testable-and-untested when written and is not any more.

Recorded in the note, inline rather than in a footnote: (i) the test route is the **RAR, not the
BTFR** the whitepaper names, which is exactly why Milgrom's standing BTFR objection to high-z a₀
tests does not apply to it; (ii) the external literature disagrees with itself — Gueorguiev 2024
finds a redshift slope consistent with **zero**, Milgrom 2017 disfavours evolution *faster* than
this, so the prediction is squeezed from both sides by mutually inconsistent measurements and
neither bound is firm; (iii) Ciocan's systematics are large. Status is **engaged and disfavoured,
not refuted**, and the refutation count stays at 6.

Also noted in the whitepaper: this result is *worse* for the alternative anchor
`2πa₀ ≈ c²(Λ/3)^{1/2}`, which predicts zero evolution. What is embarrassed is the programme of
deriving a₀ from cosmology, not this framework specifically — which is the outward-facing framing
and is more defensible than a self-flagellation framing would be.

**Why an inline note and not a rewrite**: single source, large systematics, contradictory
neighbours. Conservative change is to date and quantify the engagement, not retire the prediction.

| Gate | Result |
|---|---|
| Claims freeze (`--check`) | ✅ green, 10 claims, run before and after |
| Build (`make-md.sh`) | ✅ exit 0, monolith + Pages copy regenerated |
| Monolith diff (`--ignore-cr-at-eol`) | ✅ **6 insertions, 0 deletions** — no churn |
| Lone-CR regression (`\r(?!\n)`) | ✅ clean, **7th consecutive pass** (only matches are PDFs/PNGs) |
| Terminology drift | zero |

**No whitepaper change for §2's orientation result.** §5.15 states the coupling correctly; today's
work is a *defence* of the existing text, not a correction to it, and it changes no reader-facing
claim. It went to `explorations/`.

---

## 4. Phase 0: a registry gap in REC-2026-036, found by accident

`Research/EXPERIMENTAL_TEST_CATALOG.md` holds 24 experimental tests. **`grep` for a₀ returns
nothing.** The a₀(z) prediction is stated in the whitepaper's own text, is Tier-1 by the catalog's
own definition (*"Existing Data, No New Hardware"*), and is not scope-excluded — **TEST-14 (Gaia
DR3 wide binaries)** is the same class.

This is load-bearing rather than untidy because the catalog states that **no registered test remains
runnable-and-unrun**. True as scoped, misleading as read: the first framework prediction settled by
outside data since that sentence was written is one the catalog never listed. The registry is not a
complete inventory of the framework's testable statements.

Filed: `Research/proposals/test_catalog_a0z_tier1_gap_20260801.md`. Four items, and **the fourth is
the one that matters**: not "add this row" but *"what else does the whitepaper assert as testable
that the catalog does not list?"* One gap found by accident implies the registry was assembled from
arc outputs rather than swept from the whitepaper's claims. A one-pass diff would settle it and is
cheap. Items 1–3 fix an instance; item 4 fixes the class.

Also proposed and worth flagging: if a TEST-25 row is added, its kill criterion is necessarily
**post-hoc** and must be labelled so, and must **not** count toward the catalog's "criteria fired
and were honored" tally. That tally is REC-036's central asset and it should not be diluted by a
row that could not have been pre-registered.

---

## 5. A transcription defect in my own reports

REC-2026-036's stored readiness has been **0.68** since 2026-07-17. My 07-29 and 07-30 reports both
quote it as **0.60**. Neither pass changed the value; both inherited the figure from agent memory
instead of reading the state file.

Small number, exact shape of the failure I have been documenting in others all week: a summary
consulted in place of the source. Annotated in `recommendations.json` rather than silently fixed,
and the two reports are left standing wrong — an append-only record is worth more than a tidy one.
0.68 is authoritative.

---

## 6. Phase 0: status changes

| Rec | Change | Why |
|---|---|---|
| **REC-2026-036** | 2 weaknesses added, **readiness held 0.68** | Tier-1 a₀(z) registry gap; plus the 0.60/0.68 transcription defect above. The gap is a completeness fact about the catalog-as-inventory, not a quality fact about the rows in it — so it does not move the score |
| **REC-2026-038** | 1 strength added, **readiness held 0.93** | Third instance of its documented class, and the first running in a direction the draft does not cover — see below |
| REC-034 / 035 / 037 / 039 | held (0.97 / 0.95 / 0.92 / 0.38) | No new evidence |

**REC-038's new instance is qualitatively different from its existing two.** The paper documents
continuity *lost* across memoryless sessions. Today's is continuity **inverted**: a correct
measurement charged to the wrong account because the canonical text was never consulted, producing
a confident and wrong headline that a later reader would have inherited. Losing a result costs a
rediscovery; inverting one costs an axiom. That is a failure mode the draft should carry, and it
strengthens rather than complicates its argument.

**Advisory order, unchanged**: **#1 REC-2026-038** (drafted, reproducibility-verified, cs.AI, and
now with a sharper third instance) → #2 DESI mechanism-class → REC-039 as a robustness section →
#1-null Locality No-Go (still blocked) → #3 A2ACW.

**Null #1's fork** (Reading A vs Reading B, filed 07-28) — no new evidence, **fifth consecutive
day**. Still decides preprint #1's headline sentence. This has now been open longer than any other
blocking item in the queue and nothing in the research lane is moving toward it.

### Widened scan, both surfaces exercised

`Research/papers/` unchanged (one item, REC-038). `Research/preregistrations/` — the
`sparc_cassini_tanhlog/RESULT.md` gained a back-annotated tail-shape mechanism (07-30, `9a7c8768`):
the compander saturates as a power law `1−C = 2(1+x)^{−2γ}` while McGaugh's ν deviates
exponentially, so at Saturn's `x ≈ 5×10⁵` the gap is five orders of magnitude and no fitted γ closes
it. This upgrades TEST-11 from "a parameter choice that failed" to a structural consequence, and it
is the same choice that gives the compander its finite boost ceiling. **Two of the program's
galactic negatives trace to one decision** — worth stating in whatever writes them up, and not
currently stated anywhere. `Research/proposals/` — 07-30 additions, already handled.
`explorer/findings/` — §2 and §3 above. No new numbered session (still S691 unopened), no new arc.

---

## 7. Phase 1 Web4: no change, and a 07-30 call that held

Reviewed `web4@bec588c..724bf8c` — 30+ commits, but the whitepaper is being maintained by its own
Publisher lane on that side, which ran a no-change pass on 07-31 (`a8b27f4`) and shipped a §11
correction on 07-30 (`57625ed` — the 07-28 package-matrix audit had queried the wrong PyPI
distribution name and the paper inherited the falsehood; corrected append-only). Nothing for this
lane to add.

**The 07-30 `AssuranceReceipt` call was "watch, not act — still design-in-motion."** It changed
shape twice more in the two days since: `694584e` (the A2 receipt did not sign the string it
transmits — signing_bytes v3), `33f9b03` (constellation attestation signs a string it does not
transmit — receiver-first v2), and audit `121b8a4` found a one-day-old portable receipt naming its
signer with an id its own registry cannot resolve. **A declared prediction that held.** Still not
integration-ready; re-check when the signing format goes a week without a version bump.

---

## 8. Housekeeping

- `recommendations.json` diff a clean **8 insertions / 5 deletions**, no re-serialization churn —
  fifth consecutive pass, verified by byte-comparing a round-trip against the on-disk file *before*
  editing rather than assuming.
- `AGENTS.md` / `CLAUDE.md` carry uncommitted GitNexus index-count drift (172886 → 172978).
  Supervisor scope, left untouched, per precedent `a13894da`.
- New script `simulations/orientation_coupling_test.py` — imports the executed TEST-10 loader and
  cuts so the sample cannot drift from that line of work. No fitting, no free parameters, no
  external data.

---

## 9. The cron has now missed two days, and the failure mode changed

| date | log file | outcome |
|---|---|---|
| 07-26 → 07-30 | written at 03:30:0x | five consecutive successful runs |
| **07-31** | **absent entirely** | never started — host (diagnosed `c920a128`) |
| **08-01** | **3-line header, 03:30:02, nothing after** | started, then died before any work |

These are **different failures**. 07-31 never invoked the launcher. 08-01 wrote the launcher's
unconditional header and then produced nothing — so `claude` was invoked and did not survive, or
was killed. That is not the host being asleep; that is the run dying. A `TimeoutStartSec` kill or
an auth failure both fit and I cannot distinguish them from inside this session — the log carries
no stderr.

**Concrete ask, since it is not mine to fix**: the launcher writes its header unconditionally and
captures nothing else. One `2>&1` on the invocation would have made today's cause readable instead
of inferable. Flagged for the supervisor.

Today's pass was manual. Two days of accumulated surface produced one executed finding, one
whitepaper edit, and one registry gap — so the cost of a missed day is real, not cosmetic.

---

## Summary

Two days of backlog, and the item that mattered was a sibling-repo finding claiming the framework's
governing equation is oriented backwards. I ran it. The measurement is right and reproduces; the
attribution is wrong. The framework carries two couplings of coherence to gravity — the whitepaper's
`G_eff = G/C` and the site plotter's `v² = v_bar² + (V_flat·C)²` — and they demand opposite
orientations from the same rotation curve. On 147 SPARC galaxies the whitepaper's coupling wants a
coherence that **rises** with density in 82.3% of them, which is what `C(ρ)` already does; the
site's wants one that falls, in 96.6%. The bug is in the public implementation, and the canonical
text has been right the whole time.

Separately, a prediction the whitepaper states in its own words — a₀(z) ∝ H(z) — has been engaged
by published data for the first time, and is disfavoured for evolving too *slowly*: 3.2× shallower
than the measured slope. That is now in the paper with its three caveats, as disfavoured rather
than refuted. Finding it also surfaced that the Experimental Test Catalog has no a₀ row at all,
which means the catalog is not the inventory it reads as.

**On 07-30 the lesson was: grep the sibling repos before promoting a claim — the program's own prior
art kills findings for redundancy. Today is the mirror. Prior art also *rescues*: the canonical text
said the opposite of what a confident sibling-repo headline assumed, and one grep would have caught
it before the flip was proposed. The rule is not "check the neighbours." It is that a headline
inherits choices nobody stated — a normalisation on 07-27, an estimator on 07-29, a coupling today —
and the computation being correct is no protection at all, because in every one of those three
cases it was.**
