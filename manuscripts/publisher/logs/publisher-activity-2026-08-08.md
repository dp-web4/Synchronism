# Publisher Activity Log — 2026-08-08

**Run**: autonomous daily pass, Synchronism/manuscripts. **Closing banner at the bottom** — if it is absent,
the run did not finish (see `read-my-own-log-before-acting`: a live log is header-only because
`claude -p | Add-Content` flushes at completion, so absence of the banner is the only valid liveness read).

---

## WAKE — am I working on the right thing?

Archivist (09:30 UTC) reports **0 new Synchronism core Sessions, 4th consecutive deliberate zero**, and
7 Synchronism commits in the window: a kimi saturation-horizon synthesis waking the hive-organs arc, and a
door-#1/EFE line that registered a blocker and then **refuted its own blocker** within one day. Three of the
seven commits are authors removing their own claims.

So there is no arc to evaluate. The live question is whether any of the four new documents touches a standing
whitepaper claim. Two do not (the coarse-graining supersession converges with what this lane derived
independently on 08-07; the EFE proposal explicitly holds the count and reinforces the 08-05 withdrawal). One
does, on an axis I had not examined: it says the count of 6 is **convention-dependent**, and that
**TEST-09 and TEST-10 are one structural root**.

My own log has written "Count HELD at 6" on five consecutive days without ever stating what convention the 6
counts. That is `headlines-inherit-unstated-choices`, sixth instance, and this time the unstated choice is
mine and it is in the whitepaper. **That is the right thing to work on.**

---

## Scan (§1b, all seven surfaces — the list exists because this scan has been blind four times)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new core Sessions |
| `Research/papers/` | unchanged since 2026-07-23 (REC-038's manuscript, now 16 days untouched) |
| `Research/proposals/` | **3 new (08-07)** — all read in full |
| `Research/preregistrations/` | unchanged since 2026-07-30 |
| `explorations/` | 1 new (08-07 kimi saturation horizon) — arc-charter scope, not publication scope |
| `synchronism-site/explorer/findings/` | newest 08-07 08:13, read via the proposal that carries it |
| `synchronism-site/explorer/topics/` | nothing new; `done/` touched 08-07 |

Also read, and load-bearing today, both **outside** the §1b list: `synchronism-site/maintainer/logs/2026-08-07.md`
and `synchronism-site/visitor/logs/2026-08-07.md`. Those two are where the count figures actually originate;
the proposal only reports them. **Candidate §1b addition** — the sibling repo's *lane logs*, not just its
`findings/`, because a claim can circulate in a log for a day before any finding or proposal carries it, and
today both circulating numbers were log-only in origin.

---

## FOCUS — what was done

### 1. Verified the incoming claim rather than trusting it

The proposal's count argument traces to two primary sources, both of which I opened:

- `explorations/2026-07-15-verify-test08-test10-executed-mond-shared-class-law.md` — *"the bounded boost
  B ≤ 1/Ω_m = 3.17 caps f_DM = 1 − C at 1 − Ω_m = 68.5% … same structural ceiling as TEST-09, different
  observable."* On file **since 07-15**.
- `Research/proposals/boost_ceiling_provenance_and_class_exclusion.md` **open question 4**, added **07-28**:
  *"Should TEST-09 and TEST-10 be counted as one refutation or two? … a real methodological question (how
  much independent empirical content does a corollary carry?), not a typo to silently fix … Gates on dp."*

So the observation is 24 days old and the question is 11 days old and **registered**. It was nonetheless
re-filed as a novel P0 on 08-07 by **two independent non-author readers** within hours of each other. That
reframes the item: the failure is not retrieval, it is that the register does not reach the page doing the
counting.

I also checked the *inverse attribution*: the visitor claimed TEST-09/10 fire only under the amplitude force
law; the maintainer corrected this to the division law and I confirmed the correction at the code
(`test09_btfr_bounded_boost_real_sparc.py` codes `g_obs = g_bar/C`). Right defect, inverted attribution —
already fixed upstream, not repeated here.

### 2. Measured the whitepaper's own exposure

Four assertions of "unchanged at 6" (`00-executive-summary`, `07-conclusion`, `15-dark-matter` ×2 in body
prose plus the CHANGELOGs), **zero enumerations**. Then, by grep over `whitepaper/sections/` excluding
`meta/`:

```
TEST-10 → 0    3.17 → 0    68.5 → 0    0.927 → 0    "boost ceiling" → 0
CHSH    → 0    no-signaling → 0    Kuramoto → 0    1.98 → 0
TEST-09 → 14   TEST-11 → 3    Cassini → 6    17.95 → 1    0.0001 → 10    BTFR → 20
```

**Two of the six counted refutations are described nowhere in the whitepaper.** The `Bell` hits are §5.4's
exposition of Bell's theorem and a prior-art specificity list — checked in context, not counted blind. The
`dwarf` hits are generic ("dwarf galaxies where the effect is strongest") — also checked.

### 3. Executed the measurement (`simulations/publisher_20260808_test09_test10_independence.py`)

Reproduced both site executions first (TEST-09 returns 3.35 vs observed 3.75 catalogue / 3.79 estimator —
the ledger's 0.41 is against the catalogue V_flat; my 0.44 is against the like-for-like outer-3 estimator the
script's own docstring calls the comparable one; **both fire**, no discrepancy to report). Then three steps:

1. **Identity.** `f_DM = 1 − 1/B` ⇒ TEST-10's criterion and TEST-09's ceiling are the *same inequality*.
   `max|f_DM − (1 − 1/B)| = 2.2×10⁻¹⁶` over 123 galaxies; the two criteria select **93 and 93** galaxies with
   **0 disagreements** at the outermost point and **88 and 88** with 0 disagreements at the error-weighted
   outer-3 point. Not "corollaries of a shared root" — one inequality in two variables.
2. **Shared inputs.** Same SPARC files, same Υ = 0.5/0.7, same force law, same Ω_m = 0.315, φ = 1.618034,
   a₀ = 1.05×10⁻¹⁰. Verified by opening the scripts.
3. **Conditional kill.** Delete every TEST-10 firing galaxy, refit the BTFR. Full sample: 3.79 ± 0.10 vs
   3.35 ± 0.07 ⇒ 0.44 (3.5σ), fires. Exceeders only (N=93): 0.33 (2.6σ), fires. **Non-exceeders (N=30):
   3.39 ± 0.23 vs 3.38 ± 0.18 ⇒ 0.01 (0.0σ), does not fire.**

**And then I stated the power of my own null**, because that is the 08-07 lesson: σ on the non-exceeder
deviation is **0.29**, so the 0.3 bar sits **1.0σ** out, and the mass lever arm shortens 3.45 → 3.06 dex. The
subsample *could not* have shown independence had it existed. Verdict written as **"burden unmet, not
disproved"** — the strongest honest form.

**Both circulating figures are wrong in the same direction.** "3–4, not 6" (maintainer) and "at most four"
(visitor) both over-collapse; on the measured axis it is **6 → 5**, one pair. Note the visitor's arithmetic
independently: if two of six share one root, six becomes five, not four.

### 4. Found the axis nobody was counting

Five of six rows run on **SPARC**; three of the scripts behind them share a **fixed** Υ = 0.5/0.7; every
non-SPARC channel enters only jointly or by reference (Cosmicflows-4, Cassini, the textbook CHSH violation);
and the sixth (B1) analyses **no external dataset at all** — it runs the framework's own harness. So the
site's "6 refutations executed on external data … astronomical, ephemeris, and laboratory" is loose on its
third type, and **an M/L shift moves five rows together**. This program has already been bitten once on that
exact axis (REC-039's "28 galaxies" is 15 at Υ*_disk = 0.7). Root-counting is the less useful question.

### 5. Integrated, conservatively

Five-paragraph amendment to `15-dark-matter`'s 08-03 note + CHANGELOG entry + §6a entry in
`whitepaper/PUBLISHER_CONTEXT.md`. **Count HELD at 6, no row merged, Bucket 0 = 0** — the naming decision is
dp's, registered 07-28, and I did not pre-empt it. Re-read the new text against itself afterwards (the second
half of `phase0-scan-is-shape-blind`) and caught three of my own errors before commit: I had written "only two
channels independent of SPARC" while forgetting Cosmicflows-4 (three, and each joint); I had written a
sentence implying B1 rests on SPARC (it rests on nothing external); and I had implied 93/123 = 76% conflicts
with TEST-10's published 69% when it is the same criterion on a narrower cut.

### 6. Web4 — no change, verified not asserted

3 commits, all `docs/audits/` only, zero whitepaper scope, zero `web4-standard/` mutation. PAIRED-CHANNELS
fired by judgement on a **new** channel — C336's retired crypto tokens — and returned no action:
`P-256|ECDH|AES-(128|256)-GCM|SHA-3|W4-(BASE|FIPS|IOT)-1|cose:` over `whitepaper/sections/` = **0 hits**; the
only matches in the tree are in the frozen `archive/sections-2026-07-09-pre-rewrite/`. Fourth independent
confirmation of the structural-zero prior, first from the crypto-suite direction. Last whitepaper-source
commit remains `57625ed` (07-30).

---

## Gates

| Gate | Result |
|---|---|
| `make-md.sh` | PASS (7,529 lines, 680K) |
| `claims/render_claim_surfaces.py --check` | exit 0 — 10 claims, v1 freeze verified |
| Edit present in build | `count of test IDs` → 1/1 (build + docs monoliths) |
| Artifact churn, **both** numbers | content 30 / raw **23,412** ⇒ churn. `git checkout` restored both trees; CI builds. **9th firing** |
| Lone-CR | binary paths only (PDF/PNG) — no text path, no WSL-mount blame |
| `recommendations.json` | 14 insertions / 9 deletions, valid JSON, `ensure_ascii=True` + trailing newline |

---

## So what?

The count was never wrong. What was wrong is that it was **unwalkable** — four assertions, no list, two
members undescribed — and an unwalkable claim gets re-discovered rather than remembered. That is the same
mechanism the site's own lane named on 08-07 for the unregistered EFE = 0 row, on the same day it
instantiated it. New rule: **a bare tally is a claim about independence; publish the list or drop the
number.**

Two things I would not have predicted. First, the direction: the honest-sounding downgrade circulating in the
sibling repo ("3–4 roots, not 6") **overstates** the collapse, so today's correction runs *toward* the
framework on the root axis while running *against* it on a new axis (five of six rows share one dataset and
one M/L). Third correction in five days running toward the framework, and the first that does both at once.
Second, the base rate: REC-038's standing worry was that its instance count is an artifact of one observer
counting where they look hardest. Today's instance is **externally sourced, doubled, and against an item
already adjudicated open** — which is the first evidence that the re-discovery rate is a property of
unenumerated claims rather than of my attention. That belongs in the manuscript, which is still untouched
since 07-23.

**Candidate §1b addition for the next pass**: the sibling repo's `maintainer/logs/` and `visitor/logs/`.
Both of today's load-bearing numbers originated there and were only *reported* by the proposal I would
otherwise have read alone.

---

**CLOSING BANNER — run complete 2026-08-08.** 1 whitepaper integration (Synchronism §5.15), 1 verified
no-change (Web4), 3 recommendations updated (038, 039, 036 — all readiness HELD), 0 new recommendations,
1 new script. Count HELD at 6. Bucket 0 = 0.
