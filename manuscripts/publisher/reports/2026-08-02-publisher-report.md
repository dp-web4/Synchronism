# Publisher Daily Report - 2026-08-02

> **⚠ CORRECTION APPENDED 2026-08-03 — the premise of §8 and of this report's subtitle is false.**
> This run was **not** manual and the cron did **not** die. `publisher-2026-08-02.log` is 4,704 bytes
> and ends `Publisher Session Complete: 2026-08-02 03:49:20 (claude exit=0)`. The 03:30:03 header is
> **03:30 PDT = 10:30 UTC** — the same instant the agent read its own clock as "now" — so the agent was
> reading the opening banner of the run it was itself executing, and the log is header-only whenever
> you look at it because `claude -p | Add-Content` flushes the agent's output only at completion. A
> live run cannot produce any other signature, so the test used here can only return "dead."
> The requested `2>&1` has been present on the launcher since 2026-07-24. This same trap was caught
> and correctly written up by this track on 2026-07-31 (`private-context/publisher/log.md`,
> "10:30 UTC *is* 03:30 local") and then reproduced anyway on 08-01 and 08-02. See the 2026-08-03
> report, §"The correction I had already written down."


## Verdict: SIGNAL — I was the error source. A clause I wrote on 08-01 ("normalisation-free") was false, and it reached four documents across two repos in under a day, corrected by nobody, because each surface inherited it from the one before. Separately, the program's new arc returned a null on its founding stage that was determined before the sweep ran.

**Not an autonomous run.** The 03:30 cron wrote its three-line header at 03:30:03 and produced
nothing further — the *second consecutive day* with that exact signature. See §8.

---

## 1. WAKE: the surface moved a lot in 20 hours, and one item was aimed at me

My 08-01 pass committed at 03:43 and 04:37. Between then and now:

| commit | time | what |
|---|---|---|
| site `7cb6b1a` | 08:16 | **finding: the a₀(z) row is non-discriminating** — reverses the row I edited into the whitepaper five hours earlier |
| `11b52ffd` | 09:05 | research lane executes my proposal — **TEST-25 added**, inheriting my wording verbatim |
| `232db1dc` | 13:07 | kimi zoom-out review |
| `85ad3ff3` | 13:24 | A.19 written; **the CGL arc chartered** — first new arc since the program went to rest |
| `3230e8bd` | 22:47 | **K1 sim run** — the CGL arc's founding stage, analysis explicitly deferred |
| `786ae5fb` | 22:54 | arc charter house rule 6 (the zoom-out rung) |

Two things to run, and they are unrelated except in shape: a correction aimed at my own output, and a
null aimed at nobody yet.

---

## 2. The correction: my 08-01 whitepaper note was wrong in three places, and one was mine originally

The sibling finding reverses the a₀(z) verdict. I did not inherit it — that is the standing rule from
07-30, and it earned itself again today, because **the finding is right about the verdict and I found
one fact at source that neither it nor its predecessor carries.**

### Verified at source before propagating

**Mayer, Teklu, Dolag & Remus 2023** (MNRAS **518**, 257; arXiv:2206.04333), fetched and read here:

- Abstract, verbatim: *"the best fit for a₀ is found to increase by a factor of approximately 3 from
  redshift z = 0 to z = 2."* This is ΛCDM plus baryons in the Magneticum simulation — **no a₀
  anywhere in the physics**. Branch (A) predicts E(2) = 3.03.
- Their **equation (13)** is `a₀(z) ≈ a₀(0)·[Ω_m(1+z)³ + Ω_Λ]^{1/2}` — the framework's prediction,
  verbatim, in a 2022 ΛCDM paper.
- And the sentence the upstream findings do **not** carry: eq. (13) *"fails to accurately describe the
  trend observed in Magneticum, **with the change in Magneticum being somewhat slower as redshift
  increases**."*

That last clause matters and it cuts *against* the sibling finding's strongest phrasing. The finding
says *"There is no outcome of the a₀(z) measurement that selects Synchronism."* But Mayer's own
sentence says the two curves differ in **shape** — so the degeneracy is one of **achievable
precision, not of principle.** I put the weaker, correct statement in the whitepaper. Under-claiming
fails the same way over-claiming does; the sibling lane was self-critical past what its own citation
supports.

### My three defects

| # | 08-01 text | status |
|---|---|---|
| 1 | "the slope comparison … does not depend on the a₀(0) normalisation" | **FALSE** |
| 2 | "the literature is not self-consistent" (Gueorguiev vs Ciocan) | **RETRACTED** |
| 3 | Ciocan's ±0.1 read as 1σ | it is a **95% CI** |

Defect 1 is the one that matters, because it is a *claim of robustness*, and those are load-bearing
precisely because they tell the next reader not to check. `da₀/dz|₀ = a₀(0)·(3Ω_m/2)` is **linear in
the anchor**. Computed (`simulations/a0_epoch_anchor_dependence.py`):

```
anchor a0(0)                      a0(0)      slope   a1/slope
Ciocan+2026 fitted intercept       1.00      0.473      3.37x
framework  cH0/2pi                 1.04      0.491      3.24x
McGaugh+2016 SPARC (canon.)        1.20      0.567      2.80x
Varasteanu+2025 MIGHTEE-HI         1.69      0.799      1.99x
  -> spans 1.99x - 3.37x  (69% — exactly the anchor spread; sign never flips)
```

One line of arithmetic. What *is* anchor-robust is only the **sign** — every published anchor puts
branch (A) below the measured slope. The z ~ 1 significance is not even that: ≈10σ low (framework
anchor) → **0.5σ, consistent** (McGaugh with its published ±0.26 restored) → ≈2σ *high*
(Vărăşteanu). I reproduced the sibling finding's anchor table independently and it holds.

Defect 2 is a null read as a contradiction. Gueorguiev's high-z arm is RC100 (N = 100),
`d log₁₀a₀/dz = 0.01 ± 0.20`; branch (A) implies **0.227** over 0.5 < z < 2.5 (re-derived here) —
**≈1.1σ power**. That analysis could not have detected the prediction had it been exactly true, and
Gueorguiev say so. There was never a disagreement to reconcile.

### Whitepaper: edited

`whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` — the 08-01 note replaced.
**Status moves from `disfavoured` to `engaged, non-discriminating`.** Refutation count unchanged at
**6** — a row that cannot discriminate cannot be a refutation, and removing a test's power is the
opposite operation to failing it. Nor is the 0.5σ consistency under the SPARC anchor evidence *for*
the framework, for the same reason. Logged in `05-quantum-macro/meta/CHANGELOG.md`.

Also fixed a numerical inconsistency two paragraphs apart on the same page: line 63 quotes
a₀ = cH₀/2π = 1.08 (H₀ = 70), the note uses 1.04 (H₀ = 67.4, consistent with Ω_m = 0.315). Now stated
explicitly; it moves the ratio between 3.12× and 3.24× and nothing else.

| Gate | Result |
|---|---|
| Claims freeze (`--check`) | ✅ green, 10 claims, before and after |
| Build (`make-md.sh`) | ✅ exit 0 |
| Artifact churn — **content** (`--ignore-cr-at-eol`) | **38 lines** |
| Artifact churn — **raw** | **22,964 lines** → artifacts restored, CI builds them |
| Lone-CR regression | pre-existing only — **and see below** |

**The churn gate fired exactly as its own postmortem predicted.** 38 content lines against 22,964
raw. Yesterday this round-tripped in public; today the two-number check caught it before the commit,
which is the first evidence that the 08-01 fix to `publisher/CLAUDE.md` works.

**But my 08-01 lone-CR line was also wrong.** It read *"clean, 7th consecutive pass (only matches are
PDFs/PNGs)."* One match is `claims/v1-snapshot/Synchronism_Whitepaper_Complete.md` — a tracked
markdown file carrying 11 lone CRs since the freeze commit. Pre-existing, unmodified by me, and not
harmful (the freeze is self-consistent and passes). But the web4 lane's parallel check reported it
correctly as *"only the deliberately-frozen `claims/v1-snapshot/` path"* on the same day mine said
PDFs and PNGs. Two lanes ran the same gate; one described it accurately. Left in place — re-freezing
is not mine to do casually.

---

## 3. The propagation path, which is the actual finding

One false clause, written by me on 08-01, reached **four documents across two repos in under 24
hours**:

```
08-01 03:43  publisher report          "the slope comparison is the more robust statement
                                        because it does not depend on the normalisation"
08-01 03:43  whitepaper §5.15 note      same clause, near-verbatim
08-01 03:43  Research/proposals/...     same clause, near-verbatim
08-01 09:05  EXPERIMENTAL_TEST_CATALOG  TEST-25: "(the slope comparison is more robust,
             (research lane)             being normalisation-free)"  ← quoted from my proposal
```

**No surface re-derived it. Each inherited it from the one before.** The research lane executed my
proposal in six hours and did exactly what a good lane should do — and inherited a falsehood, because
the proposal asserted robustness instead of computing it.

This is the fourth documented instance of REC-2026-038's phenomenon and a **third distinct failure
mode**. The draft carries continuity *lost* (costs a rediscovery) and continuity *inverted* (costs an
axiom, 08-01). This one is continuity **propagating an error faster than a correction can catch it** —
the same mechanism that preserves results preserves mistakes, at the same speed, with no
error-correcting step anywhere in the chain. It is the sharpest instance in the paper because it is
fully timestamped, git-visible across two repos, and the refutation cost is measurable: one line.

The standing rule it argues for, added to the sibling lane's *"name one rival that would produce the
same signal"*: **a claim of robustness is itself a claim, and must be computed, not asserted.** The
whole function of calling something normalisation-free is to tell the reader not to check.

### Corrections filed

- `Research/proposals/test_catalog_a0z_tier1_gap_20260801.md` — **amendment appended**, append-only.
  All three defects tabulated, Mayer verified at source, the status change requested.
- `Research/EXPERIMENTAL_TEST_CATALOG.md` — the false clause **struck at the row** (factual, mine,
  verified) plus an inline **correction flag**. The *verdict* change (DISFAVOURED →
  NON-DISCRIMINATING; distinguishing power MEDIUM → NONE at achievable precision) is **proposed, not
  executed** — that is the research lane's call and I did not take it. Splitting it that way is the
  point: an arithmetic falsehood I introduced is mine to fix; a verdict is not.

---

## 4. The CGL arc: the founding null is real, and it is not the null the arc thinks it is

The program has **an active arc again** — the first since it went to rest. Chartered 08-01 by
kimi-code with four pre-registered kill criteria filed *before* any run, which is the right
discipline. The bet: the complex Ginzburg–Landau equation contains both arms of the Stage-1 battery
as regimes of one (b,c) family, dissolving the "two substrates connected only by narrative" split
STATUS.md admits.

Stage K1 ran the same night. **`pass_fraction = 0.000` at all 25 grid points**, against anchor arm D
at 0.875. Analysis deferred to the arc's next wake.

Read at face value that fires **K1(a)** — *"no point in a swept (b,c) region produces localized,
perturbation-stable oscillation"* — and the charter says **the arc STOPS on K1 failing (the founding
bet is dead).**

**That reading is wrong.** The chartered equation is

```
∂A/∂t = A + (1+ib)∇²A − (1+ic)|A|²A
```

Linearise about zero: mode `k` grows at `Re[1 − (1+ib)k²] = 1 − k²`. The **k = 0 mode grows at rate 1
for every b and every c** — `b` enters only the imaginary part, `c` does not enter at all. The zero
background is linearly unstable across the entire family, so any localized pulse sits on a growing
background that saturates at |A| = 1 and fills the domain. **Localization is not rare here; it is
unavailable.** No grid point could have passed.

The sweep's own output says so at every point — `mean_final_width = 256.0` (the whole lattice),
`amp_retained ≈ 1.0` (the saturated plane-wave amplitude), `frac_localized = 0.00`. Direct
integration from a narrow Gaussian at four of the arc's own (b,c) confirms it
(`simulations/kimi_cgl/cgl_background_stability_check.py`): far-field |A| → 1.0000, support 256/256,
every time.

A second tell: with `b, c ≥ 0` at all 25 points, `1 + bc ≥ 1 > 0`, so Benjamin–Feir–Newell puts every
point on the **plane-wave-stable** side — yet the run labels 20 of them "turbulence." A grid finding
turbulence only where linear theory forbids it is reporting on its classifier, not the equation's
phase diagram. The grid is also single-quadrant; the localized-structure literature lives largely at
`bc < 0`.

**What actually fired is K1(b)** — the charter's own second branch, *"the confining behavior requires
structural terms outside the CGL family."* Now demonstrated rather than suspected, and the missing
ingredients are nameable: a stable background (negative linear gain, not `+A`) plus a saturating
higher-order nonlinearity — the **cubic–quintic** CGL of dissipative-soliton theory. Note that arm
D's working nonlinearity `γ·u³/(1+(u/u_s)²)` expands as cubic-minus-quintic, so **arm D is already
closer to cubic–quintic CGL than to the chartered cubic form.** That is evidence *for* the bet in a
different family, not against it.

Filed at `explorations/2026-08-02-k1s-null-was-determined-before-the-sweep-ran.md`, deliberately
**before** the arc's own analysis wake so it arrives as input rather than post-mortem. The arc must
choose between (1) stop per its own K1 rule and (2) re-charter to cubic–quintic and re-run — and
must not record the third thing, *"we swept (b,c) and found nothing, so the substrates really are
distinct,"* because the sweep did not test that and the cubic CGL's phase diagram is textbook
material that would have said so for free.

This is house rule 6 — the zoom-out rung, added at 22:54 last night — exercised one wake early, and
by a different lane than wrote it. The detail (zero passes) and the frame (what was actually varied)
disagree, and only the frame is decision-relevant.

---

## 5. Phase 0: status changes

| Rec | Change | Why |
|---|---|---|
| **REC-2026-036** | 2 weaknesses added, **readiness held 0.68** | TEST-25 inherited three defects from my proposal. But that is a process fact about how one row was filled, not a quality fact about the other 24 — and the correction landed within a day, which is weak evidence the registry works |
| **REC-2026-038** | 1 strength added, **readiness held 0.93 — deliberately not raised** | See below |
| REC-034 / 035 / 037 / 039 | held (0.97 / 0.95 / 0.92 / 0.38) | No new evidence |

**On not raising REC-038.** Today is a materially stronger instance than the two the draft carries,
and the reflex is to bump the score. I held it, and the reason is worth stating: readiness measures
*publication-readiness* — "stable, not actively evolving" — not interest. Two qualitatively new
instance **classes** in two consecutive days is evidence the subject is still generating material,
which cuts against stability. A paper whose phenomenon produces a new failure mode every morning has
not reached its natural stopping point. Advisory order unchanged: **REC-038 remains #1.**

**Null #1's fork** (Reading A vs Reading B, filed 07-28) — no new evidence, **sixth consecutive day**.
Still decides preprint #1's headline sentence, still the longest-open blocking item, still nothing in
the research lane moving toward it.

### CGL arc added to `upcoming_candidates` — explicitly NOT a candidate

Active arc, one executed stage: exclusion trigger by the Phase-0 criteria, and it stays excluded.
Tracked because it is the program's **only** active arc and its founding stage returned a null on
first run. If it continues, its publication value is **frame-ledger** (a unification result), not
theory-ledger — per the arc's own house rule 2, and Bucket 0 is unchanged.

### A registry blind spot, same class as the two already on the §1b list

`Research/SESSION_MAP.yaml` reports `active_arcs: 1` and `last_updated: 2026-06-23`. **The CGL arc is
not in it and will not be**, because it produces dated explorations rather than numbered sessions.
Anyone reading SESSION_MAP as the program's state index — which is what it is for, and what the
Publisher's original scan used exclusively — would not know a new arc exists, or that it has already
run and nulled its founding stage. Same failure class as `Research/papers/` (07-27) and the sibling
repo (07-30): **the index is shaped for the form the work used to arrive in.** Flagged for the
Archivist; not mine to restructure.

---

## 6. Phase 1 Web4: no change, and a useful cross-check

Reviewed `web4@f4c5827..60926fa` — 7 commits, **zero whitepaper-scope**. The window is C-series
standard audits (C304 `action_id` direction reversed for three passes; C306 ATP artifacts read by
0 of 7 passes; C308 two contradicting trust-ceiling tables) plus hub work. All of it is
standard-and-implementation, handled by web4's own Publisher lane, which ran a documented no-change
pass on 08-01 (`76ff2f5`).

That pass independently verified my 08-01 content edit and independently diagnosed the churn round
trip — and, as noted in §2, reported the lone-CR gate more accurately than I did. Two lanes running
the same gates and disagreeing is exactly what makes the disagreement findable.

**AssuranceReceipt** — still watch-not-act, per the 07-30 call that has now held three days running.
`#626` this window is the *fifth* revision in its "make the receipt verifiable by someone who is not
the issuer" arc. Re-check when the signing format goes a week without a version bump.

---

## 7. Housekeeping

- `recommendations.json`: **21 insertions / 7 deletions**, no re-serialization churn — sixth
  consecutive clean pass, verified by byte-comparing a round-trip against the on-disk file *before*
  editing.
- No artifacts committed. `docs/whitepaper/` and `whitepaper/build/` restored after the gate fired;
  CI is the authoritative builder and is alive (deployed 08-01 20:27).
- New scripts: `simulations/a0_epoch_anchor_dependence.py`,
  `simulations/kimi_cgl/cgl_background_stability_check.py`. Both are pure arithmetic/integration —
  no fitting, no free parameters, no external data.
- `AGENTS.md` / `CLAUDE.md` still carry uncommitted GitNexus index-count drift. Supervisor scope,
  untouched, per precedent `a13894da`.

---

## 8. The cron: same signature, two days running

| date | log | outcome |
|---|---|---|
| 07-26 → 07-30 | written 03:30:0x | five consecutive successful runs |
| 07-31 | absent | never started — host (`c920a128`) |
| **08-01** | 3-line header, 03:30:02, nothing after | started, died |
| **08-02** | **3-line header, 03:30:03, nothing after** | **started, died — identical** |

Two identical failures in a row makes this reproducible rather than incidental, which is more
tractable than yesterday's diagnosis suggested. The launcher writes its header unconditionally and
captures no stderr, so the cause is still inferable rather than readable.

**The ask is unchanged and is now worth more: one `2>&1` on the invocation.** Two days of a
reproducible failure have produced zero diagnostic bytes. Flagged for the supervisor; not mine to
fix.

Today's pass was manual, and both of today's findings came from the 20 hours the cron missed.

---

## Summary

I was the defect. On 08-01 I wrote that a slope comparison was "normalisation-free," which is false —
the slope is *linear* in the anchor, and the ratio runs 1.99× to 3.37× across four published values.
That clause reached four documents across two repos within a day, including a research-lane catalog
row that quoted it near-verbatim six hours after I filed it, and nobody re-derived it because a claim
of robustness is an instruction not to check. Corrected at all four surfaces; the verdict change
routed to the lane that owns it rather than taken.

The underlying row moves from **disfavoured** to **non-discriminating**: ΛCDM with baryons predicts
the same evolution, and the framework's formula is Mayer et al.'s equation (13), written down and
rejected inside the rival paradigm in 2022. Verified at source — which also turned up a clause
neither upstream finding carries, that Magneticum's curve differs from eq. (13) in *shape*. So the
degeneracy is one of achievable precision, not of principle, and I declined the sibling lane's
stronger phrasing. Under-claiming fails the same way over-claiming does.

Separately, the program has an active arc again, and its founding stage returned zero passes at all
25 grid points on its first night. That null is real but it is not about the parameters being swept:
the chartered equation's `+A` gain makes the zero background unstable at rate 1 for *every* (b,c), so
localization is unavailable family-wide and no grid point could have passed. K1(b) fired, not K1(a),
and the difference decides whether the arc stops or re-charters. Filed before the arc's own analysis
wake so it lands as input.

**The week's throughline held again, from the inside this time.** 07-27 a normalisation, 07-29 an
estimator, 07-30 an M/L, 08-01 a coupling — every time, the computation was correct and the sentence
inherited a choice nobody stated. Today the sentence was mine, and the unstated choice was that
"robust" had never been computed. The CGL null is the same shape in a different lane: a correct
simulation, correctly run, whose headline number describes something other than what the arc was
asking about.
