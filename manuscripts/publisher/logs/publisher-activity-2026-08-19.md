# Publisher Activity Log — 2026-08-19

**Run**: autonomous daily pass, RUN-ID 21856. **Window**: 24h (08-18 03:42 PDT → now).
**Archivist context**: 2026-08-19 09:30 UTC. 0 new Synchronism numbered Sessions (**12th**
consecutive deliberate zero), counts re-derived at HEAD `64c1a0637` (strict 2888 / loose 2894).
14 new SAGE raising sessions. Two things worth carrying: the Archivist retired a threshold it had
invented *earlier in the same run* (a guessed `≥5 reuses = bank` cut, replaced by a measured
bimodal "hapax" boundary) — a chosen threshold and a measured boundary look identical in the output
and differ completely in what they license. And it avoided a **phantom Session 107**: this window's
exploration filename ends `...session107-is-a-different-refuted-mechanism.md`, which a filename
substring parse would have logged as a new session *and* as a duplicate-number anomaly against the
real #107. Settled with `--diff-filter=A`. A mention is not an instance.

---

## Sequence

1. Read `publisher/CLAUDE.md`, own launcher log header (RUN-ID 21856 — mine, so this is *this* run,
   not a prior failure), Archivist log, collective log, all three state files.
2. §1b surface scan, all seven surfaces. Sibling repos at `origin/<default>`, refs named.
3. Located the window's material: the **DE-sector locality fork**, executed.
4. Verified its identities independently in SymPy before touching anything.
5. Found the pass's larger defect while checking where the finding's consequences needed to land.
6. Four in-place corrections, two state updates, report, logs.

---

## Surface scan (§1b) — seven surfaces, refs named

| Surface | Ref | New in window | Verdict |
|---|---|---|---|
| `Research/SESSION_MAP.yaml` | HEAD | Archivist re-derivation only | 0 new sessions |
| `Research/papers/` | HEAD | 0 | — |
| `Research/proposals/` | HEAD | 1 (`de_locality_fork_perturbation_channel_factor_two_20260818.md`) | **material** |
| `Research/preregistrations/` | HEAD | 0 | — |
| `explorations/` | HEAD | 2 (DE locality fork triage; hive-organs wake) | **material** |
| `synchronism-site/explorer/findings/` | `main` @ `7e77d1e` | 1 finding + 4 scripts/outputs | **material, executed** |
| `synchronism-site/explorer/topics/` | `main` @ `7e77d1e` | 0 | — |
| `web4` | `origin/main` @ `184dd61` | 10 commits, 1 touching `whitepaper/` | **structural zero** (see below) |

**The 2026-08-09 rule earned its keep this pass.** The `web4` working tree is parked on
`cbp/concepts-normative-home` at an older commit. A bare HEAD scan would have read **zero** commits
for the window; `origin/main` carries **ten**. This is the second consecutive pass where naming the
ref changed the input rather than merely the wording.

---

## 1. The window's physics: the DE-sector locality fork, executed — and verified here

The dark-energy sector's last un-evaluated degree of freedom was that every published calculation
(Sessions #100/#101, the 08-11 covariant audit, the 08-12 direct fit) evaluated `C` at the
**background mean** density, while the framework's own one-equation postulate requires the **local**
density. The site explorer executed the fork on 08-18; the research lane read it, verified it and
routed it the same day.

**Recorded because it is the thing that used to fail:** on 2026-08-03 this file logged two halves of
one argument written the same morning in two repos by authors who could not see each other. This
time each half cites the other. That failure mode did not recur.

**Verified here before propagation** — `simulations/publisher_20260819_de_locality_identities.py`,
symbolic, exact, all `x > 0` and `γ > 0`, 4/4 assertions pass:

| # | Identity | Result |
|---|---|---|
| 1 | `C = ((1+x)^{2γ} − 1)/((1+x)^{2γ} + 1)` | exact |
| 2 | `ρ_DE/ρ_crit = 2x/((1+x)^{2γ} − 1)` | exact |
| 3 | at γ = ½: `ρ_DE/ρ_crit = 2` **exactly**, `dρ_DE/dρ_m ≡ 0` | exact |
| 4 | `δ_DE/δ_m = 1 + w_DE = ε(ln(1+x)/x − 1) + O(ε²)`, `ε ≡ 2γ − 1`; ≡ 0 at γ = ½ | exact |

(3) is the substantive one: ρ_DE is constant in **density**, therefore in **space** — the same inside
a void, a cluster, a disk and a neutron star. So at γ = ½ the branch is ΛCDM **non-perturbatively**,
not merely degenerate with it in a background fit. (4) says the clustering channel is not independent
of the expansion-history channel; both are linear in the same single ε, with no order-ε⁰ term.

### The correction I made against the source rather than relaying it

The finding states *"the locality horn improves the required precision by 7%."* It does not.

`σ_γ = (½ − γ̂)/3` is **horn-independent by construction** — it prices the γ-space separation between
the Λ-point and the MOND-point, which is exactly how §5.15 already derives it. So `0.004` at
γ̂ = 0.487 and `0.0037` at γ̂ = 0.489 are **one quantity rounded twice**, and locality buys **0%**
there. The 7% is a central-value shift wearing a locality label.

The real ×2 is in the observable, and it is exactly 2: read off the source's own Horn N / Horn L peak
table (0.223% vs 0.444%, ratio → 2.00 as γ → ½), the required **σ(fσ₈)** relaxes 0.074% → 0.15%.

**Verdict untouched** — DESI DR2 delivers ~1–3% per bin, so ~10× underpowered, and the parameter the
test would measure sits at **ε = −0.022 ± 0.220 = 0.10σ from exact ΛCDM**. This is a precision
correction, not a verdict change. It is worth making because under-claiming fails the same way
over-claiming does, and because the 7% had already propagated into `PREDICTIONS.md`.

---

## 2. The pass's larger finding: a correction that reached a banner and never reached the body

Checking where the fork's consequences needed to land took me to
`Research/EXPERIMENTAL_TEST_CATALOG.md` — the document **REC-036 is**, and the one addressed in its
own text *"For the Publisher: This document is a message."*

It carries the **S674/S675 census as a banner at line 1**: *5 executed → all collapsed*
(TEST-04/08/14/15/18), *"of the 9 untested, **0** have a verified first-principles-derived amplitude"*,
TEST-07 dimensionally inconsistent, TEST-17's numbers don't follow and are Lorentz-excluded at ~11
orders of magnitude. **Its body was never updated against it.**

The closing section, *Summary: The Decisive Tests* — the part a reader acts on:

| Rank | Test | Status in this document's own header |
|---|---|---|
| **1** | TEST-04 BAO Coherence Modulation | executed → **collapsed** |
| **2** | TEST-14 Wide Binary Density Dependence | executed → **collapsed** |
| 3 | TEST-11 EEG Anesthesia Phase Transition | untested — but in the 0-of-9 no-derived-amplitude set |
| **4** | TEST-15 GW Speed–DM Column | executed → **collapsed** |
| **5** | TEST-07 Cosmic Interference Patterns | 500 Mpc **dimensionally inconsistent** (m² not m, S632) |

**None of the five survives the header's own criteria.** And *What "Positive" Would Mean* named
TEST-04, TEST-07, TEST-14, TEST-17 — **4 of 4 disqualified**.

TEST-04 is the sharpest: it held the **#1 decisive slot for 15 months** while Session #107 — the very
session it derives from — predicts BAO matching ΛCDM at *exactly 0%* in all five bins; its 10⁻⁴
amplitude has no session-level derivation; and its 10⁻⁵ falsification threshold sits ~3,000× below
DESI's best precision, so it is unfalsifiable as written. It was withdrawn 2026-05-04.

**Fixed in place.** Each entry now carries its census verdict inline, and the four names under *What
Positive Would Mean* were **removed rather than reranked** — because the catalog holds no surviving
candidate to promote. The honest operative statement is now inscribed: *the program's
forward-looking decisive-test list is empty; it needs regeneration, not repair.*

One ledger disagreement recorded rather than smoothed: the catalog files **TEST-07 under "genuinely
untested"** while `PREDICTIONS.md` Bucket 2 files it under **refuted**. Refuted is the stronger and
better-sourced reading — a dimensionally inconsistent scale is not awaiting data.

**Why this matters beyond one file.** This is the 2026-08-14 propagation split — *corrections reach
authors, not artifacts* — observed **inside a single document**, header versus body. The correction
landed where it was cheapest to write it. Convergent evidence arrived independently this window from
the web4 lane, which spent four commits (C402/C404/C406/C408) finding the same shape in its own
specs, including *"the fix that closed two rows widened the requirement list and repaired the query,
and never touched the only data the spec has."* Two lanes, two repos, no shared instruction, one
finding.

---

## 3. Session #107's primary now carries its errata

`Research/Session107_DESI_Forecasts.md` read **"Status: COMPLETE"** with "3.1σ / 3.2σ per bin, 6.6σ
combined" and **no erratum**, against three independent corrections none of which had reached it:

1. **2026-05-05** — DESI DR1 full-shape disfavors it at every LRG bin; by the session's *own*
   falsification ladder, ΛCDM is favored at every bin and at the combined fit.
2. **2026-07-14** — the registered criterion was met at only ~1.5σ; the quoted 2.4σ is a
   GR-conditioned σ₈-amplitude statistic, and DESI's own MG analysis puts μ₀ within 1σ of zero.
3. **2026-08-18** — the mechanism is `G_local/G_global = C_cosmic/C_galactic < 1` plus an **assumed**
   σ₈(z=0) = 0.76: the door-#1 local-density gravity law applied to growth, refuted on the SPARC RAR
   2026-08-15, resting on a `C_cosmic ≠ C_galactic` distinction withdrawn 2026-08-11.

The research lane flagged this as *"an archive-lane erratum"* and left it there. The primary now
carries all three, with the genuine DE forecast (−0.22%, 0.10σ) beside the 6.6σ and the **173×**
disagreement stated — two mutually inconsistent cosmological growth predictions that sat in one
archive as two models never recorded as different, and **the large one is the refuted one.**

No forecast table was altered; they stand as the December 2025 record.

---

## 4. Whitepaper §5.15 — one word removed, because it was load-bearing in the wrong direction

The standing sentence read *"the γ = 1/2 branch is still exactly the ΛCDM **background**."* That
qualifier was not inert. The site's 08-18 visitor Pass 4 read the matching *"background-only, no
perturbation sector"* caveat as an unexplored channel and filed **"either horn closes TEST-26 in
2026"** — a proposal to convert a dp-gated 2027–28 registration into a near-term test. The execution
refutes exactly that.

Strengthened in place (not appended): the branch is ΛCDM **non-perturbatively** — same background,
same linear perturbations, same nonlinear field configuration — with the four verified identities,
the ×2-in-σ(fσ₈), the 0%-in-σ_γ correction, and the pressure-sector conditionality carried (the two
defensible δp_DE treatments disagree in **sign** at z = 0, and the covariant completion that would
settle it has no surviving member — the forks are nested, and the innermost one is already empty).

**Propagation checked by opening the matches, not counting them** (08-05 lesson): `background-only`
and `perturbation sector` occur in whitepaper sources **only** in §5.15. No executive-summary or
conclusion edit is owed.

CHANGELOG entry added to `05-quantum-macro/meta/`. Count **HELD at 6**; Bucket 0 = 0.

---

## 5. Web4 — 11th consecutive structural zero, and a live divergence watch that did not trip

`origin/main` @ `184dd61`; HEAD parked on `cbp/concepts-normative-home` (stated, per the rule).
10 commits, exactly **one** touching `whitepaper/` — `4c6e6a8`, the web4 lane's own Publisher run log
appended to `PUBLISHER_CONTEXT.md`. **Zero normative whitepaper edits.**

R7 trigger unchanged and correct: **registry publish** (crates.io / PyPI / npm), restored 08-18 after
two passes on a criterion my own record contradicted. Hub F0 completion is POC and does not fire it.

**MRH divergence watch (opened 08-18): checked, NOT tripped.** Verified rather than asserted —
`forum/mrh-relevance-contract-2026-08-17.md` untouched since its creating commit; `git grep` over
web4's `whitepaper/` at `origin/main` for the contract reading returns **nothing**; Synchronism §4.2
still carries *"Spatial Boundary: the physical distance beyond which influences become negligible."*
Trip condition unchanged: act only if one whitepaper adopts the contract reading while the other
keeps the radius reading.

---

## 6. Standing escalations

- **Site maintainer 401, day 7.** `maintainer/logs/` shows `Failed to authenticate. API Error: 401
  OAuth access token is invalid.` on 08-16, 08-17 and 08-18 (and 08-12/13/14). This is a dead
  `CLAUDE_ADMIN_TOKEN`, **not** the account outage — the explorer ran a full session on the same host
  on 08-18. **Owner action; it will not self-heal.** The backlog it blocks grew again this window:
  `/dark-energy` line 195 still publishes *"The sector is background-only. There is no perturbation
  sector"*, which the site's **own** 08-18 finding declares false (§6b). Production continues to
  outpace correction because the correction lane is offline, not because the corrections are hard.
- **REC-040**, 16 days on a sole blocker (external prior-art walk), **held at 0.55**. Priced at one
  agent / one document / one day on 08-18. Deliberately **not** re-escalated: this file restated the
  same unrun gate three times before 08-18, and a fourth restatement adds nothing the pricing settled.
- **Preprint / two-paper strategy**, and the TEST-09/TEST-10 count-collapse question: still dp's call.

---

## Verification

| Check | Result |
|---|---|
| `simulations/publisher_20260819_de_locality_identities.py` | 4/4 symbolic assertions pass |
| `whitepaper/make-md.sh` | exit 0 |
| `whitepaper/make-web.sh` | exit 0 |
| Churn gate — content (`--ignore-cr-at-eol`) | 98 ins / 500 del |
| Churn gate — raw | **11,984 ins / 12,386 del** |
| Action | **artifacts restored, not staged** — `git checkout -- docs/whitepaper/ whitepaper/build/` |

The raw number exceeding the content number by two orders of magnitude is the documented CRLF failure
mode, exactly as the rule predicts. CI is the authoritative builder; only sources are committed.

**And the same gate fired a second time, on a SOURCE — which is what the 08-10 amendment was for.**
Staging `05-quantum-macro/meta/CHANGELOG.md` normally produced **68 raw lines against 4 of content**.
The cause is not my edit: that file's committed blob is **mixed** (32 CRLF lines, 41 LF), the repo has
**no `.gitattributes`**, and `core.autocrlf=input` therefore rewrites all 32 CRLF lines to LF at
`git add` time. So a four-line append to a CHANGELOG silently carries a 32-line ending rewrite, and
the raw/content gate is the only thing that shows it. Staged with `git -c core.autocrlf=false add`;
final staged diff is **RAW == CONTENT == 543 / 14 across all 11 files, zero churn.**

This is the first *measured* source-side instance of the hazard the standing `.gitattributes` proposal
describes (filed 07-26, restated 07-31, dp's call). The 08-01 round-trip priced the artifact side at
22,410 content-free lines; the source side is small per event but silent, and it fires on every append
to any mixed-ending file. Recorded, not escalated — it remains a one-time-normalization decision that
is dp's to make.

---

## Ledger

**READMEs updated**: 0 · **Publication candidates**: 0 new · **Whitepaper proposals**: 0 Synchronism
(1 edit made directly, minor/in-place), 0 Web4 (11th structural zero) · **Refutation count**: HELD at
6 · **Bucket 0**: 0.

**Files changed**: `Research/EXPERIMENTAL_TEST_CATALOG.md` (propagation fix),
`Research/Session107_DESI_Forecasts.md` (erratum), `PREDICTIONS.md` (7% correction),
`whitepaper/sections/05-quantum-macro/15-dark-matter/dark_matter.md` (+ section CHANGELOG),
`simulations/publisher_20260819_de_locality_identities.py` (new), three state files, this log,
the daily report.

---

## So what?

The physics was a null and the lane that produced it already knew that. What this pass adds is two
things the null could not: the **7% was not a gain from locality** — a horn-independent quantity was
credited to the horn — and, more importantly, **a document can be corrected at the top and still be
wrong all the way down.** The test catalog kept a withdrawn test in its #1 decisive slot for fifteen
months *underneath a banner saying it had collapsed*, and the same day two independent lanes in two
repos found that same shape in their own specs. The rule this pass earns is narrow and checkable:
**a correction has a landing site, and the banner is not it.** Ask which section a reader would act
on, and check that one.
