# Publisher Daily Report - 2026-08-03

## Verdict: SIGNAL, twice. (1) The program proved its galaxy sector is *exactly* MOND under one substitution, and empirically bounded that substitution at ≤25%, in two repositories 48 minutes apart, with neither author able to see the other. (2) The cron has been healthy all along; the two "failures" I reported were me reading my own run's opening banner, an error I had already caught and correctly written down on 07-31.

**This is an autonomous run** — and so were 08-01 and 08-02. See §5, which is a correction to
three of my own documents.

---

## 1. WAKE: the largest window in weeks, and the two halves of one argument landed in different repos

Nine commits since my last pass, across two repositories:

| commit | repo · time (PDT) | what |
|---|---|---|
| `647edff` / `061638d` | site · 06:10 | γ=1/2 MOND identity back-annotated; **EFE mechanism retracted** |
| `581d5cb7` | Synchronism · 06:09 | same, back-annotated here |
| **`c4a7f29`** | **site · 08:18** | **RAR scatter no-go EXECUTED on real SPARC — local density carries zero information** |
| `d7ed36fb` | Synchronism · 09:06 | triage: γ=1/2 exact MOND point; EFE=0; a₀(z) NON-DISCRIMINATING |
| `06e1dbb6` | Synchronism · 09:58 | **K1 verdict: instrument repaired, re-swept, KILL — the CGL arc stops per charter** |
| `007f7bfc` | Synchronism · 13:19 | **the alignment arc chartered** (dp's common-substrate frame) |
| `098704a5` / `66bd3b35` | Synchronism · 17:19/17:38 | omega decomposition: plan, then a clean **negative** |
| `022ada55` | Synchronism · 03:22 today | publisher launcher instrumented — timezone, UTC, RUN-ID |

An arc closed and a new one opened in the same afternoon. Three registered dead ends. And the
program's sharpest galaxy-sector result of the year, which nobody has yet noticed is one result.

---

## 2. One substitution, and its exclusion, 48 minutes apart

The research lane's γ = 1/2 triage closes on this sentence:

> Every other galaxy result (BTFR slope, DM-fraction ceiling, RAR shape) is downstream of that single
> g_bar→ρ substitution failing — **which is precisely the locality no-go.**

A conditional, written in the voice of something still to be settled. **Forty-eight minutes earlier,
in the sibling repository, that antecedent had been executed on real SPARC data, form-free, with a
validated positive control.** The triage does not cite it. Nothing in the Synchronism repository did
— I grepped the finding's own distinctive strings (`informationally empty`, `partial correlation`,
`scatter no-go`) across `Research/`, `explorations/` and `PREDICTIONS.md`: nil.

### The identity — re-derived here, not inherited

At γ = 1/2 the Hill identity collapses:

```
tanh(½·ln(1+x)) = x/(x+2) = μ_simple(x/2)      identically, for all x
max |deviation| over x ∈ [1e-6, 1e6] = 5.55e-17   (machine precision)
```

Since the framework defines `f_DM = 1 − C`, **C occupies exactly MOND's μ seat**. So the whole galaxy
sector is one statement: MOND, with μ's argument swapped from enclosed-mass acceleration `g_bar` to
local density ρ. Not "MOND-like on several observables" — *one substitution.* And the SPARC
form-selection table's preferred form is Hill n = 1 exactly, so the data **select** that point.

I re-verified rather than accepted, per this lane's own 08-02 rule that a claim of exactness is a
claim. The degeneracy caveat is the part that costs the program a number it quotes: the O(x²) term
separating γ from ρ_crit carries coefficient γ(2γ−1)/2 — computed here as **−0.0319 at γ = 0.425,
−0.0054 at 0.489, exactly 0.000 at ½, +0.060 at 0.600.** It vanishes at essentially the SPARC-preferred
value, so "γ ≈ 0.489 preferred by SPARC" is degenerate with the ρ_crit prior and is not a standalone
measurement of γ. That qualifies the TEST-11/Cassini interval 0.425–0.600 wherever it is cited as a number.

### The exclusion — 2,622 SPARC points, no functional form, no fitting

| conditioner | robust σ(log B_req) |
|---|---|
| log g_bar (MOND/RAR) | **0.1178 dex** — reproduces Lelli+2017's published 0.11–0.13 |
| log ρ (this framework) | **0.1613 dex** |

The 1.37× excess **rises to 1.77×** when the framework is given 175 free per-galaxy offsets — so it is
*within*-galaxy, the radial shape of the density profile, not a calibration error. RAR residuals carry
**r = +0.0012** with local-density residuals at fixed `g_bar`, 95% CI **[−0.076, +0.081]**, ≤0.7% of
variance, against a pipeline that recovers injected r = 0.10/0.20/0.30. And the constructive bound:
interpolating `log u_α = (1−α)·log Σ + α·log g_bar` gives α ≥ 0.75 at 95% — **SPARC admits at most 25%
weight on a local-density variable. This sector sits at 100%.**

The escape hatch is narrowed with a reason rather than asserted shut: no exponential kernel on Σ at any
range ≤30 kpc reaches `g_bar` performance (best 1.21× at λ ≈ 7 disk scale lengths, then degrading),
because `g_bar = GM(<r)/r²` carries a 1/r² that no convolution of Σ can generate. That turns "make the
coupling differential" — the published BCM 2017 symmetron escape this program found on 07-27 — from a
free dial into a specific, checkable obligation.

### The scoping correction, which is mine and not the finding's

The finding's propagation note tells the maintainer this is novel because *"Lelli tested surface
columns, never volumetric ρ — this is that test, run."* But its own **§2** says the headline needs no
thickness model because at constant h `log ρ = log Σ − const`, and its own **Scope** block says the
volumetric claim is not tested independently of the surface-density one. And this program's prior-art
audit (`explorations/2026-07-23-locality-nogo-prior-art-audit.md`, row 13, read from the arXiv PDF)
records that **Lelli et al. 2017 already tested RAR residuals against stellar surface density at R and
found no significant correlation.**

| § | status |
|---|---|
| §4 decisive null | **quantification of published prior art.** The CI, the attenuation bound and the positive control are real additions. Discovery, no. |
| §6 smoothing scan / 1/r² obligation | **new — and the most valuable piece.** |
| §7 α-admixture bound ≤ 25% | **new, reusable, framework-independent.** |

Nothing here touches the verdict. It changes what may be cited as novel, and a draft that led with §4
would be scooped on arrival by a 2017 paper on the same dataset. **Refutation count held at 6** — the
sibling finding says so itself and gives the right reason (same underlying failure, executed form-free,
not an independent refutation). Bucket 0 unchanged (0). Filed at
`explorations/2026-08-03-publisher-one-substitution-and-its-exclusion.md`; flagged for the site lane,
not edited there.

---

## 3. Whitepaper: the "breakthrough discovery" is now a theorem, and that is the problem

`05-quantum-macro/15-dark-matter/dark_matter.md` has carried this since Sessions #86-89:

> A breakthrough discovery: MOND and Synchronism are the same physics in different parameterizations.

**That sentence is now exactly true, which is precisely why it is no longer a discovery.** Every prior
statement of the reparametrization verdict arrived from an audit or an external critic, in the hedged
form "reduces to MOND in the testable regime." This one is an *identity*, and it arrives from the
framework's own SPARC-preferred parameter value. Let the compander float and it lands on simple μ.
There is no residual to argue about.

Integrated as a status note under that heading, with the empirical half (≤25% admixture) beside it,
the scope paragraph above, and the two attached corrections: γ ≈ 0.489 is not a standalone γ
measurement, and the **EFE over-claim is retracted** — the galaxy sector has no field equation, so the
"0.3–0.4× MOND" figure was never derived from one; the actual structure gives **EFE = 0 exactly** (SEP
satisfied by construction), sharper and in tension with Chae et al. 2020's ~4σ SPARC detection.
Logged in `05-quantum-macro/meta/CHANGELOG.md`.

| Gate | Result |
|---|---|
| Claims freeze (`--check`) | ✅ 10 claims, v1 freeze verified, exit 0 |
| Build (`make-md.sh`) | ✅ exit 0 |
| Churn — **content** (`--ignore-cr-at-eol`) | **40 lines** |
| Churn — **raw** | **23,168 lines** → artifacts restored, not committed; CI builds them |
| Lone-CR | exactly one path, the frozen `claims/v1-snapshot/` — **10th consecutive pass** |

Third consecutive pass in which the two-number churn check produced the right action.

---

## 4. Phase 0

| Item | Change | Why |
|---|---|---|
| **CGL Arc** | `active` → **`closed_killed`** | K1 KILL, arc stopped per its own charter. See below. |
| **Alignment Arc** | **new**, `active` | Chartered 08-02, zero stages. Not a candidate (active arc). |
| **REC-2026-040** | **OPENED at 0.45 / MEDIUM** | The α-admixture bound as a framework-independent methods note. |
| REC-2026-039 | weakness added, **held 0.38** | Superseded as the lead instance of its own genre. |
| REC-2026-038 | strength + weakness added, **held 0.93** | Fifth instance, fifth failure mode — and it is §5 below. |
| REC-2026-035 | weakness added, **held 0.95** | An open question against its own stated blocker. |
| REC-034/036/037 | held (0.97 / 0.68 / 0.92) | No new evidence. |

**The CGL arc closed correctly, and my 08-02 input was engaged with rather than absorbed.** The arc
lane read `explorations/2026-08-02-k1s-null-was-determined-before-the-sweep-ran.md`, cited my
`cgl_background_stability_check.py` by commit, agreed the v1 instrument was contaminated by the rival
arm's vacuum-background ontology, **built a CGL-native replacement**, and re-swept. Repaired best point
0.042 against a 0.10 bar and a D-arm anchor of 0.875 in the same harness. Verdict K1 KILL — and rather
than declaring INDETERMINATE to stay alive, the arc stopped. Residual loopholes (the negative-c
half-plane crossing the Benjamin–Feir line, topological phase-slip ICs, the 1/24 weak spot) are named
and costed at ~2 h for a *new* charter, not used to keep this one breathing. That is the discipline
working, and it is worth saying so plainly: **a lane with no memory of writing the rule followed it
against its own founding bet.**

**REC-2026-040, opened at 0.45 and not higher, deliberately.** The α-admixture bound is the strongest
new publishable item in the window and the only *framework-independent* one in the whole file — it
survives Synchronism being wrong, which nothing else here does. But REC-2026-039 was opened at
0.72/HIGH on 07-29 and downgraded to 0.38 the next day when its result turned out to have been executed
in-repo eight weeks earlier. So: opened low, with the external prior-art gate named as weakness #1 and
the §4-is-prior-art scoping stated in the summary rather than discovered later. Blocking actions, in
order: run the external gate on the α-interpolation construction *specifically* (not on the null, which
is prior art); instantiate on a second variable (TEST-12 is one function call away); resolve authorship
and which repo owns the draft.

**Null #1's Reading A/B fork — seventh consecutive day** with no new evidence, still the longest-open
blocker on preprint #1's headline. Partially overtaken by events: the fork asks how to define locality
operationally, and ≤25% admixture is an operational answer in the one place it was blocking. It does
not settle the headline sentence.

---

## 5. The correction I had already written down

**The cron never failed. I did.**

| date | log | ends |
|---|---|---|
| 08-01 | 3,514 bytes | `Session Complete: 03:44:57 (claude exit=0)` |
| 08-02 | 4,704 bytes | `Session Complete: 03:49:20 (claude exit=0)` |

Both healthy. Both autonomous. I reported both as failures, in the daily report, the activity log,
`PUBLISHER_CONTEXT.md` §6 and the collective log, and escalated an OWNER-ACTION to the supervisor
asking for a `2>&1` **that has been on the launcher since 07-24.**

The mechanism: the agent reads its own clock as **10:30 UTC**, reads the header's bare **03:30** as a
local time, and infers a run that died seven hours ago. But 03:30 PDT *is* 10:30 UTC — it is the
opening banner of the run currently executing. And the log is header-only *whenever you look at it*,
because `claude -p ... | Add-Content` flushes the agent's output only at completion. **A live run
cannot present any other signature.** The test I used can only return "dead."

The 08-02 log's first content line, written by the cron-launched agent into the file its own live run
was producing:

> `[10:30] Manual pass begins (cron header at 03:30:03, no work — 2nd consecutive identical failure)`

A death certificate, signed by the deceased, in the deceased's own handwriting.

**And I had already caught this.** On 2026-07-31, in `private-context/publisher/log.md`:

> *"I opened this pass by declaring the cron a silent failure based on a 135-byte log, then found its
> mtime was one minute old: 10:30 UTC is 03:30 local... Anyone auditing Publisher liveness by log size
> should use the closing banner, not the byte count."*

Correct, complete, durable, on this track's own designated coordination surface — and it reached
nobody, including its own author, twice. The pre-07-24 header-only logs *were* real failures; the
`2>&1` fix on 07-24 made the signature non-diagnostic, and nobody noticed the detector had gone blind.
The same trap is recorded as misread four times in the Archivist track (07-26/27/29/30). The launcher
was instrumented at 03:22 today (`022ada55`) with timezone, UTC and a RUN-ID, which kills the ambiguity
at the source — I read that header, and it is why this pass caught it.

**Corrections filed** (append-only banners, not rewrites): `reports/2026-08-01-…`,
`reports/2026-08-02-…`, `logs/publisher-activity-2026-08-02.md`, `PUBLISHER_CONTEXT.md` §6, and the
collective log.

### Why this is the day's second finding and not just an apology

Put beside 08-02, it measures an asymmetry the program has been asserting:

- **08-02**: one false clause I wrote reached **four documents across two repositories in under 24
  hours**, unaided, corrected by nobody.
- **08-03**: one true correction I wrote reached **nobody in 72 hours**, including me, twice.

Same track, same surfaces, same week, opposite directions. **Errors propagate; corrections must be
chased.** That is a fifth distinct failure mode for REC-2026-038 — continuity that *exists, is correct,
is durable, and is never consulted* — and it is the sharpest instance the manuscript could carry,
because it is first-person, fully timestamped, and self-refuting on its face. Added there. Readiness
still held at 0.93: four qualitatively new modes in three days is evidence the subject is still
evolving, which is what readiness measures. My recommendation to that draft is now to stop enumerating
instances and write the asymmetry structurally — the instance list is what keeps moving; the asymmetry
is what is stable.

A sixth-mode candidate arrived the same day and points the other way: the CGL instrument repair (§4) is
cross-lane continuity performing **error correction** within twelve hours. Every instance the draft
carries is continuity failing. That would be its first positive one.

---

## 6. Phase 1 Web4: no change

`60926fa..e5f221a`, 3 commits, zero whitepaper-scope. C310 (`t3-v3-tensors`, 8th delta) reports the
spec **byte-frozen and substantively clean against itself for the 8th consecutive pass**, with all four
findings routed to owners and zero spec mutation — no canonical-terminology exposure (T3/V3 definitions
untouched; the findings are in the artifacts the spec *cites*). Plus a hub admin fix and web4's own
Publisher no-change log. Handled by that lane. AssuranceReceipt still watch-not-act, fourth day.

---

## 7. Housekeeping

- `recommendations.json`: **55 insertions / 9 deletions**, raw diff equals content diff — **seventh
  consecutive pass with no re-serialization churn**, round-trip byte-verified before editing.
- No artifacts committed; `docs/whitepaper/` and `whitepaper/build/` restored after the churn gate fired.
- `whitepaper_sync.json` was stale (`last_review` 07-30); updated.
- `AGENTS.md` / `CLAUDE.md` still carry uncommitted GitNexus index-count drift. Supervisor scope,
  untouched, per precedent `a13894da`.
- **Retired watch item**: the cron. **Retired watch item**: the CGL arc's stop-vs-recharter choice.
- **Standing, dp's call, fourth deferral**: `.gitattributes` (`docs/whitepaper/** text eol=lf`), now
  priced at ~23,000 lines per occurrence.

---

## Summary

The program established, on one morning, that its galaxy sector is *exactly* MOND under a single
substitution — γ = 1/2 makes the coherence compander identically MOND's simple μ, verified here to
machine precision, and the SPARC data select that value — and, in a different repository 48 minutes
earlier, that the substitution is permitted at most 25% of the way by the same SPARC rotation curves,
form-free, with a validated positive control. Those two statements are the premise and the conclusion
of one argument. Neither author could see the other. The sibling-repo scan added three days ago after
it cost a full pass is the only reason they are now in the same document. I also declined the source
finding's novelty claim for its own headline section: at constant scale height that test is Lelli et
al. 2017's, and the genuinely new content is the 1/r² obligation on differential couplings and the
admixture bound. A draft led by the wrong half would be scooped by a nine-year-old paper.

And the instrument I have been using to check my own liveness cannot return anything but "dead." Two
healthy autonomous runs were reported as failures across four documents and escalated to the
supervisor, after I had caught the identical error on 07-31 and written the correct rule into the
collective log. That record was complete, durable, correctly stated, and inert. Set against 08-02 —
where one false clause of mine reached four documents in a day with no help at all — the week has now
measured both directions of the same channel, and they are not symmetric. The repository preserves
errors at the speed of reading and corrections at the speed of someone deciding to look.
