# Publisher Daily Report - 2026-07-29

## Verdict: SIGNAL — the research lane's new capstone framing is over-strong on one word, and the headline number it prescribed to fix a real problem is the weakest-supported point on the curve. Both checked by computing on them. One of my own guesses was refuted in the process.

Autonomous cron run, **fourth consecutive successful start**.

---

## 1. Phase 0: the research lane produced a capstone, twice over, after yesterday's run

Two commits landed on 07-28 *after* my 03:42 pass (`56cbb1b0` 06:18, `8b4a7e66` 09:04). They inscribe
the **nested-submodel reframe** — presented as "the sharpest statement of the whole door-#1 cage":

> If the bounded boost is the framework's structural difference from MOND, the galaxy sector is
> **Synchronism(galactic) = MOND ∩ {B ≤ 1/Ω_m}** — a strictly nested submodel. It has two possible
> outcomes against its parent: statistically indistinguishable, or refuted. **It cannot win** — a
> priori, from structure, no SPARC data required.

Plus a self-caveat: the TEST-10 headline "69% of SPARC exceed the 68.5% ceiling" is
**convention-dependent** (a boost is a dynamical/baryonic ratio, whose cosmic value is Ω_m/Ω_b = 6.39,
not 1/Ω_m = 3.17), with `f_DM,max = 0.927 ⇒ B ≥ 13.7` prescribed as the definition-free replacement.

The self-caveat is a good catch of the same class the project has been catching all week — a number
outliving the computation that justified it. The reframe is where I want to spend this report.

### Propagation status: it hasn't, and this time that is the right answer

`git grep -i "nested.submodel"` returns hits in `PREDICTIONS.md`, `SESSION_FOCUS.md`, and the two
triage/proposal docs. **Zero** in the whitepaper, the docs build, or the preprint strategy proposal.

Yesterday the identical situation was a defect. Today it is not, and the difference is worth stating
because I nearly swept it on reflex:

| | 07-27 locality retraction | 07-28 nested-submodel reframe |
|---|---|---|
| What the whitepaper carries | the **falsified** universal sentence | S684's fit-XOR-discriminate fork — subsumed, **not contradicted** |
| Truth-status of published text | wrong as of that morning | unchanged |
| Correct action | propagate immediately | **do not propagate yet** |

Additive narrative integration follows batch cadence or dp's packaging decision, and the reframe's own
count questions (TEST-09+10 = one refutation or two? does TEST-11 join the headline?) **gate on dp**.
Inscribing my correction of a framing the research lane hasn't seen yet, into the reader-facing
document, would be yesterday's failure with the sign flipped again. **No whitepaper section edit
today** — the propagation surfaces I did sweep are the proposals and the state files.

---

## 2. The correction: "a nested submodel cannot win" conflates fit with selection

The nesting fact is right. The inference is wrong, and wrong in a standard way.

**A nested submodel cannot achieve a better *fit*. It can absolutely be *selected*.** That is the
entire basis of parsimony. The canonical example is in the same field as the claim: ΛCDM is a strictly
nested submodel of wCDM (w = −1 is a restriction of free w). It never fits better. It wins anyway, on
every evidence-based criterion, because the restriction costs nothing and sharpens the prior
predictive. "ΛCDM cannot win because it is nested" would be recognised instantly as a category error.

The boost ceiling is exactly that kind of restriction: **B_max = 1/Ω_m is fixed by cosmology, not
fitted to rotation curves.** It adds no free parameter. So what nesting buys a priori is a
**dichotomy**, not a verdict:

| ceiling is… | outcome | adjudicated by |
|---|---|---|
| slack across the data | tie on fit, **win on selection** | model comparison — a real, favourable result |
| binding somewhere | **refuted** | the constraint is falsified |

Which branch obtains is **empirical**. It is answered by SPARC and by nothing else.

**The verdict does not move** — the ceiling binds, hard, and this pass re-verified it. What moves is
what the executed work *was*. The triage concludes the result "was never in doubt" and that every
execution was "a corollary confirming where the ceiling binds." Under the corrected reading, the SPARC
runs were **the decision procedure between two genuinely different outcomes, one of them favourable**.
That is the strongest available answer to "why run it at all," and the reframe as written throws it away.

**A dependency that cuts both ways, and I state it because a referee will find it.** The nesting
premise needs the galactic law to be a restriction with no compensating free parameter. `C(a)` has Ω_m
fixed and a₀ free (as MOND's is) — but **φ, which this archive's own 2026-07-17 provenance audit found
is fitted-then-named, not derived**. If φ is genuinely free, the law is a two-shape-parameter family
that happens to be bounded, not a restriction of a one-parameter parent — and then the nesting premise
fails outright. Either way **"it cannot win" is not established by nesting**: under fixed φ the nesting
is real and the conclusion is wrong; under free φ the conclusion may hold and nesting is not why.

Filed as `Research/proposals/nested_submodel_fit_versus_selection.md` **with a pre-registered
falsifier**: if the program's operative comparison turns out to be fit-only with no complexity term
anywhere, "cannot win" is true under that convention and my objection is terminological. *Evidence
against that reading:* the claim itself uses the word "**select**" ("0 tests could select Synchronism
over MOND+EFE+ΛCDM"), and the archive's own comparisons apply complexity-penalised selection (S559:
"AIC/BIC both favor gradient").

---

## 3. The executed part: the prescribed headline statistic is the weakest point on its own curve

The triage was right that "69%" is convention-dependent and right to want a definition-free
replacement. I did not think `f_DM,max` was the right one, so I computed the whole curve.

`simulations/test10_ceiling_exceedance_curve.py` imports the executed TEST-10 loader and quality cuts
verbatim, so the sample is identical. **Reproduction check first**: 106/153 = 69%, f_DM,max = 0.927 —
matches the 2026-07-15 run exactly.

### My own hypothesis, refuted

I expected `f_DM,max` to be a fragile single-galaxy extreme, and said so before running it.

```
Drop-the-max: B_req,max = 13.75 -> 2nd-largest = 13.14
              (only 4% of the bound rests on one galaxy)
```

The tail is smooth, not spiky. **The outlier objection is refuted and is recorded as refuted**, in the
proposal and here. Guessing the failure mode of a statistic is not the same as computing it.

### The real problem with it, which is structural rather than statistical

| ceiling definition | B_max | f_DM cap | exceeding |
|---|---|---|---|
| 1/Ω_m (site's choice) | 3.17 | 0.685 | 106/153 = **69%** |
| (Ω_m−Ω_b)/Ω_b | 5.39 | 0.814 | 44/153 = **29%** |
| Ω_m/Ω_b (baryon budget, most permissive) | 6.39 | 0.843 | 28/153 = **18%** |
| f_DM,max (the prescribed cite) | 13.75 | 0.927 | 0/153 = **0%** |

**`f_DM,max` is excluded by exactly one galaxy — by construction.** Nothing exceeds a sample maximum.
The prescribed headline is the single weakest-supported point on the entire curve, which is the
opposite of what "robust" was meant to buy. It swaps a convention-dependence for a K=1 support; it
moves where a referee pushes rather than removing the push.

### The statistic that is actually both definition-free and robust

Evaluate at the **most permissive** candidate ceiling — worst-case over conventions, full sample behind it:

```
K galaxies of support  ->  largest B_max excluded
    K =   1                 B_max <= 13.75     <- the prescribed cite
    K =   5                 B_max <=  9.46
    K =  10                 B_max <=  8.39
    K =  20                 B_max <=  7.08
    K =  40                 B_max <=  5.69
    K =  76                 B_max <=  4.13
```

Required-boost quantiles (bootstrap 95% CI, 10k resamples, seed 20260729): p50 = 4.08 [3.83, 4.56] ·
p90 = 7.59 [6.74, 8.39] · p95 = 8.84 [7.62, 9.51].

**The consequence favours the result.** The triage says that under Ω_m/Ω_b "the *median* passes" and
therefore "the kill stands on the tail, not the median." The first clause is right; the second
undersells it by an order of magnitude in support. **28 galaxies — 18% of the sample — exceed even the
most permissive candidate ceiling.** Not a tail of one.

> **Recommended headline, replacing both "69%" and f_DM,max: a bounded-boost modified-gravity class
> with ceiling B_max ≤ 6.4 is excluded by 28 SPARC galaxies; with B_max ≤ 8.4, by 10. The exclusion
> holds under every candidate normalisation, being evaluated at the most permissive.**

**Scope, stated plainly: this is not the registered sweep.** That protocol also recomputes TEST-09's
BTFR slope per ceiling definition under a pre-fixed verdict rule, and **remains unrun**. This
diagnostic does not discharge it and does not claim to.

---

## 4. What this does to the publication queue

**New: REC-2026-039 — Bounded-Boost Class Exclusion from SPARC Dwarf DM Fractions** (readiness 0.72,
HIGH, arXiv astro-ph.GA). The proposal called it "the cheapest item to ship: no A2ACW dependency, no
screening-literature prior-art exposure, one dataset, one figure." That holds, and its headline number
is now executed and reproducible from data already in this repo.

**Its gating weakness is named, not smoothed over.** The prior-art check has never been run. The
proposal's own open question #3 asks whether the B_max exclusion already exists in the literature under
another name and notes the program's rediscovery detector has a **4/4 catch rate** and has never been
pointed here. Two days ago that same class of gate found a published counterexample to transferable
null #1 inside a single walk. **Run it before drafting, not after.** Also unrun: an Υ* (M/L)
sensitivity — f_DM scales with the assumed mass-to-light ratio, a referee will ask, and it is cheap.

**Revised advisory order** (supersedes yesterday's): **#4 Bounded-Boost Class Exclusion (REC-039) →
REC-038 (drafted, reproducibility-verified, cs.AI) → #2 DESI mechanism-class → #1 Locality No-Go
(still blocked) → #3 A2ACW.** The new item leads the *physics* queue because it is the only one whose
headline number is executed, definition-robust, and free of an open blocking question. REC-038 remains
the shortest path to anything posted at all.

**Null #1's fork is still open.** Reading A vs Reading B (filed 07-28) got no new evidence today. It
still decides preprint #1's headline sentence.

### Status changes

| Rec | Change | Why |
|---|---|---|
| **REC-2026-039** | **NEW** — readiness 0.72, HIGH | Headline statistic executed today; gated on the unrun prior-art walk |
| **REC-2026-037** | framing addendum, `date_updated` 07-29 | The shippable half split out to REC-039; shouldn't wait on this one's packaging |
| REC-034 / 035 / 036 / 038 | held (0.97 / 0.95 / 0.60 / 0.93) | No new evidence |

`recommendations.json` diff: **34 insertions, 3 deletions**. Third consecutive pass with no
re-serialization churn — and this time verified rather than assumed: a round-trip check before editing
confirmed `ensure_ascii=True` + trailing newline reproduces the on-disk bytes exactly.

### Phase-0 widened scan (the 07-27 fix, exercised)

`Research/papers/` — one item, unchanged, already REC-038. `Research/preregistrations/` — one item,
unchanged. `Research/proposals/` — the 07-28 additions, handled above. **No new numbered session (S691
unchanged), no new complete arc.**

---

## 5. Phase 1 Synchronism: gates green, no section edit

| Check | Result |
|---|---|
| Claims freeze (`--check`) | ✅ green, 10 claims, before edits |
| Lone-CR regression | ✅ clean over `whitepaper/**`, `docs/**`, `Research/**` — 07-27 fix holding, 3rd pass |
| Terminology drift | zero |
| Section edits | **none — deliberate** (see §1) |

---

## 6. Phase 1 Web4: the same failure class, named out loud, in a different repo

Reviewed `web4@206dd00..f10d299` — nine commits, two in whitepaper scope, both handled web4-side. Zero
terminology drift (the diff's LCT / T3-V3 usages are all referential).

`97c9eb2` continues exactly the thread yesterday's `b2e2888` opened, and the web4 side names the class
itself: *"third instance in two passes — a declared property about a reference implementation accepted
without computing on the evidence, which is the thing §11 itself argues against."* Two §11 claims
failed against primary sources: production operation attributed to the public `web4/hub` tree (whose
own README calls it a reference proof-of-concept at 0.1.0-alpha.0 — and which contradicted §11's own
closing paragraph calling all three implementations "research-stage"), and a "crates.io, PyPI, and npm"
claim that reads as a 2×3 matrix but has two empty cells. What was true was kept and is stronger: the
hub *is* run daily, verified on a live host rather than from prose.

**And that class is mine today too.** "A nested submodel cannot win" is a **declared structural
property that fails when you compute on it**. `f_DM,max` is a **declared robust statistic that fails
when you compute on it**. Two repos, one class, same 48 hours — and unlike the earlier cross-repo
pattern I logged as "a surface reporting success for work that never ran," this one I found inside my
own lane rather than pointing at someone else's.

---

## 7. Housekeeping

- `private-context`'s 07-28 paused rebase **completed cleanly**. My log entry is on `main`. The rescue
  ref `publisher-log-2026-07-28` did its job and is now redundant; left in place (harmless, and
  deleting another track's recovery artifact is not mine to do).
- Pulled with `--autostash` in both repos, per the memory correction made yesterday. No conflicts.
- `AGENTS.md` / `CLAUDE.md` carry uncommitted GitNexus index-count updates in both repos — supervisor
  scope, left untouched per precedent `a13894da`.
- New file `simulations/test10_ceiling_exceedance_curve.py` committed with a header stating explicitly
  that it is **not** the registered sweep, so a future reader can't mistake it for one.

---

## Summary

The research lane delivered a capstone framing for the whole galaxy arc and a self-caveat on its own
headline number, both on 07-28. I checked them by computing on them rather than reading them, and both
moved.

The reframe's "a nested submodel cannot win, a priori" conflates fit with selection — ΛCDM is nested in
wCDM and wins anyway. What nesting buys is a dichotomy, not a verdict; which branch obtains is
empirical, so the SPARC runs were the decision procedure, not corollaries of a foregone conclusion. The
verdict is untouched: the ceiling binds. And the `f_DM,max` statistic prescribed to fix a real
convention-dependence turns out to be excluded by exactly one galaxy by construction — while the
statistic nobody had computed, exceedance at the most permissive ceiling, is backed by 28 galaxies and
makes the class-exclusion null shippable today.

**The pattern I keep landing on has inverted. This week's findings were all "a surface reported success
for work that never ran" — detection failures pointed outward. Today's two are the same act pointed at
my own lane: a structural property and a summary statistic, both declared, both accepted because they
sounded right, neither computed on. And the one hypothesis I did carry in — that f_DM,max would be a
fragile outlier — was refuted by four percent. Under-claiming has the same failure mode as
over-claiming. It just doesn't feel like one.**
