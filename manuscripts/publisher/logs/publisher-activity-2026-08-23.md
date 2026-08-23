# Publisher Activity — 2026-08-23

**RUN-ID**: 26069 · autonomous · window 2026-08-22T09:30Z → 2026-08-23T10:30Z
**Refs scanned**: `Synchronism` origin/main `d2ccb960` (HEAD main, same) · `synchronism-site`
origin/main `ccacc2d` (HEAD main, same) · `web4` origin/main `dbfd2747` (**HEAD parked on
`cbp/concepts-normative-home`** — bare-HEAD scan would again have read the wrong tree)
**Refutation count**: HELD at 6 · **Bucket 0**: 0

---

## The day in one sentence

A rule this lane put into the whitepaper on 08-20 was refuted on 08-22 by the exact falsifier the
whitepaper itself had named, and the refutation arrived asserting three times that nothing had
propagated — while the refuted text sat in three body sections and ten built files on the deployed
Pages surface.

---

## 1. What fired, and why it was mine to fix

On **2026-08-20** this lane propagated an 08-19 site-explorer result into
`executive_summary.md`, `05-quantum-macro/15-dark-matter/dark_matter.md` (§5.15) and
`conclusion.md`:

> the required non-locality is specifically **inward-cumulative**; a *symmetric* kernel `f(|r−r′|)`
> fails at **any** range; **the discriminating axis is symmetry, not range**; therefore every
> natural repair of a pointwise multiplier lies in the dead branch and the live branch reconstructs
> `g_bar` itself.

It was propagated **with its hedge**, and the hedge was unusually good — it named the falsifier and
pre-registered the consequence of failing it:

> its own cheapest falsifier — project a genuine 3-D Yukawa `e^(−mr)/r` through an exponential disk
> onto a radial kernel and re-run — is **unexecuted**, and **if a screened *linear* scalar lands in
> the live branch the sorting rule is wrong and the escape taxonomy reopens**.

On **2026-08-22** (`298a950e`, site script `yukawa_symmetric_kernel_self_check.py`) that falsifier
ran and returned exactly that branch. So the whitepaper's own conditional retraction fired. It does
not self-apply; walking it was this pass's obligation.

**Named unrun for three consecutive passes.** 08-20 recorded "Yukawa first, then the prior-art
walk"; 08-21 and 08-22 both logged it unexecuted. Another lane ran it. That is fine — but it means
the item I was tracking as *my* blocker was discharged by someone else while I tracked it, which is
worth noticing about how blockers get named here.

---

## 2. Verified before propagating — and the analytic half is the stronger half

The refutation needs **no SPARC number**. `(∇²−m²)h = 4πGρ` has Green's function `e^(−mr)/r`, which
is two-sided and isotropic — *symmetric* — at **every** screening length, and recovers Poisson as
`m → 0`. Re-derived in closed form in this lane (`sympy`):

| `λ = 1/m` (units of `r`) | `F_Yukawa / F_Newton` |
|---|---|
| 1 | 0.7358 |
| 10 | 0.9953 |
| 10³ | 1.0000 |

matching the source's 0.736 / 0.995 / 1.000 exactly. **The branch the rule declared dead contains
`g_bar` itself**, so the rule was refutable by construction on the day it was written — no execution
required. That is the form I put in the whitepaper, because it does not depend on anyone's pipeline.

Empirically it also fails, and I re-ran the source's own bootstrap rather than accept the report:
`yukawa_addendum_bootstraps.py` → **exit 0**, tables reproduced value-for-value.

- validation gate: unscreened Yukawa **0.1080 dex** vs `g_bar` **0.1107** (cost −0.0026)
- symmetric field `g_Y` **OVERLAPS** `g_bar` at **every** `λ_s ≥ 1 R_d`
- head-to-head at matched range: symmetric **beats** causal on every row (1.23/1.19/1.14/1.09/1.00/0.96/0.95/0.98 vs 1.31/1.27/1.20/1.13/1.07/1.04/1.02/1.03)
- 08-19's "separated at every finite λ ≤ 4 R_d" **does not reproduce**: +0.0075 [−0.0040, +0.0223]
  masked, +0.0096 [−0.0046, +0.0215] on 08-19's own unmasked point set — **OVERLAPS both ways**
  (its CI lower edge had been +0.0003 dex)

---

## 3. Two limits I carried that the source's triage note does not

Not relayed — read off the source's own printed bootstrap and re-derived here.

### (a) The replacement axis is range-conditional. It is not yet a ladder.

The source replaces "symmetry" with *which derivative of Φ keys the modification* — prior art
(Joyce, Jain, Khoury & Trodden, Phys. Rep. **568**, 2015), which is the right move and is the
genuine contribution (the *measurement*, not the classification). But the ladder is tabulated as if
range-independent, and in its own §(4) bootstrap **two of the three rungs swap verdict** between the
only two ranges tested:

| functional | `λ_s = 4 R_d` | `λ_s = ∞` |
|---|---|---|
| `⟨Σ⟩_Y` (normalised; the `∇²Φ` side) | **SEPARATED** +0.0190 [+0.0046, +0.0330] | OVERLAPS +0.0123 [−0.0094, +0.0415] |
| `\|Φ_Y\|` (the `Φ` rung) | **OVERLAPS** +0.0171 [−0.0023, +0.0376] | SEPARATED +0.0327 [+0.0118, +0.0578] |
| `g_Y` (the `∇Φ` rung) | OVERLAPS −0.0011 | OVERLAPS −0.0016 |

They cross. And the published ladder's `∇²Φ` row (1.40×) is **pointwise Σ, scored at `λ_s = 0`**,
while its other two rows are scored at `λ_s = ∞` — so the three rungs are not a matched comparison
at all. **Only `∇Φ` is range-robust.**

This is the same defect the document diagnoses — *conclusion wider than its own test* — committed
again in the replacement, one section later. I did not say so in the whitepaper in those words; I
just wrote the claim at the width its evidence supports.

### (b) The correction shed its own caveat in one hop

`Research/proposals/yukawa_selfcheck_...md` carries, in its own words:

> **My own going-in hypothesis (normalisation is the operative factor) is NOT established.** The
> normalised-mean member at `λ_s = ∞` is not resolved from `g_bar` (+0.0123 [−0.0094, +0.0415]) …
> The ladder claim rests on the **Φ-vs-∇Φ** contrast, which is bootstrap-clean, not on a
> normalisation axis.

The triage note this repository committed (`explorations/2026-08-22-yukawa-…md`) and the commit
message both state the normalisation mechanism **flatly, as the explanation**: *"a NORMALISED
smoothing of Σ (intensive, degrades with range)"*. Same author, same day, one hop — and the caveat
did not travel. `[[a-gate-is-an-event-not-a-property]]`, second clause, cleanest instance yet: the
conclusion propagated and the disclaimer stayed behind.

(It is also checkable against the numbers: the normalised `⟨Σ⟩_Y` **improves** monotonically with
range, 1.39× → 1.12×. It does not "degrade with range.")

---

## 4. REC-2026-038, instance 19 — a retraction that under-scopes its own blast radius

The finding says, twice, and the commit message a third time:

- *"The hedge worked — the rule never propagated as cited."*
- *"so it never propagated."*
- *"The failure was three days from a public page."*

All three are **false**, and one grep settles it. At the moment those words were written the refuted
rule was live in:

- 3 whitepaper body sections (`executive_summary.md:191`, `dark_matter.md:92`, `conclusion.md:79`)
- 10 built files under `docs/whitepaper/**` and `whitepaper/build/**` — **the deployed GitHub Pages
  surface** (`Synchronism_Whitepaper_Complete.md`, `section_0/5/7/11.html`, and the `web-clean`
  mirrors)

It was **on** a public page for three days, not three days from one.

**The new sub-genus**: the author checked the surface *they authored* (site findings) and inferred
the blast radius from **the hedge they had attached** rather than from a search. That is a different
mechanism from the sixteen instances before it, and it is more insidious, because the hedge is
evidence of *care* — it reads as though containment had been arranged. It had not.

The correct lesson is one word off the one recorded: **a hedge makes a correction cheap; it does not
make it unnecessary.** A hedged wrong claim on a public page is still a wrong claim. Here the hedge
did earn its keep — it named the falsifier, pre-registered the consequence, and made today's edit a
five-minute decision instead of a re-litigation. That is a large benefit. It is not containment.

Now a standing check (watch item): **for any incoming retraction, grep `whitepaper/sections/` AND
`docs/whitepaper/` for the retracted claim before accepting the source's blast-radius statement.**
Cost: one grep. Two lanes have now made this inference in eight days.

---

## 5. REC-2026-038, instance 20 — prior art in the research layer, absent from the publication layer

Fell out of the prior-art gate (§7). Counted with `git grep -il` over this repo vs the built
monolith:

| cited author | whitepaper mentions | repo files |
|---|---|---|
| Milgrom | 17 | 128 |
| Chae | 29 | 160 |
| McGaugh | 15 | 240 |
| Lelli | 14 | 515 |
| Burrage | 14 | 33 |
| **Famaey** | **0** | **26** (+45 in `synchronism-site`) |
| **Sanders** | **0** | 16 |

The whitepaper cites the sources it **argued with** and none of the **reviews that would situate the
argument**. Famaey & McGaugh 2012 — the canonical MOND Living Review, the document that would tell a
reader what is already known about μ-function families and `a₀` — is used 71 times across two repos
and appears **zero** times in the public paper. Previous sub-genus (08-20) was sector-wise: a
citation in-house in one sector, absent from the adjacent one. This one is **layer-wise**: research →
publication.

Today's edit fixes one instance of this in passing (Joyce+2015 was 0 in the whitepaper; it is now
cited in all three corrected sections).

---

## 6. New probe run — clean negative, and the cleanliness is the finding

Borrowed same-window from web4 `62f41005`, whose framing is the whole point:

> *"Eighteen Publisher passes have asked whether the corpus's news reaches the paper. None asked
> whether what the paper already prints still resolves."*

That lane found three dead citations in its own whitepaper. Run against Synchronism:

- **19 unique URLs, 19 resolve.** 17 are self-references, checked against `git ls-tree -r
  origin/main` (**0 dead** — exact, no network); the other 2 return **200**.
- **8 arXiv identifiers, 8 resolve**, all to topically-correct papers (MOND `a₀` evolution, DESI
  full-shape, MUSE-DARK RAR, the RAR identifiability audit).

**But the probe cannot fail here**, and that is the result: **0 of the 19 URLs leave the `dp-web4`
org. 0 DOIs. No References section in 7,802 lines.** The ~40 author-year literature citations are
non-resolvable prose. So the clean sweep measures the absence of an external link surface, not the
maintenance of one — `[[a-zero-streak-measures-the-probe]]` applied to a single probe rather than a
streak. Logged as a whitepaper-wide publication-readiness defect; **not** scored under REC-036,
which is a catalog container.

**Two instrument bugs caught, both of which would have produced spectacular false positives:**

1. `http://export.arxiv.org/...` returned **0 bytes silently** for all 8 IDs. Read naively that is
   *"all eight arXiv citations are fabricated"* — in a repository that carries an actual
   fabricated-consensus retraction (2026-07-10). `https://` works. I noticed only because
   `2206.04333` is a 2022 ID that is obviously real.
2. Three arXiv phrase queries returned 0 hits because of a **truncated word inside quotes**
   (`abs:"degenerac"`); arXiv does not prefix-match inside a phrase.

New rule, added to watch items: **a probe whose result is uniformly extreme is the suspect —
validate it against a known-positive before reporting.** This is `[[a-check-that-contradicts-a-proof-is-the-suspect]]`
generalised from proofs to instruments.

---

## 7. The prior-art gate — FIRST EXECUTION in five passes

Priced at "one pass" since 08-18 and unrun ever since, across REC-040 and REC-041. It has never been
the most urgent thing on any given day, which is precisely why it never ran — the correction work is
always more pressing. That is the efficiency attractor, not an accident, and holding it a sixth time
would have been perseveration. Partially run today.

**REC-2026-041** (γ is a relabelling of `a₀` in the tanh-log compander):

- **Nearest neighbour found, already in-house**: arXiv **2608.08945**, *"An Identifiability Audit of
  One-Parameter Structural Corrections to the Radial Acceleration Relation in SPARC"* (2026-08-09).
  It is already cited by this program — **it is the source of the 0.106 dex per-galaxy nuisance
  floor that REC-041's own result quotes**. Related but **distinct**: it establishes *empirical*
  non-identifiability (zero-point/nuisance freedom absorbs the gain); REC-041 derives an *exact
  analytic* degeneracy (`a₀` and γ enter only as `a₀/γ`). REC-041 is plausibly the analytic
  explanation for one class of that paper's empirical finding, and must be framed **against** it,
  not beside it.
- **The general claim is NOT cleared, and I am not recording a clean negative.** arXiv abstract
  search is the wrong instrument — the μ-function-shape/`a₀` degeneracy lives in paper bodies and
  reviews. Three targeted queries returned 0–1 hits, one of them from the bug in §6.
- **Named next step, now specific rather than generic**: Famaey & McGaugh 2012
  (arXiv:**1112.3960**), μ-function families section — does the deep-limit normalisation convention
  `μ(x) → x` already make *any* deep-limit prefactor degenerate with `a₀` by construction? If yes —
  the likely outcome — REC-041 narrows from "new degeneracy" to "this framework's γ is an instance
  of a known-degenerate normalisation." Smaller, still publishable, and exactly the
  degenerate-vs-productive reparametrisation question the project charter poses.

---

## 8. Changes made

**Whitepaper — 3 body sections, in place** (`executive_summary.md`, `05-quantum-macro/15-dark-matter/dark_matter.md`,
`conclusion.md`). Withdrawn: *"symmetric kernels fail at any range"*, *"the discriminating axis is
symmetry, not range"*, *"every natural repair … lies in the dead branch"*, *"the live branch
reconstructs precisely the variable `C(ρ)` was introduced to replace"*, and 08-19's *"separated at
every finite λ ≤ 4 R_d"*. Added: the analytic refutation, the validation gate, the non-reproduction,
the constructive bound `λ_s ≳ 1 R_d ≈ 2.4 kpc`, and the replacement axis **range-conditional with
both limits from §3**. The 73.0%-raw over-refutation guard preserved verbatim.

**Whitepaper — 3 CHANGELOG entries added**, and **3 historical entries marked `[SUPERSEDED]`**
(00-exec :112, 07-conclusion :100, 05-quantum-macro :79 partial — item (1) Brans–Dicke/Cassini
stands, only item (2) Kernel symmetry is superseded). Needed because the build monolith **includes
CHANGELOGs**, so those entries were still asserting the refuted rule as fact to a reader of the
public document.

**Verification**: `make-md.sh` exit 0, `make-web-clean.sh` exit 0. Line-level check of the rebuilt
monolith: **all 11 occurrences of the refuted phrase sit inside a refutation or supersession marker;
0 unguarded.** Artifacts **restored, not staged** — CI is the authoritative builder.

**Churn gate — fired once and was caught.** First write of
`05-quantum-macro/meta/CHANGELOG.md` gave **RAW 166 / CONTENT 4**: the file is **mixed-ending
(81 CRLF of 136 lines)** and a text-mode rewrite normalised it. Restored from the blob, re-appended
byte-wise with matching CRLF. Final state **RAW == CONTENT on every touched file** (20 ins / 8 del).

**State**: `recommendations.json` (4 records), `whitepaper_sync.json`, `whitepaper_web4.json`.

---

## 9. Recommendation decisions

| REC | before | after | why |
|---|---|---|---|
| **040** — form-free admixture bound | 0.50 | **0.58 ▲** | The 08-20 second blocker is **discharged by a refutation, not an addition**: the result that had superseded its Section 6 is the one that died, so Section 6 is the lead again. The same execution supplies a *citable* companion it lacked — the two-sided reading (pointwise admixture bounded **and** the smoothing length at which the bound dissolves, `λ_s ≳ 1 R_d`, measured). Prior-art gate is sole blocker again. Its `weaknesses` field carried the refuted rule verbatim and was corrected — my own state file was a landing site. |
| **041** — γ↔`a₀` degeneracy | 0.55 | **0.55 =** | Gate no longer wholly unrun (§7). Held rather than moved: the nearest neighbour found is in-house and non-fatal, but the general claim is uncleared and the named next step is likely to *narrow* the contribution. Moving on a half-run gate would be the mistake the gate exists to prevent. |
| **038** — prior art goes unfound | 0.93 | **0.93 =** | Two new instances (19, 20), one opening a sub-genus. External sourcing 6/20. |
| **036** — test catalog | 0.45 | **0.45 =** | No new evidence; 08-22 reconciliation-table restatement stands. §6's readiness defect is whitepaper-wide, explicitly **not** scored here. |

No new recommendations opened. No arcs closed; **16th consecutive deliberate zero** on numbered
sessions (Archivist), fourth consecutive pass where a zero session count sat over a non-zero content
window — and the **first** in which the window's content refuted live whitepaper text.

---

## 10. Web4

7 commits at `origin/main` `dbfd2747`; **1 touches `whitepaper/`** (`62f41005`, §6). No action
required from this lane — fixed there, artifacts rebuilt, the one author call routed to dp in their
log. The 12-zero streak stays ended; two consecutive non-zero windows.

Noted, not acted on: `fb63c51b` *"forum: add epistemic type and provenance vocabulary"* introduces
vocabulary in a repo with protected canonical terms. Nothing in it collides with
LCT/MRH/T3/V3/ATP/ADP/R6/R7 this pass, but vocabulary commits are the shape that produces drift —
re-read next window.

The 08-22 handback (`WITNESS_PURPOSE = "witnessing"`; operational-key delegation absent from
`web4-standard/` entirely) is unanswered there. No new keys appear to have been minted against it.

---

## Watch items

1. **Self-descriptive probe — class 3 of 4 run, clean** (§6). Probe stands at **2 hits / 3
   outings**. Remaining unrun: (a) negative-existence claims the whitepaper makes about its own
   repository; (b) **arc-status labels** — *"All prior research arcs closed as of Session #616"*
   (`conclusion.md:81`) against `SESSION_MAP`'s `complete_arcs: 43` and the post-closure sub-arc
   extensions the same sections describe running to #666.
2. **Hedge-vs-containment** (§4). Standing grep before accepting any retraction's blast-radius
   claim. Two lanes in eight days.
3. **Instrument-returns-uniform-null** (§6). Validate against a known-positive first.
4. Carried: five values for one session count (Archivist-owned, referred);
   `executive_summary.md:169`'s Gnosis γ = 2 cross-domain-validation clause;
   `the-s-curve-is-an-axis-artifact.md`, queued on a premise now known to be `p ≡ 1`-conditional.

---

## Not done, named

- **REC-041's prior-art gate is half-run.** Famaey & McGaugh 2012 §μ-function-families is the named,
  specific next action — no longer "do a prior-art walk."
- **REC-040's prior-art gate: still unrun**, now sole blocker again, 4th consecutive pass. Priced at
  one pass since 08-18 and that price is now demonstrated (one lane, one document, one day).
- **No References section / 0 DOIs in the whitepaper** (§6). Mechanical, cheap, and a hard blocker
  for any preprint route. Not started — it is a structural edit and wants a decision on citation
  style first, which is dp's call.
- Site lane: `/mond-unification` satellite closure text and the `~950×` headline at source.
  Maintainer unreachable **day 11** (dead `CLAUDE_ADMIN_TOKEN`, owner action, will not self-heal).
  The site also now carries the refuted symmetry rule at source and I cannot reach it.
- Cross-repo `TEST-NN` renumber — declined again (two disjoint namespaces, maintainer unreachable).

---

## So what

The whitepaper contained a **conditional retraction that had already specified its own trigger**:
name the falsifier, name what failing it would mean, mark the claim uncitable until run. The
falsifier ran, the condition fired, and walking it took one pass with no judgement calls left over.
That is the hedge working, and it is worth more than the claim it was protecting.

What the hedge did **not** do is what its author thought it did. "I marked it uncitable, so it never
propagated" is an inference about the world drawn from an intention, and it was wrong by three body
sections and ten deployed files. The gap between *I flagged it* and *therefore nobody relied on it*
is the same gap as `[[flagging-is-not-gating]]`, one level up: not flagging-instead-of-computing, but
**flagging-instead-of-checking-who-already-read-it**.

And the replacement axis, adopted the same day to fix a conclusion that was wider than its test, is
itself stated wider than its test — a two-range bootstrap in which two of three rungs cross, tabled
as a range-independent ladder, with one rung scored at a different screening length from the other
two. The failure class does not get retired by being named. It gets retired by someone re-reading
the new claim against the same standard, which is a different act from writing the lesson down.

Cheapest thing learned today, and it cost nothing: web4 asked whether its paper's own citations still
resolve, and found three dead. Synchronism's all resolve — because 0 of 19 leave the org, there are
no DOIs, and there is no References section. A probe that cannot fail is not a passing grade.
