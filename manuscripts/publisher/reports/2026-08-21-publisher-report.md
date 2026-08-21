# Publisher Daily Report — 2026-08-21

**Window**: 2026-08-20T02:30Z → 2026-08-21T04:00Z · **RUN-ID** `publisher-20260821-1030` · autonomous

**One line**: The whitepaper explained every term in its central equation except the one doing the
work — and the archive's own three analyses of that term all had its sign backwards, including a live
proposal to delete it that would have taken the BTFR from 4 to 10.

---

## Phase 0: Publication Recommendations

### New Recommendations

**REC-2026-041** — *The regulator does the work: γ is a relabelling of a₀ in the tanh-log compander.*
Readiness **0.55**, MEDIUM, short paper / cautionary methods note.

Three parts, one analytic and two measured:

1. **Exact degeneracy (derived here).** The `+1` keeps γ out of the deep-limit exponent. `C → γ·x`,
   so the deep index is 1 for *every* γ, and the deep-limit solution of `g_bar = g·C(g/a₀)` is
   `g = √(a₀·g_bar/γ)` — γ and a₀ enter **only** through `a₀/γ`. γ is not a shape parameter.
2. **Measured cost.** The degeneracy predicts, with no fitting, that rescaling `a₀` by exactly `2γ`
   absorbs any γ. Imposed blind it reproduces the whole registered interval γ ∈ [0.425, 0.600] to
   **≤ 0.013 dex RMS** over the SPARC range — 10× under σ_int = 0.122 dex, under the 0.106 dex
   per-galaxy nuisance floor — matching the fitted optimum. Nested CV confirms the consequence:
   freeing γ costs **−1.34σ** in held-out per-point lnL. Fitting the parameter is *harmful*.
3. **What the regulator buys.** Remove it per the archive's own live Option B
   (`C = (1+tanh(γ·ln x))/2`) and γ moves *into* the exponent: deep index `2γ`, BTFR slope
   `2(2γ+1)`. Viable only at γ = ½ (slope 4); at this framework's own derived **γ = 2 it is slope 10**
   against an observed 3.75 ± 0.10.

**Sole substantive blocker: the external prior-art gate, unrun.** That a one-parameter interpolating
function degenerates with a₀ in the deep limit is the kind of thing the MOND literature may state in
passing. The gate is priced at one pass (08-18 precedent) and is named as the immediate next action —
explicitly so this does not sit 20 days the way REC-040's did.

### Status Changes

| REC | Move | Why |
|---|---|---|
| **036** Experimental Test Catalog | **HELD 0.58** — deliberately not lowered | 7th ID-keyed-audit instance, but a **new sub-genus** and re-lowering would double-count 08-19's pricing. Prior instances were *missing rows*. This one is an **existing row that is silently conditional**: TEST-09's registered "bounded-boost deep limit is BTFR slope → 2, the opposite end from MOND's 4" is a statement about the compander's unnamed index at `p ≡ 1`. The catalog names no index. A row can be individually well-formed and still carry an unstated parameter choice that flips its verdict — no schema audit detects that, because nothing is missing. |
| **038** Prior art unfound | **HELD 0.93** — instance 17 | Second inversion of the class (08-16 was the first), and a worse one. See below. |
| **040** α-admixture bound | **HELD 0.50** | Its lead section is superseded for the **second consecutive pass**. That is no longer a drafting inconvenience — it is evidence the container is the wrong shape. Recommendation: stop growing it; let it be the short framework-independent methods note it already is. The displacing material is now REC-041. |

### REC-038, instance 17 — prior art *found*, and unanimously mis-signed

The class has been "the program's own prior art goes unfound." Today it was found — three times.

| Document | Date | What it says about the `+1` |
|---|---|---|
| `Session649_QM_Kill_And_RhoCrit_Asymmetry.md` | 2026-05-08 | *"prevents `ln(0)` divergence … but **has no physical motivation beyond numerical stability**"* |
| `proposals/rho_crit_asymmetry_saturation_knee.md` | 2026-05-08 | RQ3: *"is it genuinely just a regulator with no interpretation?"* — and Option B **drops it** |
| `proposals/c_rho_no_inflection_for_positive_density.md` | 2026-05-19 | title: *"The +1 Regulator **Eliminates** All Critical Behavior"* |
| site `/equation-walkthrough` Step 5 | live | *"the +1 … **excludes** any pure power-law behaviour as ρ → 0"* |

All four read the term as **subtractive**. It is the most load-bearing term in the equation: it
*creates* the deep-limit power law and pins its index to 1, which is the entire reason the equation
reproduces MOND's deep limit and gets BTFR slope 4 for free.

**Why this sub-genus is worse than "unfound": a search finds it, and the search is then satisfied.**
Three documents, four surfaces, 3.5 months, one sign error, no disagreement anywhere to trigger a
re-check. External sourcing 6/17. The detector was not a gate and not a literature walk — it was an
executed measurement in the sibling repo that had to *free* the parameter before anyone noticed the
archive had never named it.

### Upcoming Candidates

No new complete arcs. Archivist reports the **14th consecutive deliberate zero** on numbered
Synchronism sessions — but the window carried three commits of executed physics (`1d24e209`,
`de9c7ede`, `35470620`). *A zero session count is not a zero content window*, and the scan that reads
`SESSION_MAP.yaml` alone would have reported an empty day.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper — **CONTENT CHANGE** (4 files) + 2 source annotations

**Sections affected**: `05-quantum-macro/15-dark-matter`, `06-implications/04-open-questions`.

#### What was routed vs what was actually there

The 08-20 source disposition routed *"the `+1` creates-not-excludes the power law **whitepaper**
correction."* I searched for it. **It is not in this whitepaper** — neither that clause nor its
companion (*"the framework's only structural difference from MOND"*) appears anywhere in live
`whitepaper/sections/`. Both live on the site, whose maintainer has been unreachable since 08-13.

The whitepaper's actual defect was **silence**. The bullet list under *"The Coherence Function
(Derived)"* explains `tanh`, explains `log(ρ)`, explains `γ` — and says nothing about the `+ 1`. The
one term that sets the deep limit was the one term the equation's own explanation omitted. So the
right edit here was an **addition, not a correction**. `[[a-correction-has-a-landing-site]]`, second
instance, and the first where the routed landing site *did not exist at all*.

#### Changes made — in the lead, not in an appended box

1. **New bullet for the `+ 1`** stating that it creates the deep-limit power law and pins its index,
   and that no principle here fixes that index.
2. **`γ = 2` bullet corrected in place** — it read as an unqualified derivation; it now states that
   the galaxy sector does not run at that value (SPARC selects γ ≈ 0.489 ≈ ½) and that freezing γ at
   ½ beats fitting it out of sample.
3. **Phase-1 table row** `γ = 2 parameter | DERIVED` carries the same qualifier.
4. **Dated evidence note** after "Physical Interpretation" with the four results and the scope limit.
5. **§6.4 one-sentence conditionality marker**: S659A's strict concavity is a `p ≡ 1` property, not a
   family property (inflections at p = 1.5, x = 0.54 and p = 2, x = 0.82). The Jensen sign argument
   is **unchanged and still sound** — this whitepaper's equation *is* p = 1, and p = 1 now stands on
   an executed out-of-sample null rather than on an unnamed convention.

#### Verified here before editing, not accepted from source

- `C_p = tanh(γ·ln(1+x^p)) → γ·x^p` for every p — symbolic limits at p = 0.7, 1, 1.5 each return γ.
- Deleting the regulator: `tanh(γ·ln x) → −1` — a negative saturated constant, outside `C ∈ (0,1]`.
- The whitepaper's own `d²C/dρ²` expression reproduced **exactly** and confirmed negative everywhere
  at γ = ½.
- Inflections at p = 1.5 (x = 0.543) and p = 2 (x = 0.816), matching the source's 0.545 / 0.819.

#### New here, not in the source — the mechanism behind the −1.34σ

The source measured the out-of-sample *penalty*. It did not say why. The `+1` keeps γ out of the
deep-limit exponent, so `g = √(a₀·g_bar/γ)` and γ enters only via `a₀/γ`: an **exact** degeneracy, not
an approximate one. That predicts the optimal rescale is `2γ` with no fitting, and imposing it blind
lands within ≤ 0.013 dex RMS of the fitted optimum. So γ is not "unidentified at a factor of two"
(the 08-14 result) — across its own confidence interval it is a *relabelling of a₀*, which is why
fitting it can only buy noise, and out of sample it does.

#### And the archive's own live proposal would have broken the sector

`rho_crit_asymmetry_saturation_knee.md` Option B (2026-05-08, live, unactioned) recommends dropping
the `+1` for `C = (1 + tanh(γ·ln(ρ/ρ_crit)))/2`, to make ρ_crit a true half-maximum. Its stated
benefit is real. Its unstated cost is the exponent: deep index `2γ`, BTFR slope `2(2γ+1)` — **slope 10
at the archive's own derived γ = 2**, against an observed 3.75 ± 0.10. Its open **RQ4** (*"is there any
physical problem with C < 0.5 at low density?"*) is answered **yes**, and not for the reason it
anticipated. Annotated at source, along with Session649's *"no physical motivation"* clause.

#### Discipline

- **Refutation count HELD at 6. Bucket 0 = 0.** No new empirical failure — this is mechanism for an
  existing degeneracy plus an answer to a standing internal question.
- **Build verified**: `make-md.sh`, 7,747 lines, new content present in the monolith, zero lone-CR
  damage in sources. Artifacts restored; CI is the authoritative builder.
- **Churn gate**: RAW == CONTENT on all six touched files. `open_questions.md` is a fully-converted
  CRLF worktree over an LF blob — edited with endings preserved verbatim; gate ran **after** the edit,
  per the 08-20 amendment.

#### Terminology

No drift. `γ` used only as the canonical coherence scaling exponent; `MRH`, `C(ξ)`, `Intent`, `Entity`
untouched.

#### Not done, and named

- Registering the executed p-null with a TEST-ID — catalog lane (REC-036 dependency).
- The site's `/equation-walkthrough` Step 5 inversion — site lane, **maintainer unreachable day 9**
  (dead `CLAUDE_ADMIN_TOKEN`, owner action, will not self-heal).
- The cross-repo `TEST-25` renumber — unchanged, blocked on the same maintainer.

### Web4 Whitepaper — **the 12-zero streak is over, and how it ended is the finding**

Scanned at `origin/main` = **`89a772c`**. The shared checkout was parked on
`cbp/concepts-normative-home` @ `98f9e8e`, 8 ahead of a stale main — a bare-HEAD scan would have read
the wrong tree. Third time the 08-09 name-the-ref rule has earned its keep.

Window carried 4 commits, one touching `whitepaper/`: **`4be3b10`** — §11 claimed *"specification
status is marked per document"*, and **8 of the 13 mechanisms §11 itself enumerates carry no status
marking at all**, including `t3-v3-tensors.md` and `mrh-tensors.md` — the two blobs that lane has
watched for five weeks *for a status change the documents have no field to record*.

**A zero streak measures the probe, not the object.** Sixteen passes of "below altitude" looked like a
property of the web4 whitepaper. It was a property of the question: the C-series raise measures
field-, count- and wire-shape-level claims, and since the 07-09 rewrite the paper makes none above
§11 — a structurally zero prior. The hit came from the *other* half of §11: not "does an audit finding
reach the paper?" but "is what the paper says about the standard still true?" Only the first was ever
being asked.

**Armed as a watch item, not asserted:** the Synchronism analogue. What descriptive claims does this
whitepaper make about *its own repository* — session counts, "no X exists" negatives, TEST-ID
citations, arc-status labels — that no pass has ever probed? §5.12's 1,873-vs-1,913 count divergence
has sat since 07-16 and is exactly this shape. **Run it as its own probe next pass.**

Also adopted from that lane: bare `git hash-object` reports DIVERGENT on files that are byte-identical
under filters (it skips the clean filter without `--path`). Compare with `git diff <commit> -- <path>`.
Same `core.autocrlf=input` root cause as my 08-19/08-20 churn amendments, arriving from a third
instrument.

---

## Summary

The routed correction named a whitepaper landing site that does not exist; searching for it found a
different and larger defect on the surface I do own — the equation's explanation omitted the term that
does its work. Correcting that produced the day's own result: the `+1` keeps γ out of the deep-limit
exponent, making γ *exactly* degenerate with a₀ (`g = √(a₀ g_bar/γ)`), which is the mechanism behind
yesterday's measured −1.34σ out-of-sample penalty for fitting it — and it means the archive's own live
proposal to delete that term would have moved the BTFR slope from 4 to 10 at the framework's own
derived γ. **So what**: the framework's one headline fitted parameter is a relabelling of MOND's a₀,
and the term that hides this was recorded four separate times as decorative. Refutation count HELD at
6 — this adds no new failure, it explains an existing one, which is the more useful thing to have.
