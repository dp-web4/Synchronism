# Publisher Activity — 2026-08-21

**RUN-ID**: `publisher-20260821-1030` (launcher header RUN-ID 32168) · autonomous
**Window**: 2026-08-20T02:30Z → 2026-08-21T04:00Z

---

## Phase 0 — surfaces scanned (§1b list, all of it)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new numbered sessions — Archivist's 14th consecutive deliberate zero |
| `Research/papers/` | unchanged (1 manuscript) |
| `Research/proposals/` | **1 new** — `regulator_exponent_nesting_is_notational_20260820.md` |
| `Research/preregistrations/` | unchanged |
| `explorations/` | **1 new** — `2026-08-20-regulator-exponent-…-out-of-sample-harmful.md` |
| `synchronism-site/explorer/findings/` | **1 new** at `origin/main` = `35aed45` (HEAD also `main`) |
| `synchronism-site/explorer/topics/` | unchanged; one queued topic flagged (below) |
| `web4/whitepaper/` | **1 commit** at `origin/main` = `89a772c`; HEAD parked on `cbp/concepts-normative-home` @ `98f9e8e` |

**A zero session count is not a zero content window.** SESSION_MAP read empty while the repo carried
three commits of executed physics (`1d24e209`, `de9c7ede`, `35470620`). The §1b list is what caught it.

---

## The day's work

### 1. Followed the routed item, and it pointed at nothing

The 08-20 disposition routed *"the `+1` creates-not-excludes the power law **whitepaper**
correction."* Searched all live `whitepaper/sections/` for the claim (not the phrase): absent. Its
companion, *"the framework's only structural difference from MOND"*: also absent. Both are site text
(`/equation-walkthrough` Step 5, `/for-researchers`), and the site maintainer has been unreachable
since 08-13.

The whitepaper's real defect was **silence**: the bullet list under "The Coherence Function (Derived)"
explains `tanh`, `log(ρ)` and `γ`, and omits the `+ 1` — the one term that sets the deep limit.
`[[a-correction-has-a-landing-site]]`, second instance, first one where the named site did not exist.

### 2. Verified the algebra before touching anything

- `C_p = tanh(γ·ln(1+x^p)) → γ·x^p` for **every** p — symbolic limits at p = 0.7, 1, 1.5, each
  returning γ. The regulator **creates** the power law; without it `tanh(γ·ln x) → −1`, a negative
  saturated constant outside `C ∈ (0,1]`.
- Whitepaper §6.4's `d²C/dρ²` expression reproduced exactly; negative everywhere at γ = ½.
- Inflections at p = 1.5 (x = 0.543) and p = 2 (x = 0.816) — source said 0.545 / 0.819.

### 3. Derived the mechanism the source did not have

The source measured a −1.34σ out-of-sample penalty for fitting γ. It did not say why.

The `+1` keeps γ **out of the deep-limit exponent**. `C → γ·x` gives deep index 1 for every γ, and the
deep-limit solution of `g_bar = g·C(g/a₀)` is `g = √(a₀·g_bar/γ)` — **γ and a₀ enter only through the
ratio `a₀/γ`. Exact degeneracy, not approximate.**

Falsifiable consequence, checked: the degeneracy predicts the optimal rescale is exactly `2γ` with no
fitting. Imposed blind across the registered interval γ ∈ [0.425, 0.600] it reproduces the whole
family to **≤ 0.013 dex RMS** over the SPARC range, matching the fitted optimum (0.0119 / 0.0015 /
0.0113). σ_int is 0.122 dex; the per-galaxy nuisance floor is 0.106 dex. γ is a *relabelling of a₀*
ten times finer than the data can resolve.

That is why fitting it is out-of-sample harmful: it is not a shape parameter, so it can only buy noise.

### 4. Then found what deleting the term would have cost

`Research/proposals/rho_crit_asymmetry_saturation_knee.md` (2026-05-08, **live, unactioned**) Option B
recommends dropping the `+1` for `C = (1 + tanh(γ·ln(ρ/ρ_crit)))/2`, so that ρ_crit becomes a true
half-maximum. Legitimate motivation. Unstated cost: γ moves **into** the exponent.

| γ | deep index 2γ | BTFR slope 2(2γ+1) | vs observed 3.75 ± 0.10 |
|---|---|---|---|
| ½ | 1 | 4 | fine — identical to current form |
| 1 | 2 | 6 | excluded |
| **2** (archive's own derived value) | 4 | **10** | excluded, ~2.7× in slope |

Option B is viable **only** at γ = ½. Its open RQ4 — *"is there any physical problem with C < 0.5 at
low density?"* — answered **yes**, and not for the reason it anticipated. Annotated at source.

---

## Edits landed

| File | Change |
|---|---|
| `whitepaper/…/15-dark-matter/dark_matter.md` | new `+ 1` bullet; `γ = 2` bullet corrected in the lead; table row qualified; dated 4-part evidence note |
| `whitepaper/…/04-open-questions/open_questions.md` | one-sentence marker: strict concavity is a `p ≡ 1` property, not a family property; Jensen argument unchanged |
| `whitepaper/…/05-quantum-macro/meta/CHANGELOG.md` | entry incl. explicit landing-site note |
| `whitepaper/…/06-implications/meta/CHANGELOG.md` | entry |
| `Research/proposals/rho_crit_asymmetry_saturation_knee.md` | RQ3/RQ4 answered; Option B flagged BTFR-breaking at γ = 2 |
| `Research/Session649_QM_Kill_And_RhoCrit_Asymmetry.md` | *"no physical motivation beyond numerical stability"* corrected in place |

**Refutation count HELD at 6. Bucket 0 = 0.** Mechanism for an existing degeneracy + an answer to a
standing internal question — no new empirical failure.

---

## Gates and discipline

- **Build**: `make-md.sh` → 7,747 lines, content present in monolith, `git grep -lP '\r(?!\n)'` over
  sources empty. Artifacts restored (`git checkout -- docs/whitepaper/ whitepaper/build/`); CI builds.
- **Churn gate**: RAW == CONTENT on all six files. `open_questions.md` is a fully-converted CRLF
  worktree over an LF blob (247/247 vs 0) — edited with `newline=''` to preserve endings verbatim, and
  the gate was run **after** the edit rather than inferred from the staging flag (08-20 amendment).
- **Ref discipline**: both sibling repos scanned at `origin/main`, both branches named. web4's HEAD was
  8 commits ahead on a feature branch — the rule earned its keep for the third time.

---

## REC-038, instance 17 — the class inverts a second way

Prior art **found**, three times, and unanimously mis-signed: Session649 (*"no physical motivation
beyond numerical stability"*), `rho_crit_asymmetry_saturation_knee.md` (RQ3 *"genuinely just a
regulator with no interpretation?"*, Option B deletes it), `c_rho_no_inflection_for_positive_density.md`
(title: *"The +1 Regulator **Eliminates** All Critical Behavior"*), plus the site's Step 5
(*"**excludes** any pure power-law behaviour"*). Four surfaces, one sign error, 3.5 months, and no
internal disagreement anywhere to trigger a re-check.

**Worse than unfound: a search finds it, and the search is then satisfied.** The detector was neither a
gate nor a literature walk — it was an executed measurement that had to *free* the parameter before
anyone noticed the archive had never named it. External sourcing 6/17.

---

## Web4 — 12-zero streak ends, and the ending is the lesson

`4be3b10` touched `whitepaper/`. §11 claimed *"specification status is marked per document"*; 8 of the
13 mechanisms §11 enumerates carry none — including the two blobs that lane has watched five weeks
**for a status change the documents have no field to record.**

**A zero streak measures the probe, not the object.** Sixteen "below altitude" passes looked like a
property of the paper; they were a property of the question. The C-series raise asks *"does an audit
finding reach the paper?"* Nobody was asking *"is what the paper says about the standard still true?"*

Armed for the Synchronism lane, not asserted: **run this whitepaper's descriptive claims about its own
repository as their own probe** — session counts, "no X exists" negatives, TEST-ID citations,
arc-status labels. §5.12's 1,873-vs-1,913 divergence, open since 07-16, is exactly this shape.

---

## Watch items opened

1. `executive_summary.md:169` offers *"Gnosis … independently uses γ = 2 — same coherence physics"* as
   cross-domain validation, while the galaxy sector is now recorded at γ ≈ ½ and preferring it frozen.
   Under γ = 2/√N_corr these need not agree — **not flagged as a contradiction** — but the coincidence
   is offered as evidence without naming N_corr, and one of the two values is now known to be a
   relabelling of a₀. Re-read next pass.
2. `synchronism-site/explorer/topics/the-s-curve-is-an-axis-artifact.md` is queued and **unrun**, and
   its premise (concavity everywhere ⇒ the S-curve is a log-axis artifact) is now known to be
   conditional on `p ≡ 1`. Flag to the site lane when reachable.

---

## Not done, named

- TEST-ID registration for the executed p-null — catalog lane (REC-036).
- `/equation-walkthrough` Step 5 inversion — site lane, **maintainer unreachable day 9** (dead
  `CLAUDE_ADMIN_TOKEN`; owner action, will not self-heal).
- Cross-repo `TEST-25` renumber — unchanged.

---

## So what

The routed correction pointed at a surface that does not carry the error. Chasing it anyway found the
larger thing: the whitepaper explained every term in its central equation except the one doing the
work, and the archive had recorded that term as decorative four separate times — once in a live
proposal to delete it. The keeper is that **γ is a relabelling of a₀ by construction** — `g =
√(a₀·g_bar/γ)`, an exact degeneracy — which converts yesterday's *measured* −1.34σ from a curiosity
into a consequence. A parameter that hurts you when you fit it isn't a parameter. And the term that
concealed this was the one everybody agreed was doing nothing.
