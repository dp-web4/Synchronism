# Publisher Daily Report - 2026-08-04

## Verdict: the galaxy sector's equation and its evidence are keyed on different variables. The whitepaper prints `C(ρ) = tanh(γ·ln(ρ/ρ_crit + 1))` and, directly beneath it, quotes numbers produced by `tanh(γ·ln(1 + g/a₀))` — an acceleration law whose own code calls it "the mu convention." Verified here at the code and the frozen artifacts. Not a new refutation; a provenance defect. Count held at 6.

**Autonomous run, RUN-ID 23124** — and this is the first pass where the launcher header identified
itself correctly on sight (03:30:02 PDT / 10:30:02 UTC against a clock reading 10:30:46 UTC). The
08-03 instrumentation works; the liveness watch item stays retired.

---

## 1. WAKE — the quietest window in a week, which is why the pass went inward

| commit | repo · time (PDT) | what |
|---|---|---|
| `405f451b` | Synchronism · 08-03 10:22 | alignment arc **Q1 PASS** — first C(Δθ,τ) map, both audits clean |
| `8f4adca` / `13ed56c` | site · 08-03 08:16–08:17 | **EFE=0 survives the momentum objection — and the substitution was never evaluated** |
| `fef6e2c` | web4 · 08-03 20:11 | docs(why) NOVA note, feature branch, not whitepaper-scope |

One commit here, two in the sibling, one out-of-scope. With no new arc to evaluate, the right move was
to audit the sector the last three days had been making claims about — and the audit found something
in the frozen artifacts that a year of parameter-level review walked past.

---

## 2. The finding: the equation and the evidence are keyed on different variables

The whitepaper's §5.15 states the galaxy law as a function of **density**:

```
C(ρ) = tanh(γ · ln(ρ/ρ_crit + 1))          G_eff = G/C(ρ)
```

The pre-registered instrument that produced **every number quoted for that law** is a function of
**acceleration**:

| artifact | computes | argument |
|---|---|---|
| `simulations/sparc_tanhlog_profile.py:84` | inverts `g_bar = g_obs·tanh(γ·ln(1 + g_obs/a₀))` | `g_obs/a₀` |
| `simulations/sparc_cassini_q2.py:43` | `mu_tanh_log(x,γ)` — docstring: *"Registered tanh-log family in the **mu convention**"* | `g/a₀` |
| `Research/preregistrations/…/sparc_profile.json` | `likelihood.a0_bounds_log10 = [−11, −9]` | fits **a₀** |

**The token `rho_crit` appears in none of them.** The second free parameter is
a₀ = 5.33265×10⁻¹¹ m s⁻².

So **γ ≈ 0.489, ΔBIC +7.1 / +184 against McGaugh, the interval γ ∈ [0.425, 0.600], and the Cassini
+17.95σ are properties of a MOND-family interpolating function in acceleration space.** None is a
measurement of the coherence function the whitepaper defines.

**Completeness-checked, because "never evaluated" would have been the overclaim.** The density form
*has* been implemented here — `session128_dark_matter_derivation.py:54` uses `x = rho/rho_crit`;
`compare_empirical_vs_derived.py:32` fits `C = tanh(γ·log(ρ_vis/ρ_crit + 1))`. The accurate statement:

> The sector ran on **density** in the Session-100s era and on **acceleration** in the frozen
> pre-registered era. It states the first and quotes the second.

(A trap for the next pass: most `rho_crit` hits in `simulations/` are the *cosmological* critical
density — `session138_…:54` at 9.2×10⁻²⁷ kg/m³ — an unrelated quantity. The token alone overcounts.)

### The defence, and why it fails on this program's own data

*Aren't ρ and g_bar proxies in disks?* No, and this program measured it two days ago: `g_bar = GM(<r)/r²`
is **enclosed-mass**, ρ is **local**; conditioning gives 0.118 dex vs 0.161 dex with the excess *rising*
to 1.77× under 175 free per-galaxy offsets; and α ≥ 0.75 at 95% caps a local-density variable at 25%.
The two are separable on this very dataset and the acceleration one wins.

### What this changes, and what it does not

**Not a new refutation. Count held at 6, Bucket 0 unchanged (0).** The site lane reached the neighbouring
result from the other direction on 08-03 — evaluating the stated density law moves predicted velocities
by **2–5 orders of magnitude** and fails by *functional form* (a flat curve needs boost growing ~linearly
in r; ρ falls exponentially in a disk so `1/C` grows exponentially; no (γ, ρ_crit) reconciles them) — and
correctly declined to bump the count for the same reason. Concurred.

**What changes is the shape of the verdict, and it amends my own framing from yesterday.** 08-03 called
the sector "MOND under one substitution, `g_bar → ρ`." That puts Synchronism on ρ and MOND on g and
treats the substitution as the framework's live content.

> **The substitution is not between Synchronism and MOND. It is between Synchronism's prose and
> Synchronism's own instrument.** The empirical support claimed for `C(ρ)` was produced by `C(g/a₀)`,
> which *is* MOND. The framework's own variable sits at the far end of the region its own data exclude.

That is not a reparametrization claim. A reparametrization is a choice of coordinates; this is a number
computed under one law being quoted as evidence for another.

---

## 3. Two arithmetic corrections, both against myself

Found by recomputing rather than re-reading. Both originated in the research lane's 08-02 triage and
were inherited here on 08-03 unchecked.

**(i) The separating coefficient was wrong.** The 08-03 report and whitepaper note say the term separating
γ from the MOND point *"carries the coefficient γ(2γ−1)."* Expanding:

```
D(x,γ) = tanh(γ·ln(1+x)) − x/(x+2) = (γ − ½)·(x − x²/2) + O(x³)
```

First order in (γ − ½), not second order. At x = 10⁻⁵:

| γ | numeric | (γ − ½) | published γ(2γ−1)/2 |
|---|---|---|---|
| 0.425 | −0.075000 | −0.075000 | −0.031875 |
| 0.489 | −0.011000 | −0.011000 | −0.005379 |
| 0.500 | +0.000000 | +0.000000 | +0.000000 |
| 0.600 | +0.100000 | +0.100000 | +0.060000 |
| 2.000 | +1.499992 | +1.500000 | **+3.000000** |

~2× low near ½, 2× **high** at γ = 2. **The vanishing point γ = ½ is correct**, which is why the
conclusion survived — but a first-order departure is *easier* for data to resolve than "the O(x²) term"
implied, which cuts against the degeneracy argument it was introduced to support.

**(ii) The degeneracy partner was wrong.** The note says γ ≈ 0.489 *"is degenerate with the ρ_crit
prior."* **There is no ρ_crit in that fit.** In the deep-MOND limit `tanh(γ·ln(1 + g/a₀)) → γ·g/a₀`, so
the family depends on γ and a₀ **only through γ/a₀**:

| γ | a₀ (m s⁻²) | γ/a₀ |
|---|---|---|
| 0.2445 | 2.666×10⁻¹¹ | 9.171×10⁹ |
| 0.489 | 5.333×10⁻¹¹ | 9.170×10⁹ |
| 0.978 | 1.067×10⁻¹⁰ | 9.170×10⁹ |

μ agrees to <1% across a 4× span in γ. And the fingerprint was sitting in the registered artifact all
along: its profiled **a₀ = 5.33×10⁻¹¹ is 2.11× below its own reference McGaugh a₀ = 1.128×10⁻¹⁰** — the
compensation γ ≈ ½ requires. **The conclusion stands and is now correctly grounded.**

### The self-observation, which is the part worth keeping

My 08-03 report says, of the γ = ½ identity: *"I re-verified rather than accepted, per this lane's own
08-02 rule that a claim of exactness is a claim."* True — of the **headline**, checked to 5.55×10⁻¹⁷.
False of the **caveat in the next sentence**, transcribed unchecked and then carried to five further
surfaces by me.

**Verifying the striking claim and inheriting the boring one adjacent to it is its own failure mode**,
and it is not covered by any rule this lane currently holds. The rule that would have caught it is not
"verify exactness claims" — it is *the sentence you are least tempted to check is the one travelling on
your verification of its neighbour.*

**Propagation done properly this time** — grepped my own phrasing first, per the standing rule, and
**read the matches rather than the file list**: 20+ files matched `2γ−1` across three repos, but
`SESSION_FOCUS.md:233` is an unrelated *exponent* `(1+y)^2γ−1`, not this coefficient. Four real in-repo
instances corrected (`PREDICTIONS.md`, `SESSION_FOCUS.md:362`, the 08-02 triage, the 08-02 proposal),
plus my 08-03 report and `whitepaper_sync.json` — six surfaces. The site's **live `/coherence-function`
page** carries it too: **flagged, not edited** (sibling scope, per the 08-03 precedent).

---

## 4. Phase 1 — whitepaper

Three amendments to §5.15's status note, plus a two-surface propagation.

**The propagation is the non-obvious half, and it was two-signed.** Both `executive_summary.md` and
`conclusion.md` carry this sentence:

> "The one rung where C is measurable (**galaxy: C = g_bar/g_obs**) is the prediction target itself, and
> the **C(ρ)** prediction of it is refuted (S661, ΔBIC=+184)"

That single clause calls C measurable as an **acceleration** ratio and attributes the refutation to the
**density** law. The discrepancy is inside the sentence. `[SCOPED 2026-08-04]` markers added to both
(parity preserved, nothing deleted), stating both signs: it loosens a recorded refutation *and* voids
the sector's only quantitative support, and the stated law fails far worse when actually run. Net
effect is **unfavourable** to the framework. *"The galactic sector closed by execution" stands; the
execution that closed it is not the one cited.*

**Checked at the matches, not the file list** — the 08-03 lesson, reapplied and load-bearing again:
`0.489` occurs **0 times** in both published surfaces, so the γ figure carried no obligation at all.
`ΔBIC` did (4 hits in exec, 3 in conclusion). Sweeping on the γ figure alone would have found nothing
and closed the check.

**Back-annotation filed** — `Research/proposals/dielectric_completion_and_efe_linearity_equivalence_20260804.md`.
The site lane's 08-03 finding closes with an explicit request to register its constructive result here;
it was 24 h old and unactioned. A momentum-conserving field equation exists — `∇·[C(ρ)∇Φ] = 4πGρ` — it
answers Felten (1984) by Noether, reproduces `g = g_N/C(ρ)` under Gauss, and is **linear in Φ**. So
**EFE = 0 is the signature of that linearity, not an artifact of absent dynamics.** This *withdraws* the
08-03 premise "the galaxy sector has no field equation" while *strengthening* that retraction's
conclusion. And the same linearity forces `C(0) = 0` exactly — zero permittivity in vacuum, divergent
exterior field. **The sector's one surviving structural prediction and its worst pathology are one
property.**

| Gate | Result |
|---|---|
| Claims freeze (`--check`) | ✅ 10 claims, v1 freeze verified, exit 0 |
| Build (`make-md.sh`) | ✅ exit 0, 7,448 lines |
| Churn — **content** | **44 lines** |
| Churn — **raw** | **23,084 lines** → artifacts restored, CI builds them — **5th consecutive correct firing** |
| Lone-CR | one path, the frozen `claims/v1-snapshot/` — **12th consecutive pass** |
| `recommendations.json` | 15 insertions / 8 deletions, **raw == content** — 8th consecutive pass, no re-serialization churn |

Content verified in the built monolith before restoring, with counts **reconciled to located sources**:
`mu convention` ×4 = §5.15 + 2 markers + 1 CHANGELOG; `SCOPED 2026-08-04` ×4 = 2 markers + 2 CHANGELOGs;
`AMENDED 2026-08-04` ×3 = the three §5.15 amendments; `9.17` ×2.

---

## 5. Phase 0

| Item | Change | Why |
|---|---|---|
| **REC-2026-038** | +1 strength, +1 weakness, **held 0.93** | Sixth failure mode — and the first of a new *genus*. |
| **REC-2026-040** | +1 strength, +1 weakness, **held 0.45** | Better motivating example; blocking gate untouched. |
| **Alignment Arc** | Q1 PASS recorded, `active`, 1 of 4 stages | Explicitly **not** a publication candidate. |
| REC-033/034/035/036/037/039 | held | No new evidence. |

**REC-038 gains a mode that is qualitatively unlike its other five.** Every instance the draft carries
is continuity *failing* — a correction that does not propagate. Today's is continuity **succeeding**:
γ ≈ 0.489 travelled roughly a year across a preregistration, the whitepaper, the executive summary, the
conclusion, the site and this lane's own reports, **correctly transcribed at every hop**, having detached
from the equation it was a fit of. Nothing is erroneous at any step; what is lost is the *warrant*. That
is the more dangerous direction, and no existing instance covers it.

It also **amends my own 08-03 advice** to that draft. I recommended it stop enumerating instances and
write the asymmetry ("errors propagate; corrections must be chased") structurally. That asymmetry does
not cover provenance detachment. The stable object may be a **three-way split — content, correction, and
warrant — each with its own propagation rate.** That reframing is one day old and unstable, which is the
argument for **holding 0.93 rather than raising it**: six modes in four days, the newest a new genus, is
a subject still enumerating rather than converging.

**REC-040 held at 0.45, deliberately.** Today's finding is a *better motivating example* for the
admixture note — a program quoted a fit under the wrong variable for a year, and "how local may the
organizing variable be?" is exactly the check that forces the argument to be named. But the blocking
gate is unchanged and still weakness #1: the **external prior-art search on the α-interpolation
construction specifically remains UNRUN, six days after opening.** A stronger anecdote does not clear a
prior-art gate, and recording that explicitly is the correction to how REC-039 was opened at 0.72 on
07-29.

**Alignment arc Q1 PASS, recorded and explicitly not promoted.** Both registered kill criteria cleared
on real audits: relabeling (static — no expression takes θ_A and θ_B jointly; dynamic — outcome_A
bitwise identical 200/200 under post-hoc θ_B swap, while outcome_B follows its own setting 95/200) and
conspiracy (randomized per-run detector phases, settings drawn after emission + travel, no sensitivity
at K_SUB ×½ and ×2). S(τ) rises 0.19 → 1.29 monotonically against clean nulls (0.12, 0.11). It is not a
candidate: **S ≤ 2 throughout** (no Bell violation, as any local construction must give), the full-lock
angular form is **unresolved** (triangle RMSE 0.0796 vs cosine 0.0803 — tied at 9 points, peaks rounded
by residual own-phase), and the arc's own honesty block says Kuramoto lattices and classical triangle
constructions are known physics. Its lane scored it a **frame**-ledger event, Bucket 0 unchanged — and
that self-scoring is why it can be recorded without hedging.

**Null #1's Reading A/B fork — eighth consecutive day**, no new evidence.

---

## 6. Phase 1 — Web4: no change

One commit, `fef6e2c`, a `docs(why)` NOVA alignment-as-organogenesis note revised per CBP review — not
whitepaper-scope. It sits on the unmerged feature branch `kimi/purpose-is-relational`, which is also the
local checkout state the Archivist flagged this morning as an anomaly. No canonical-terminology exposure
(LCT/MRH/T3/V3/ATP/ADP/R6/R7 untouched). AssuranceReceipt still watch-not-act, fifth day.

---

## 7. Housekeeping

- Six surfaces corrected for the coefficient error; the **site's live `/coherence-function` page** is
  flagged, not edited — sibling scope. Carried as a new watch item so it is not lost.
- No artifacts committed; `docs/whitepaper/` and `whitepaper/build/` restored after the churn gate fired.
- `AGENTS.md` / `CLAUDE.md` still carry uncommitted GitNexus index-count drift. Supervisor scope, untouched.
- **Standing, dp's call, sixth deferral**: `.gitattributes` (`docs/whitepaper/** text eol=lf`), priced at
  ~23,000 lines per occurrence — and it has now fired five consecutive times, which is the argument.
- **Watch items**: locality operational-definition fork (10th day); dp's preprint decision (**38 days**);
  REC-036 item-4 sweep; the S654 "zero active discriminators" markers unswept on the *site* surface.

---

## Summary

The pass began with the quietest window in a week — one in-repo commit — and so turned the audit on the
sector the last three days had been making claims about. The whitepaper prints a coherence function of
**density** and, directly beneath it, quotes γ ≈ 0.489, ΔBIC +7.1/+184 and a Cassini +17.95σ that were
all computed from a compander keyed on **acceleration**, by code whose own docstring calls it "the mu
convention" and whose registered likelihood profiles a₀ rather than ρ_crit. The two laws are written
identically — `tanh(γ·ln(1 + ·))` — and differ only in an argument nobody named, which is why this
survived a year and a great deal of parameter-level auditing. It is not a new refutation and the count
holds at 6; it is a provenance defect, and it relocates the reparametrization verdict: **the substitution
is not between Synchronism and MOND, it is between Synchronism's prose and Synchronism's own instrument.**

Two numbers in yesterday's integration were also wrong — the separating coefficient (first order in
(γ−½), not γ(2γ−1); 2× off in both directions) and the degeneracy partner (a₀, not ρ_crit, since that
fit has no ρ_crit at all). Both conclusions survive, correctly re-grounded. Both were inherited from the
research lane and transcribed unchecked in the same pass where I re-verified the identity beside them to
machine precision and wrote that a claim of exactness is a claim. The rule I was following protected the
sentence I was proud of and not the one next to it.

**So what?** The program's most-audited artifact — hash-stamped, pre-registered, grid-registered in
advance, honest enough to return a refutation against its own family — is a *good* MOND
interpolating-function study. The defect is not in the instrument. It is that the equation printed above
the result is not the equation the result came from. The check that would have caught it costs one line
and this lane did not have it until today: **name the argument, not the function.**
