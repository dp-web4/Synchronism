# Publisher Daily Report - 2026-08-05

## Verdict: EFE = 0 is **not evaluable** against Chae et al. 2020. The clause I sharpened yesterday — that it "stands in tension" with a ~4σ SPARC EFE detection — is withdrawn, along with the opposite-signed claim in the research lane that the framework *fails* there. Both die on one division neither performed: the framework misses the rotation curve by **+3.1 to +4.2 dex** at the very radii where the effect is measured, and the effect is **0.046–0.083 dex**. Refutation count held at 6. "Zero remaining active discriminators" is restored, unqualified.

**Autonomous run, RUN-ID 20144** — 03:30:02 PDT / 10:30:02 UTC. Header self-identified correctly; the
liveness watch item stays retired. Today's log was header-only when I opened it, which is the expected
signature of a live run, not a failure.

---

## 1. WAKE — the in-repo window was quiet, and the movement was one repo over

| commit | repo · time (PDT) | what |
|---|---|---|
| `28251dce` | Synchronism · 08-05 02:47 | Archivist: 0 new sessions, recorded as a deliberate zero |
| `82c6c3ac` | Synchronism · 08-04 09:06 | research lane **accepts and independently verifies** my 08-04 γ↔a₀ correction |
| `591d7bd8` | Synchronism · 08-04 06:18 | back-annotation of the 08-04 dielectric-completion proposal |
| **`ca1ca95`** | **site · 08-04 08:17** | **EFE=0 is NOT refutable by Chae+2020 — today's visitor P0 reversed on execution** |

Two things worth naming before the finding.

**My 08-04 correction propagated correctly and was checked rather than accepted.** The research lane's
note (`explorations/2026-08-04-accept-gamma-a0-correction…`) records: *"I checked it myself rather than
accept on trust… Confirmed at source."* That is the loop working in the direction it is supposed to.

**And the sibling lane reversed its own P0.** A visitor persona filed a P0 claiming the EFE was "the
strongest refutation on the site"; the explorer lane executed it and concluded the opposite. A track
that kills its own headline finding on execution is the thing this program says it wants, so it was
worth spending the pass on rather than filing.

---

## 2. The finding: a residual argument evaluated against a baseline that isn't there

The framework predicts **EFE = 0 exactly** — `C = C(ρ_local)`, a uniform external field leaves local ρ
unchanged, so the Strong Equivalence Principle holds by construction. Chae et al. 2020 (ApJ 904:51)
report an External Field Effect in SPARC at ~4σ sample-level, 11σ on NGC5055, 8σ on NGC5033.

Read at that level the inference is immediate, and this repo drew it **twice, in opposite directions,
within the same 08-03 window**:

| direction | claim | surfaces |
|---|---|---|
| overclaim of *testability* | EFE=0 is "sharper and more falsifiable," and "stands in tension with" Chae's ~4σ | §5.15, executive summary, conclusion — re-sharpened by me on 08-04 |
| overclaim of *refutation* | "EFE=0 IS a genuine MOND-discriminator that the framework **FAILS**" | `PREDICTIONS.md` B2, `SESSION_FOCUS.md` |

**Chae's EFE is a residual on the outer rotation curve.** Its size is **0.046 dex (NGC5055)** to
**0.083 dex (NGC5033)** — 10–17% in velocity. The stated density law, evaluated on that same curve at
those same radii, misses by **+3.1 to +4.2 dex**.

| galaxy | Chae `e` (σ) | framework baseline error | EFE signal | baseline / signal |
|---|---|---|---|---|
| NGC5033 | 0.104 (8σ) | +3.14 to +3.71 dex | 0.083 dex | **38× – 45×** |
| NGC5055 | 0.054 (11σ) | +3.62 to +4.18 dex | 0.046 dex | **80× – 92×** |

Minimum anywhere on the grid (γ ∈ {2, 0.489} × h ∈ {0.3, 1.0} kpc): **37.8×**. MOND on the same points
and the same pipeline: **+0.002 to +0.084 dex** — *inside* the signal, which is exactly what makes
Chae's measurement meaningful for MOND and meaningless for us.

**You cannot refute "zero environmental modulation" with a 10–17% velocity deficit measured on a curve
your model already misses by 10³–10⁴×.** Correct badge: **`not-evaluable`.**

### Re-derived here, and the re-derivation is what caught the trap

Per this lane's standing rule that a claim of exactness is a claim, I rebuilt the test on the **in-repo**
SPARC mass models (`simulations/sparc_real_data/MassModels_Lelli2016c.mrt`) rather than rerunning the
sibling script.

MOND reproduced the finding's per-galaxy numbers to **≤0.014 dex** on the first attempt — which validated
the pipeline. The framework arm did not, and it was off by **exactly 1.5 dex**: the signature of a
factor of 1000. Cause: `ρ_crit = 0.029·V_flat²` is in M☉/pc³ **unscaled** (site script line 68); my
assumed 10⁻³ scaling shifts every baseline error by `½·log₁₀(1000)`.

**The error ran toward the framework.** It understated the ratio ~3×, and 20× still supports the
conclusion — so it would have survived review while being wrong. This is the third time in ten days an
unstated normalization choice has moved a headline, and the first where the drift flattered the thing
being audited.

Corrected, this lane gives **38×–92×** against the finding's **37×–53×** — agreeing at the minimum,
running higher elsewhere. **The conclusion is not sensitive to the difference**, which is the honest way
to report a reproduction that doesn't match exactly.

---

## 3. Why this is worth more than its sign

Every correction this lane has shipped for weeks has moved against the framework. This one moves
**for** it:

- **Removed**: a claimed empirical tension between EFE=0 and Chae et al. 2020 — in both directions.
- **Restored**: **"zero remaining active discriminators"** stands *unqualified*. The 08-03 text had
  qualified it ("the EFE magnitude is one, in the refuting direction"); that qualification is withdrawn.
- **Unchanged**: refutation count **6**; Bucket 0 = **0**; the constructive result that EFE=0 follows
  from the completion's linearity in Φ.

It does **not** rescue `C(ρ)`. The reason EFE=0 is unreachable is that the density law misses the
rotation curve by 3–4 dex — the *same* failure the ≤25% admixture bound and the 2–5 OOM mean-relation
result already record. **The prediction is untestable because the model is badly wrong on the underlying
observable.** That is a worse position to be in than being refuted on it; it is just not a refutation.

---

## 4. The mechanism, which is the part that transfers

The executive-summary copy of this claim carried its own warning label:

> *"Magnitude not independently re-derived here; see §5.15."*

The warning was **accurate**. It was **self-authored**. It was **never removed**. And the magnitude was
the **entire load-bearing content of the claim** — the only thing that could have decided it. The claim
reached four surfaces in two days anyway.

**Flagging is not gating.** A caveat naming the unverified quantity does not stop the claim from
travelling; it just makes the resulting artifact a correctly-hedged sentence that is wrong. If a claim's
decisive quantity is uncomputed, the options are compute it or don't propagate it.

This is the companion to the mode I recorded on 08-04 — *"the sentence you are least tempted to check is
the one travelling on your verification of its neighbour."* Together:

> **Verification does not transfer across sentences, and a hedge is not a substitute for it.**

Registered as a standing gate, filed at the source finding's explicit back-annotation request (its
action #7) — `Research/proposals/baseline_signal_gate_before_badging_a_refutation_20260805.md`:

> Before badging any published result as confirming or refuting a structural prediction, compute the
> model's own error on the raw observable, at the measurement's own radii, in the measurement's own
> units. If the baseline error exceeds the signal, the badge is `not-evaluable`.

It is one division, it is available before any literature search, and it is direction-neutral — it
killed a claimed tension and a claimed refutation in the same pass.

---

## Phase 0: Publication Recommendations

**New recommendations**: 0. **Status changes**: 0. **Updated**: 1.

**REC-2026-038** (Repository-Mediated Continuity) — +1 strength, **readiness HELD at 0.93**.

Today is the **seventh instance** of the manuscript's documented class and the first that isolates the
mechanism. The paper's thesis is that append-only narrative memory *"preserves framing, lets evidence be
bypassed, and does not propagate corrections."* Prior instances showed a warrant **silently dropped** in
transit (the 6th mode, 08-04). This one shows the warrant's absence **stated explicitly in the travelling
text** — and ignored anyway.

Timing measured, not inferred: the error reached four surfaces within hours of being written; the
reversal was executed 08-04 08:17 PDT one repo over and reached this repo on **08-05, only because a lane
went looking**. That is the asymmetry the manuscript argues for, dated and quantified.

Readiness **held rather than uplifted**, per this lane's own precedent (S662/S663): a new instance of a
documented class strengthens the evidence base without being a new publication trigger.

**Upcoming candidates**: none new. The alignment arc (Q1 PASS, 08-03) remains explicitly **not** a
publication candidate — S ≤ 2 throughout, angular form unresolved.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper — **Updated**

- **Sections reviewed**: through S691. **Changes made**: 3 sections + 3 CHANGELOGs.
  - `05-quantum-macro/15-dark-matter/dark_matter.md` — 2 clauses withdrawn, `[AMENDED 2026-08-05]` block added
  - `00-executive-summary/executive_summary.md` — `[RE-AMENDED 2026-08-05]`, restores "zero remaining active discriminators"
  - `07-conclusion/conclusion.md` — same, parity preserved
- **Also corrected (research surfaces, same claim)**: `PREDICTIONS.md` B2 (both directions), `SESSION_FOCUS.md` (2 entries)
- **Proposal filed**: `Research/proposals/baseline_signal_gate_before_badging_a_refutation_20260805.md`
- **Terminology concerns**: none. No canonical term touched.
- **Flagged, not edited (sibling scope)**: `/mond-unification` carries the same *"sharper than MOND's
  observed ~4σ EFE detection"* sentence; `/galaxy-rotation:210` and `/key-claims:545` still lead with the
  **retired** MOND-shared tie, **21 days** live.

### Web4 Whitepaper — **Current**

Repos checked: web4, hardbound-core. One commit in window (`fef6e2c`, docs(why) NOVA note, 08-03) on a
feature branch and not whitepaper-scope. **Proposals**: none. **Changes**: none.

---

## Gates

| gate | result | streak |
|---|---|---|
| claims freeze | `checked 10 claims; v1 freeze verified`, exit 0 | — |
| `make-md.sh` | exit 0, **7,476 lines** | — |
| **churn (two numbers)** | **content 66 / raw 11,594** → fired, artifacts restored, CI builds | **6th consecutive correct firing** |
| lone-CR | exactly one path, `claims/v1-snapshot/…` (frozen v1 evidence) | **13th consecutive pass** |
| `recommendations.json` | 4 insertions / 3 deletions, raw == content | **9th consecutive pass, no re-serialization churn** |

The churn gate is worth one line: content **66**, raw **11,594**. Reporting the CR-stripped number alone
would have read as "66 lines, clean." Both numbers, every time.

---

## Summary

The sibling lane executed a check that reversed its own P0, and it reversed one of mine with it: EFE = 0
is **not-evaluable** against Chae et al. 2020, not "in tension" with it, and not refuted by it. I
re-derived it on in-repo SPARC data rather than accepting it — which is how I found a normalization error
in my own reproduction that ran in the framework's favour. Six surfaces corrected, refutation count held
at 6, and **"zero remaining active discriminators" restored unqualified**. The transferable finding is
not about the EFE at all: a claim propagated across four surfaces carrying an accurate note that its
decisive quantity had never been computed, and the note stopped nothing. **Flagging is not gating.**

## So what?

This advances discovery in a way the last three passes did not: it is the first correction here in weeks
that runs *toward* the framework, and it demonstrates the ledger is not a ratchet. A program that only
ever finds itself wrong is not obviously more honest than one that only ever finds itself right — both
are single-signed. Today the same discipline that produced six refutations declined to produce a seventh,
and the reason was arithmetic available to anyone at any point in the last two days.

The open question I'd put to dp: the differential completion `∇·[C(ρ)∇Φ] = 4πGρ` (registered 08-04) is
now the **only** route by which the EFE question becomes askable at all — if it closes the 3-dex baseline
gap on SPARC, a real test exists; if it doesn't, the galaxy sector has no reachable structural prediction
left. That is a sharper and cheaper decision than the preprint question that has been open **40 days**,
and it is currently unowned.
