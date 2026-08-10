# Publisher Daily Report — 2026-08-10

**Headline**: The site lane retracted the γ~1 rationale on 08-09 and concluded, from a grep it ran and
honestly cited, that nothing in `Research/` needed correcting. The corpus carries its own γ~1
rationale in different words — and it is circular in one document and wrong-signed in the other.
**Refutation count HELD at 6. Bucket 0 = 0.** This is a demotion, not a refutation.

---

## Phase 0: Publication Recommendations

### New Recommendations
None. No new session arcs; SESSION_MAP shows the 6th consecutive deliberate zero on Synchronism core
sessions.

### Status Changes
| ID | Change | Readiness |
|---|---|---|
| REC-2026-038 | +1 strength (**12th instance, new genus**), +1 weakness | **HELD 0.93** |
| REC-2026-036 | +1 weakness (**5th ID-keyed-audit instance**) | **HELD 0.68** |

**REC-2026-038 — the 12th instance is the sharpest in the ledger, and it inverts the remedy.** Every
prior instance involved a search not run, run in the wrong repo, or run in the lane's own working root
unrecognised. This one is different in kind: **the search was run, cited, and correct, and the
conclusion drawn from it was still false.** The site lane minted the invariant *"any sentence asserting
an object does not exist should cite the grep that failed to find it"*, applied it to itself in the
same document, ran the grep, reported its output accurately — and concluded the γ~1 rationale was
"site-originated" and that "Nothing in `Research/` needs correcting for this." The corpus states the
same claim in different words, so a **phrase**-grep returning zero was read as evidence about a
**claim**. The discipline was necessary and insufficient. That bounds what procedural invariants can
buy, which is a stronger result for the manuscript than another skipped-check instance would be.
Externally directed rather than self-sourced → self-sourcing tally now **4 of 11**.

Countervailing, and now the dominant risk on this recommendation: the artifact is **untouched since
2026-07-23 (18 days)**, so none of the last eight instances are in it — and the draft's corrective arc
points toward exactly the procedural discipline this instance shows to be insufficient.

**REC-2026-036 — the catalog has no slot for a voided mechanism.** P15.4 remains a registered
prediction with an untouched "Falsified if" clause, correctly so. But the catalog cannot express
"live, falsifiable, and with no derivation behind it" — such an entry reads identically to a
well-motivated one. Mirrors the 08-09 definitional-fork gap; worse here, because the un-derived
prediction organises a 1,873-phenomenon-type sector.

### Upcoming Candidates
No change. Two-paper strategy (REC-034 ALFALFA 0.97 / REC-037 0.98) still pending dp.

---

## Phase 1: Whitepaper Review

### Synchronism Whitepaper — **Updated**

**Sections affected**: `05-quantum-macro/12-chemistry` (+ meta/CHANGELOG)

#### What was found

The 2026-08-09 site-lane proposal established that *"at γ ≈ 1 the coherence function has maximum
curvature"* is false. **Confirmed here** by independent symbolic re-derivation — `dC/dx|₍ₓ₌₀₎ = γ`
(monotone, unbounded), `max_x dC/d ln x` monotone-saturating with γ=1 at 72% of ceiling, `d²C/dx² < 0`
for all `x ≥ 0`. Used sympy rather than numerical differentiation, per the 2026-08-06 lesson. **This
check cut toward the site lane's finding, not against it.**

That proposal deferred one check — *"worth confirming against the chemistry session archive (sessions
134–2660), which this grep did not cover"* — and concluded the claim was site-originated. Running it:

- The **phrase** is genuinely absent from the chemistry archive. That much holds.
- The **claim** is present in two load-bearing corpus documents, worded differently:
  - `MASTER_PREDICTIONS.md:574` (P15.4): *"Mechanism: Correlation length diverges at γ = 1 (N_corr = 4)"*
  - `Session20_Complexity_Emergence.md:36`: *"Peaks at γ = 1.0 where N_corr = 4."*

**Defect 1 — Session #20's derivation is circular.** In `C_eff ∝ (2/γ)·(γ/2)·exp(−(γ−1)²/σ²)`, the
stated *"balance between order and disorder"* factor `(2/γ)·(γ/2)` is **identically 1**. The expression
reduces to a Gaussian centred at γ=1 by construction; the peak is the hand-written constant in
`(γ−1)²`. All five rows of the published table are reproduced by the Gaussian **alone** to a worst
deviation of **0.0060**, with σ² ≈ 0.7944 calibrated on the γ=0.5 row and the rest predicted. The
table's own symmetry is the tell: `C_eff(0.5) = C_eff(1.5) = 0.73` exactly, which a real
`(2/γ)·(γ/2)` modulation could preserve only by being identically 1.

**Defect 2 — P15.4's mechanism is wrong-signed under the framework's own formula.** `γ = 2/√N_corr`
inverts to `N_corr = 4/γ²`, so at γ=1 `N_corr = 4` — finite and small. Correlation length is monotone
in N_corr under any reading, so divergence sits at **γ → 0**, the corpus's own *"fully correlated /
rigid order"* limit. The corpus is also inconsistent with itself: P15.4 says *"diverges"*, Session #20
Part 4.1 says *"comparable to system size"*, and `N_corr = 4` supports neither.

#### Second, independent finding: the caveats propagated upward and never landed

§5.12 carried **zero** audit caveats, while the Executive Summary and Conclusion both carry extensive
S647 (method-gap) and S651 (null-model-gap) caveats **about §5.12's own headline number**. A reader who
navigates to the chemistry content saw the uncaveated 2026-02 framing — *"89% Validated"*, *"The
framework succeeds quantitatively for … 89% prediction success rate"*. This is correction-propagation
running in the unusual direction: corrections reached the summaries and never reached the material
being summarised.

#### Changes made

1. `12-chemistry/chemistry.md` — `[CAVEAT — added 2026-08-10]` block after the phase table recording
   both voided derivations, plus the S647/S651 caveat carried into this section for the first time.
2. `12-chemistry/chemistry.md` — annotated the two "succeeds quantitatively for" bullets that depend
   on them (89% rate; 1,873 types at γ~1).
3. `05-quantum-macro/meta/CHANGELOG.md` — full entry with rationale.
4. `Research/Chemistry/MASTER_PREDICTIONS.md` — `[CORRECTION]` block at P15.4; mechanism struck,
   **prediction retained**, "Falsified if" untouched.
5. `Research/Chemistry/Session20_Complexity_Emergence.md` — `[CORRECTION]` block at §1.2.
6. New: `simulations/publisher_20260810_gamma1_mechanism_audit.py` (4 checks, all reproduced above).
7. New: `Research/proposals/gamma1_corpus_rationale_is_circular_and_wrong_signed_20260810.md`.

Conservative throughout: no figures altered, no table rows removed. *"Universal coherence boundary"*
retained but explicitly relabelled as a name for the regularity rather than a derived result.

#### What this does and does not cost

**Nothing empirical.** The γ~1 clustering survives as an observation; what is gone is any derivation of
it, on either side. The chemistry sector now holds a per-sector fitted γ with N_corr never
independently measured, an unexplained clustering, a null model where a 2-parameter polynomial in Z
matches to |Δr| ≤ 0.07, a sign problem against sound velocity, and **no mechanism from either
direction**. **Count HELD at 6** — voiding a rationale removes support for a claim without
contradicting it.

#### Build verification (both numbers, per the 2026-08-01 gate)

| check | result |
|---|---|
| `make-md.sh` | **OK** — 7,591 lines, caveat present in monolith |
| content diff (`--ignore-cr-at-eol`) | **76 lines** |
| raw diff | **23,648 lines** |
| churn gate | **FIRED (11th)** — artifacts restored, CI builds |
| lone-CR gate | binary-only (PDF/PNG), no text paths |

### Web4 Whitepaper — **No change (5th consecutive structural zero)**

Scanned at `origin/main` per the §1b amendment. **4 commits** in window; **exactly 1** touched
`whitepaper/`, and it was `whitepaper/PUBLISHER_CONTEXT.md` — the Publisher's own log. The other three
are `docs/audits/` deltas (C342/C344/C346). No protocol element, no specification change, no
whitepaper scope.

**HEAD is still parked on `kimi/purpose-is-relational`** — unchanged since 08-09, so the shared
checkout remains stale and the `origin/main` discipline remains load-bearing rather than incidental.

### Terminology Concerns
None net-new. `γ = 2/√N_corr` is used here exactly as canonically defined — the finding *depends* on
the canonical definition rather than straining it.

---

## The transferable item

**A phrase-grep proves a phrase, not a claim.** The 08-09 invariant ("cite the grep that failed to
find it") is right and should stand — but it was followed exactly, self-applied, and still produced a
false provenance verdict. Before concluding "X is not in this corpus", search at least one
**paraphrase**, and prefer the claim's **load-bearing symbols** over its English phrasing: here
`N_corr = 4` finds both hits immediately, `"maximum curvature"` finds neither.

Corollary for §1b: the surfaces list needs a *how*, not only a *where*. Naming the repository is
insufficient when the query is lexical and the target is semantic.

---

## Open, not adjudicated

The 08-09 proposal flagged that `Coupling_Coherence_Experiment.md`'s acceptance test
(`p* = argmax |d²C/dp²|`, accept if `|p* − p_crit|/p_crit < 0.15`) may be degenerate, since C is
concave everywhere so `p* → 0` regardless of `p_crit`. **Still open** — it turns on whether the
*simulated* C(p) inherits the analytic concavity, which nobody has checked. Unrun; cheap; recorded so
it does not decay into folklore.

---

## Summary

The chemistry sector's γ~1 boundary now has **no surviving derivation on either side**: the site's
"maximum curvature" rationale was retracted 08-09, and the corpus's own rationale — which the site's
phrase-grep could not see — is circular in Session #20 and wrong-signed in P15.4 under the framework's
own `γ = 2/√N_corr`. Whitepaper §5.12 amended, both corpus documents back-annotated, prediction
retained. Separately, §5.12 turned out to carry none of the S647/S651 caveats that the summary sections
carry about §5.12's own headline number; that gap is now closed. **Count HELD at 6, Bucket 0 = 0.**
Web4 no-change, 5th structural zero, scanned at `origin/main` with HEAD still parked elsewhere.
