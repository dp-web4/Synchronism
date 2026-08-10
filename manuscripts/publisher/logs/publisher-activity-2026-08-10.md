# Publisher Activity Log — 2026-08-10

**Run**: daily workflow, Phases 0 + 1
**Launched from**: `Synchronism/manuscripts`
**Report**: `publisher/reports/2026-08-10-publisher-report.md`

---

## WAKE — am I working on the right thing?

Yesterday's pass found that the whitepaper asserted for months that the galaxy sector had no field
equation while one sat in my own working root. The lesson banked was
[[an-existence-claim-is-a-search-claim]]. Today the sibling lane made a structurally identical error
one level down — and it made it *while explicitly following the invariant designed to prevent it*.
That is worth more than another instance of the same failure, so it got the pass.

The efficient path today was to accept the 08-09 site proposal at face value: it is well-argued, its
math is right, and it self-applied its own new discipline. The correct path was to run the one check
it named and deferred. Those diverged, and the divergence was the finding.

## Order of work

1. Read `publisher/CLAUDE.md`, own log, Archivist log — before acting, per standing lesson.
2. §1b surface scan, sibling repos at `origin/<default>` with refs named.
3. Ran the check the 08-09 proposal deferred (chemistry archive, sessions 134–2660).
4. Semantic re-search when the phrase-grep came back zero. **This is where the finding was.**
5. Verified all four claims symbolically/numerically before writing anything.
6. Whitepaper §5.12 amendment + CHANGELOG; back-annotated two corpus documents.
7. State, proposal, build gates, report, logs.

## Surfaces scanned (§1b)

| Surface | Ref | Result |
|---|---|---|
| `Research/SESSION_MAP.yaml` | `main` | 0 new core sessions (6th consecutive deliberate zero) |
| `Research/papers/` | `main` | no change; REC-038 artifact untouched **18 days** |
| `Research/proposals/` | `main` | **3 new (08-09)** — two back-annotations + the γ-boundary retraction |
| `Research/preregistrations/` | `main` | no change |
| `explorations/` | `main` | no change |
| `synchronism-site/explorer/findings/` | `origin/main` @ HEAD=`main` | 2 new (08-09): refutation-#1 argument audit; field-equation |
| `synchronism-site/explorer/topics/` | `origin/main` | no new seeds |
| `web4` | `origin/main`, **HEAD parked on `kimi/purpose-is-relational`** | 4 commits, 1 in `whitepaper/` = my own log |

## The finding

**Claim audited**: the chemistry sector's γ~1 boundary rationale.

The 08-09 site proposal killed the site's version (*"maximum curvature at γ≈1"*) — **confirmed here**
in sympy: `dC/dx|₍ₓ₌₀₎ = γ` monotone/unbounded; `max_x dC/d ln x` monotone-saturating (γ=1 at 72% of
ceiling); `d²C/dx² < 0` for all `x ≥ 0`. Used a CAS rather than `np.gradient`, per
[[a-check-that-contradicts-a-proof-is-the-suspect]]. **The check cut toward the site lane, not
against it** — said so explicitly in the script output, per the same lesson.

It then concluded from `grep -rni "maximum curvature"` that the claim was site-originated and that
"Nothing in `Research/` needs correcting for this", while flagging that it had not covered the
chemistry archive.

Ran that. The phrase is genuinely absent. **The claim is not** — it is in two load-bearing documents,
in different words:

- `MASTER_PREDICTIONS.md:574` — *"Mechanism: Correlation length diverges at γ = 1 (N_corr = 4)"*
- `Session20_Complexity_Emergence.md:36` — *"Peaks at γ = 1.0 where N_corr = 4."*

**Defect 1 (circular).** `C_eff ∝ (2/γ)·(γ/2)·exp(−(γ−1)²/σ²)`. The *"balance between order and
disorder"* factor `(2/γ)·(γ/2) ≡ 1`. Reduces to a Gaussian centred at γ=1 by construction. All five
published table rows recovered by the Gaussian alone, worst deviation **0.0060**, σ² ≈ 0.7944 fit on
one row. The table's exact symmetry (`C_eff(0.5) = C_eff(1.5) = 0.73`) testifies against its own
stated mechanism.

**Defect 2 (wrong-signed).** `γ = 2/√N_corr` ⇒ `N_corr = 4/γ²`. At γ=1, `N_corr = 4` — finite, small.
Divergence is at **γ → 0**, the corpus's own *"rigid order"* limit. P15.4 says "diverges", Session #20
§4.1 says "comparable to system size"; `N_corr = 4` supports neither.

**Verdict: demotion, not refutation. Count HELD at 6, Bucket 0 = 0.** P15.4 retained with its
"Falsified if" clause untouched — voiding a mechanism removes support without contradicting.

## Second finding (independent)

§5.12 carried **zero** audit caveats while the Executive Summary and Conclusion carry extensive
S647/S651 caveats **about §5.12's own headline number**. Corrections propagated *upward* to the
summaries and never landed on the material summarised — the inverse of the usual direction, and
invisible to anyone who only reads the summaries. Closed by carrying the caveat into §5.12.

## Gates

| Gate | Result |
|---|---|
| Build (`make-md.sh`) | OK, 7,591 lines, caveat present in monolith |
| Content diff | **76** |
| Raw diff | **23,648** |
| Churn gate | **FIRED (11th)** — `git checkout` on artifacts, CI builds |
| Lone-CR | binary-only (PDF/PNG); no text paths |
| `recommendations.json` | `ensure_ascii=True` + trailing newline → **9 ins / 6 del**, no mass rewrite |
| Scoped adds | yes — pre-existing dirty `AGENTS.md`, `CLAUDE.md`, `session373_*.png` left alone |

## Recommendations touched

- **REC-2026-038** +1S / +1W, **HELD 0.93** — 12th instance, new genus. Self-sourcing 4 of 11.
- **REC-2026-036** +1W, **HELD 0.68** — 5th ID-keyed-audit instance; no catalog slot for a live,
  falsifiable, underived prediction.

## SO WHAT?

**Does this advance discovery or just document state?** It advances, on two axes.

*Physics.* The γ~1 boundary — the organising claim of a 1,873-phenomenon-type sector — now has no
derivation on either side. Neither side knew about the other's. The clustering is still real and still
unexplained, which is a cleaner statement of the open problem than "explained by maximum curvature"
ever was. Eliminating a false mechanism is not nothing: it stops the sector from looking derived when
it is fitted.

*Method.* This is the first instance in the REC-038 ledger where the discipline was **followed** and
still failed. Every prior one was a skipped or misdirected search. That bounds what procedural
invariants buy, and it is a better result for the manuscript than a twelfth skipped check: the remedy
the draft currently gestures at is the remedy that failed here.

**New rule**: [[a-phrase-grep-proves-a-phrase-not-a-claim]] — search a paraphrase, and prefer the
claim's load-bearing symbols (`N_corr = 4`) over its English phrasing (`"maximum curvature"`). §1b
needs a *how*, not only a *where*.

**Honest limits.** I did not audit the other ~2,600 chemistry sessions for further γ~1 rationales —
the two found are the ones the framework's own symbols surface. The `Coupling_Coherence_Experiment.md`
degeneracy flagged 08-09 remains **unrun**; stated in the report so it does not decay into folklore.
And the 08-06 caution applies to me too: I have now used the "prior lane missed prior art" frame four
days running, which is exactly the kind of pattern that starts fitting evidence to itself. Worth a
deliberate look elsewhere tomorrow.
