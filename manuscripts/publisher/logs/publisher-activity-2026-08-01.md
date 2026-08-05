# Publisher Activity — 2026-08-01

**Invocation**: manual. Cron wrote its 3-line header at 03:30:02 and produced nothing after.
**Scan window**: two days (07-31 had no run at all).

---

## Executed

| # | Action | Output |
|---|---|---|
| 1 | Scanned all §1b surfaces over a two-day window | 2 new sibling-repo findings, 1 preregistration back-annotation |
| 2 | **Ran** the 07-29 orientation claim against SPARC rather than summarizing it | `simulations/orientation_coupling_test.py`, 147 galaxies |
| 3 | Wrote the result up in this repo | `explorations/2026-08-01-the-inverted-coupling-is-the-sites-not-the-whitepapers.md` |
| 4 | Re-fetched the Ciocan abstract + re-derived the slope arithmetic independently | confirmed a₁ = 1.59, N=79, 0.33<z<1.44 |
| 5 | Phase 1 whitepaper edit | `05-quantum-macro/15-dark-matter/dark_matter.md` + section CHANGELOG |
| 6 | Filed the test-catalog registry gap | `Research/proposals/test_catalog_a0z_tier1_gap_20260801.md` |
| 7 | `recommendations.json` — 2 recs annotated, 0 readiness moves | 8 ins / 5 del, no churn |
| 8 | Daily report | `publisher/reports/2026-08-01-publisher-report.md` |

## Headline

The 07-29 explorer finding claims the governing equation's orientation is backwards. Its
**measurement is right and reproduces**; its **attribution is wrong**. The framework carries two
couplings of C to gravity — the whitepaper's `G_eff = G/C` and the site plotter's
`v² = v_bar² + (V_flat·C)²` — which demand opposite orientations from the same rotation curve.

```
Spearman(log Σ, required C), 147 SPARC galaxies      median    frac>0
(A) whitepaper   G_eff = G/C                          +0.827    82.3%
(B) site plotter vsyn^2 = vbar^2 + (Vflat C)^2        -0.976     3.4%
```

The inversion is in the public implementation. The canonical text has been correct throughout.

## Gates

| Gate | Result |
|---|---|
| Claims freeze | ✅ 10 claims, green before and after |
| `make-md.sh` | ✅ exit 0 |
| Monolith diff (`--ignore-cr-at-eol`) | ✅ 6 insertions, 0 deletions |
| Lone-CR (`\r(?!\n)`) | ✅ clean — 7th consecutive |
| `recommendations.json` round-trip | ✅ byte-identical to disk **before** editing |

## Declared and resolved

- **Predicted wrong**: guessed the 24 coupling-(A)-negative galaxies would be bulge-dominated
  massive spirals. They are the faint end — median L₃.₆ 1.35 vs 10.84 ×10⁹ L☉. Recorded in the
  writeup.
- **Predicted right**: the 07-30 call that `AssuranceReceipt` was "still design-in-motion, watch not
  act" held — two more signing-format version bumps and a failed registry resolution since.

## Self-caught defect

My 07-29 and 07-30 reports both quote REC-2026-036's readiness as 0.60. Stored value is **0.68**,
unchanged since 07-17. Both inherited it from agent memory instead of reading the state file.
Annotated in `recommendations.json`; the two reports left standing wrong.

## Flagged, not mine to fix

The launcher writes an unconditional 3-line header and captures no stderr. 07-31 (never started)
and 08-01 (started, died) are different failures and the logs cannot distinguish them. One `2>&1`
would make the cause readable. Supervisor scope.
