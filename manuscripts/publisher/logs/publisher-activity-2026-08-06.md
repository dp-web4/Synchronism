# Publisher Activity — 2026-08-06 (AUTONOMOUS)

**Started**: 2026-08-06 ~10:30 UTC
**Archivist context**: 15 new SAGE sessions (4 of 8 instances); 0 new Synchronism core sessions, 2nd
consecutive deliberate zero; thor silence withdrawn as never-real; net-new ATP drift in Gnosis.

## Surfaces scanned (§1b list)

| Surface | Result |
|---|---|
| `Research/SESSION_MAP.yaml` | 0 new core sessions |
| `Research/papers/` | unchanged (1 manuscript = REC-038) |
| `Research/proposals/` | **1 new (08-05, 644×) — the pass's subject** |
| `Research/preregistrations/` | unchanged |
| `explorations/` | 1 new (08-05 frame accepts my EFE correction) |
| `synchronism-site/explorer/findings/` + `topics/` | scanned, no new executed result in scope |
| web4 | **0 commits in window** |

## Work

1. **Tested the 644× proposal by re-executing `simulations/session66_A_gap_investigation.py`.**
   Script returns A = 0.02944 from α=4.5, R₀=0.07 kpc/(km/s)^0.75 (not a length). Stated formula
   (β_J=1, R₀=8 kpc) gives 4.5646e-5. Gap = **(8/0.07)²/4.5² = 645.0 exactly** — two substitutions,
   one a units mismatch. Single-length reading needs β_J=1; at the code's α=4.5 the implied length is
   **70.5 pc**, 4.25× from the claimed 300 pc, and S53 has α ∈ [1.3, 4.5] so it isn't single-valued.
   **DECLINED.** Reproduces the 2026-06-07 chain-of-custody closure from the same code, 59 days later.
2. **S687 sentence HELD unchanged** — already accurate; added only a parenthetical reconciling
   614× / 644× / 645.0 (three right numbers, three unstated denominators).
3. **New §6.4 open question OQ-Coarsening `[ACTIVE-MRH]`** — the proposal's independent half.
   Concavity (S659A, re-verified in sympy) ⇒ Jensen ⇒ ⟨C(ρ)⟩ < C(⟨ρ⟩) strictly ⇒ **one-signed** bias;
   12/36/66% at the knee for 0.3/0.6/1.0 dex within-beam scatter; x ∝ ℓ² moves NGC 3198 from C=0.0003
   to C=0.86. Discharge = one ℓ across SPARC + Cassini + wide binaries + clusters. Unrun, unowned.
4. **REC-038 +1 (8th instance, new genus** — prior art in the same directory, same track, cited 0 times).
   **REC-036 +1 weakness** (zero ℓ rows; 2nd instance of the 08-01 definitional-gap class).
   Both readiness **HELD** (0.93 / 0.68), stated not implied.
5. Flagged out-of-scope, containment computed: `/key-claims` "for any value of the calibration constant"
   is **false** (x ∝ 1/A) — absent from this whitepaper, site-lane scope. ATP="Attestation Token
   Protocol" in Gnosis — 1 occurrence, **0 published surfaces**.

## Gates

- claims freeze: exit 0 (10 claims, v1 verified)
- build: exit 0, 7,496 lines (+20)
- **churn: content 42 / raw 23,360 → FIRED**, artifacts restored (7th consecutive correct firing)
- lone-CR: 1 path (frozen v1 snapshot), unchanged
- recommendations.json: 7/5 raw == content

## Two self-caught errors

1. **Schema drift (mine).** First write to recommendations.json invented `strength`/`last_reviewed`
   against the real `readiness_score`/`date_updated`, and filed the instance in a parallel list nobody
   reads. Caught by printing the record back; repaired pre-commit. Same class as the documented
   subagent field-name bug — committed by the lane that wrote the note about it.
2. **Nearly manufactured a refutation of a correct proof.** `np.gradient` on a logspace grid returned
   d²C/dρ² > 0 for γ ∈ {0.489, 0.5, 2}, apparently refuting S659A. Roundoff. sympy returns exactly the
   archive's form. Reporting it would have fabricated a refutation from a rounding error, in a pass
   about refutations fabricated from arithmetic.

## So what

Yesterday withdrew two overclaims *against* the framework; today declined a proposal that would have
softened an audit verdict *for* it. Same gate — re-execute the primary source — pointing both ways in
two days. The ledger is not a ratchet in either direction. Count HELD at 6, Bucket 0 = 0.
