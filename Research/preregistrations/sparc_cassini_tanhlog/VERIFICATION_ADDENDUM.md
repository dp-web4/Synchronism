# Heterogeneous verification addendum

**Verifier:** Claude (Anthropic), acting as pre-merge heterogeneous reviewer<br>
**Date:** 2026-07-23<br>
**Verified tree:** PR #2 head `649d4528` (fresh git worktree; no working-tree state reused)<br>
**Runtime:** Python 3.12, numpy 1.26.4, scipy 1.11.4 (matching the recorded execution runtime)

This addendum records an independent re-execution of both registered pipeline
stages from the committed tree alone, performed before merge. It adds to the
frozen record; it modifies nothing Codex committed.

## Reproduction

Both commands from `RESULT.md` were re-run unmodified:

```bash
python3 simulations/sparc_tanhlog_profile.py --output <fresh profile>
python3 simulations/sparc_cassini_joint.py \
  Research/preregistrations/sparc_cassini_tanhlog/sparc_profile.json \
  --output <fresh joint result>
```

Results:

- **SPARC profile:** every numerical field of the regenerated profile is
  identical to the committed `sparc_profile.json` — all 79 grid rows, the
  2,807-row selection, `gamma_best = 0.489` with profiled
  `a0 = 5.33265e-11 m s^-2`, fixed `gamma = 2` Delta BIC `+184.0445`,
  free-family Delta BIC `+7.1069`, and the profile-convention
  `Delta BIC <= 10` interval `gamma = 0.425` through `0.600`.
  The only differing field is the data hash discussed below.
- **Joint execution:** the regenerated joint result is identical to the
  committed `joint_result.json` in every field — all rows, quadrature
  convergence histories, benchmark tables, sensitivity summaries, outcome
  checks, and the branch-A verdict — except `execution_parent_commit`
  (this rerun was executed from the PR merge head rather than `e05e3582`).
  Exit status 0 (robust empty intersection) reproduced.
- Registration-fidelity checks passed: the executed grid satisfies the
  registered spacing and inclusion rules; the Cassini interval matches the
  registered `[-1.928, 5.128]e-27 s^-2`; all eight mandatory sensitivities
  are present; both instrument amendments changed numerical machinery only,
  with branch D honored on both pre-amendment failures; the amendment-2
  short-circuit is one-directional (it can only record SPARC exclusions,
  never create or remove an intersection member).

## Data-hash provenance note

The registered SPARC source-data hashes in `PREREGISTRATION.md` section 6 and
in the committed `sparc_profile.json`
(`MassModels_Lelli2016c.mrt = 9108994b…`, `SPARC_Lelli2016c.mrt = 5aa0501f…`)
are hashes of the execution machine's working-tree bytes, which carry CRLF
line endings on that checkout. The committed git blob for
`simulations/sparc_real_data/MassModels_Lelli2016c.mrt` hashes to:

```text
fd56b08f93ee6be38bb5b2305a16a2dc57d98ee29034223c5d14ecf8f50759d0
```

A fresh clone with LF checkout therefore regenerates a profile whose
`data.sha256` differs from the committed value. **This is a line-ending
byte-variant of the same data, not a provenance failure**, and it is
scientifically inert here: the loader tokenizes on whitespace, so both
variants parse to the identical 2,807-row selection — proven by the
bit-identical reproduction above, which consumed the LF variant. Branch E
(SPARC provenance failure) is **not** triggered. Future registrations should
hash committed git blobs (see `claims/v1-freeze-manifest.json` practice after
its 2026-07-23 correction) so that recorded hashes are checkout-independent.

## Ledger propagation

The claim ledger gains `physics.sparc-cassini-tanhlog-qumond` (refuted
realization; umbrella untested) in the same change set as this addendum, and
the `methodology.archival-confirmation` registration obligation is recorded
as discharged. Public surfaces are regenerated from the ledger by
`claims/render_claim_surfaces.py`.
