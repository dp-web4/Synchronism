# Repository-mediated autonomous science

This package accompanies the paper, “Repository-Mediated Continuity in
Memoryless Autonomous Research: A Longitudinal Case Study of Claim Formation,
Correction, and Epistemic Drift.”

The paper studies a fixed public snapshot of `gnosis-research` and related
cross-repository corrections in Synchronism. It treats commit-subject counts
as a descriptive corpus profile, not as measures of scientific validity.
Named claim lineages were manually audited and are explicitly separated from
the lexical profile.

## Reproduce

At the source repository snapshot:

```bash
git log --reverse --format='%H%x09%aI%x09%s' > gnosis_commit_log.tsv
python3 analysis/corpus_profile.py \
  --log-file gnosis_commit_log.tsv \
  --repository-head 7c1c16d03ef0b1ac2b1c36468244b9bc54366913 \
  --tracked-file-count 600
```

Alternatively, run directly against a checkout:

```bash
python3 analysis/corpus_profile.py --repo /path/to/gnosis-research
```

The first form is preferred for archival reproduction because the TSV can be
hashed and reviewed independently. `reproducibility-manifest.json` records the
snapshot, expected aggregate output, coding limitations, and files needed for
manual lineage verification.

## Contents

- `paper.md` — manuscript.
- `analysis/corpus_profile.py` — standard-library lexical profiler.
- `analysis/claim_lineages.csv` — manually audited claim/adjudication pairs.
- `reproducibility-manifest.json` — snapshot and expected-results contract.

The package does not archive the companion repository itself. Commit hashes
are the stable locators; anyone reproducing the audit must obtain the public
repository and verify the head before interpreting the results.
