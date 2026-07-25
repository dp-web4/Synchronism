# Claim lineage

This directory is the machine-readable source for Synchronism's **v2 public
claim surfaces**. It is deliberately separate from the historical whitepaper:
the current whitepaper remains a v1 record of how the project reasoned,
including its reversals and propagation failures.

`claim_ledger.json` begins the migration from the hand-maintained
`PREDICTIONS.md`. It is authoritative for the claims it contains, but it does
not claim full coverage of that larger ledger yet. New v2 public prose must
not promote a claim absent from this file.

Each record separates:

- the umbrella ontology from a concrete realization;
- empirical status from novelty status;
- the current statement from its lineage and evidence;
- positive findings from nulls, refutations, and open obligations; and
- historical source prose from generated public prose.

Generate the public surfaces:

```bash
python3 claims/render_claim_surfaces.py
```

Check that committed surfaces are current without modifying them:

```bash
python3 claims/render_claim_surfaces.py --check
```

The renderer validates the ledger before writing. JSON is used so validation
and rendering require only the Python standard library plus Git. The v1
preservation hashes live in `v1-freeze-manifest.json`.

Freeze hashes are computed from blobs at the manifest's `frozen_at_commit`,
not from checkout bytes. This makes verification independent of line-ending
conversion and excludes untracked build products. A separate `git diff`
check ensures that the current tree has not changed any frozen path.

Frozen paths must therefore be paths nothing still writes to. `PREDICTIONS.md`
is the exception that proves it: it is the live prediction register, so its v1
state is preserved as a snapshot copy under `v1-snapshot/` and the live file
stays free to move. Freezing the live path instead made the gate fail on the
next research commit, and because `verify_v1_freeze()` runs before rendering,
that failure aborted the renderer rather than just reporting drift. Freeze
snapshots, not working files.
