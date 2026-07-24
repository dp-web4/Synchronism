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
and rendering require only the Python standard library. The v1 preservation
hashes live in `v1-freeze-manifest.json`.
