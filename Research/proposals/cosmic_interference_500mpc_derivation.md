# Proposal: Derive the ~500 Mpc Cosmic Interference Scale

**Filed by**: Site maintainer, 2026-04-25  
**Trigger**: Four-persona visitor review (Pass 4 researcher) flagged this as the one test they would actually run — and the one that currently lacks a derivation. TEST-07 claims ~500 Mpc cluster-separation oscillations but links to a 404 page and provides no derivation of the scale.

---

## The Gap

TEST-07 in the Tier-1 test catalog states:

> *Galaxy cluster separations show oscillatory modulation at λ ~ 500 Mpc. Kill: No oscillations above 3σ out to 2000 Mpc.*

This is listed as the most distinctive novel prediction — not shared with MOND, not a reparametrization, and in principle detectable in existing survey data (SDSS, DES, DESI). A leading-edge researcher visiting the site specifically said:

> *"If derivable from a coherence wavelength formula, this is the test that matters most for my field."*

And yet:
- The dedicated page `/cosmic-interference` returns 404
- No derivation of the 500 Mpc scale appears anywhere on the site
- The scale is not attributed to any session number
- BAO is at ~150 Mpc (sound-horizon scale at recombination); 500 Mpc is well above that and in no known cosmological harmonic

Without a derivation, "500 Mpc" is a number someone wrote down, not a prediction.

---

## The Research Question

**What in the Synchronism framework picks 500 Mpc as a coherence wavelength?**

More specifically:
1. Does the coherence framework have a natural large-scale oscillation length derived from ρ_crit, γ, H₀, or c?
2. The MRH is the "coherence horizon" for a system — what is the analogue for large-scale structure? A cosmological MRH would set the scale at which coherence correlations between galaxy clusters oscillate.
3. Could the 500 Mpc scale emerge from the Jeans-criterion derivation of ρ_crit applied at the cluster scale rather than the galaxy scale?
4. Is there a standing-wave argument: if the coherence field has a dispersion relation ω(k), what k corresponds to λ ~ 500 Mpc, and is it derivable from the framework's parameters?

---

## Why This Matters

If the 500 Mpc scale is derivable, it becomes:
- A **symmetric falsifier** (the researcher's term) — positive detection at the predicted scale confirms, null kills
- Currently it's an **asymmetric falsifier** — any oscillation anywhere in the cluster regime would confirm, only a global null kills
- The asymmetric design is the single most common complaint from the expert visitor across multiple sessions

If it's not derivable from the framework, the test should be removed or labeled "exploratory" rather than "Tier-1 existing data."

---

## Suggested Approach

1. **Check existing sessions**: Search the archive for "500 Mpc" or "coherence wavelength" or "large-scale structure oscillation" — is there a session where this scale was derived or motivated?

2. **Attempt the derivation from ρ_crit**: ρ_crit = A·V²_flat sets the galaxy-scale coherence threshold. For clusters, V_flat → σ_cluster (velocity dispersion). Does applying the same Jeans-criterion formula to cluster-scale dynamics give a characteristic wavelength of ~500 Mpc?

3. **Check the MRH at cluster scale**: The MRH horizon for a galaxy cluster is set by the coherence parameter γ and the cluster's dynamical time. Does the MRH of a typical cluster (10^14 M☉, σ ~ 1000 km/s) correspond to ~500 Mpc when translated to a correlation length?

4. **Failure mode to document**: If the scale cannot be derived and is instead a free parameter, that is a productive negative result — it tells us the test is "exploratory" and should be labeled accordingly. A labeled exploratory prediction is more honest than an unlabeled one.

---

## Expected Output

- A session in the archive: `SessionNNN_500Mpc_Coherence_Wavelength.md`
- Either: a derivation that anchors 500 Mpc from framework parameters (with ±uncertainty), or a documented failure showing the scale is a free parameter
- A site update: either a `/cosmic-interference` page with the derivation, or a badge downgrade to "Speculative" on TEST-07

---

## Connection to Other Open Issues

- The BAO ~10⁻⁴ shift (TEST-04) has the same derivation gap: the order of magnitude is stated but not derived
- The wide-binary density coefficient (TEST-02) lacks a predicted magnitude
- All three are in the same asymmetric-test design class — fixing one develops the pattern for the others
