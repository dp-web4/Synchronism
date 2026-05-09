# Proposal: TEST-04a Mechanism-Class Failure — Sign Reversal in fσ₈

**Filed**: 2026-05-09 (back-annotation from synchronism-site maintainer track)
**Source**: Visitor log 2026-05-09, Pass 4 (leading-edge researcher persona)
**Related proposal**: `session107_disfavored_by_desi_dr1.md` (filed 2026-05-05)

---

## The New Framing

The existing proposal (`session107_disfavored_by_desi_dr1.md`) documents the magnitude of
the TEST-04a failure. This proposal articulates a sharper claim: **the DESI DR1 failure is not a
parameter miss — it is a mechanism-class failure.**

This distinction matters for the research program:

- A **parameter miss** (e.g., predicted 0.418, observed 0.45) can be repaired by retuning
  σ₈(z=0) or the C_cosmic/C_galactic normalization.
- A **mechanism-class failure** (predicted suppression below ΛCDM, observed enhancement above ΛCDM)
  cannot be repaired within the same mechanism class. Any G_local/G_global suppressor predicts fσ₈
  below ΛCDM at all z. DESI DR1 observes fσ₈ *above* ΛCDM at LRG1 (z=0.51), LRG2 (z=0.71), and
  at the combined σ₈(z=0) level. The sign of the leading-order effect is reversed.

---

## The Mechanism and Its Prediction

Session 107's mechanism: G_local/G_global = C_cosmic/C_galactic, where C_cosmic/C_galactic < 1
during structure formation (lower-coherence cosmic background vs. higher-coherence galactic
environments). This **suppresses** the effective gravitational coupling → suppresses structure growth
rate f(z) → predicts fσ₈ below ΛCDM.

Numeric forecast from Session 107:
- fσ₈(z=0.51) ≈ 0.418 vs ΛCDM's 0.474 (Synchronism 12% below ΛCDM)
- σ₈(z=0) = 0.76 vs ΛCDM's ~0.81

DESI DR1 measurement (arXiv:2411.12021, Table 9):
- fσ₈(z=0.51) ≈ 0.55 ± 0.06 (observed ~17% **above** ΛCDM, not 12% below)
- σ₈(z=0) = 0.841 ± 0.034 (combined ShapeFit+BAO, 2.4σ above Sync's 0.76)

---

## Why This Is a Mechanism-Class Failure, Not a Parameter Miss

The literature context makes the classification clear. Current probes (KiDS-1000, DES Y3)
often find low-z σ₈ *below* Planck prediction — the "S₈ tension" — which is why
structure-growth suppression mechanisms have traction. Synchronism's mechanism is formally
in this class: it predicts suppression at low z, converging to ΛCDM at high z.

DESI DR1 RSD (full-shape) finds the opposite: fσ₈ is **consistent with or above ΛCDM** at
all low-z bins (LRG1, LRG2), and converges to ΛCDM at high z (ELG2, QSO are more
consistent). The redshift pattern is *inverted* from Session 107's predicted pattern.

A suppression mechanism that observes enhancement cannot be saved by retuning:
- Changing σ₈(z=0) from 0.76 to 0.84 might fix the combined constraint, but doesn't change
  the predicted sign of G_local/G_global suppression at low z.
- The redshift pattern reversal (Session 107: suppress at low z, converge at high z;
  DESI DR1: enhance at low z, converge at high z) is inconsistent with cumulative-suppression dynamics.

---

## Taxonomy of Failures (proposed addition to framework vocabulary)

The framework's current failure taxonomy (from `/honest-assessment`) does not distinguish these:

| Type | Example | Can repair by retuning? |
|------|---------|------------------------|
| **Magnitude miss** | Melting points 53% error | Yes — refine the density-to-property mapping |
| **Universality miss** | Critical exponents 2× off | Partially — need to change functional form or universality class |
| **Mechanism-class failure** | TEST-04a: sign-reversed fσ₈ | No — the mechanism class (suppressor) predicts the wrong sign |

The proposal is to add this taxonomy to:
1. The site's `/research-philosophy` validation-badge section
2. The honest-assessment page (TEST-04a section)
3. The tier-1-existing page (TEST-04a entry)

---

## The Redshift Pattern as an Additional Diagnostic

Beyond the sign failure, the redshift pattern provides a diagnostic for Session 107's
diagnosis branches (identified in the original proposal):

**Branch 1 (Sign error in G_local/G_global)**: If the mechanism *enhances* rather than
suppresses local gravitational coupling, the predicted pattern would be fσ₈ enhanced at
low z (where galactic overdensities are most concentrated) and converging to ΛCDM at high z.
This matches the DESI DR1 pattern *qualitatively*.

This is the most interesting diagnostic: it means the framework's mechanism might have the
right *structure* (local coupling coherence) but the wrong sign. Specifically: if
G_local/G_global > 1 in collapsed-structure environments (because C_galactic/C_cosmic > 1,
not < 1), the prediction would flip and become consistent with the data direction.

**Implication**: The question to investigate is whether C_galactic (coherence in galactic
halos, measured at ρ >> ρ_crit) is higher or lower than C_cosmic (coherence in the
large-scale overdensity field). Session 107 assumed C_cosmic/C_galactic < 1. DESI DR1's
direction suggests C_galactic/C_cosmic > 1 may be more accurate — i.e., dense halos are
*more* coherent than the cosmic average, not less.

---

## Recommended actions

### Research repo
1. Update Session 107's status header to distinguish: "REFUTED (sign-reversed)" from
   "DISFAVORED (magnitude)" — the former forecloses the mechanism class; the latter allows retuning.
2. Open a diagnostic session on the sign of C_cosmic vs. C_galactic in the actual DESI
   redshift range (z = 0.5–1.5). Is C_galactic(ρ_halo) > C_cosmic(ρ_LSS)? This is
   computable from the coherence function with ρ values from simulations.
3. Investigate whether a growth-*enhancement* interpretation is consistent with the rest of
   Session 107's assumptions.

### Site (handled by maintainer track)
- TEST-04a label: upgrade from "DISFAVORED at 2.4σ" to "Failed — mechanism-class: sign reversed"
- Honest-assessment: add a note explaining the distinction between magnitude miss and mechanism-class failure
- Tier-1-existing: same

---

## Why this matters for the framework's epistemic program

The S₈ tension (KiDS/DES finding low σ₈ vs Planck) is the observational context that
makes structure-growth suppression mechanisms *interesting*. If DESI DR1's high-z measurement
is correct (consistent with ΛCDM and above KiDS), then the S₈ tension may be resolving
in ΛCDM's favor — in which case suppression mechanisms as a class are losing traction.

Synchronism's mechanism is a suppressor. If suppressors are being disfavored by the data,
and Synchronism's suppressor additionally predicts the wrong sign, the framework faces two
strikes against its cosmology sector:
1. The mechanism class (suppressors) is losing ground in the literature
2. The specific mechanism (G_local/G_global suppression) predicts the wrong direction

This makes the `sign-reversal diagnostic' (Branch 1 above) the highest-priority research
question for the framework's cosmology arc — not because it might rescue the prediction,
but because understanding *why* the sign is wrong would illuminate the C_cosmic/C_galactic
relationship in a way that would generate a testable, correctly-signed replacement prediction.

Productive failure > safe summaries.

---

*Filed by maintainer track, synchronism-site, 2026-05-09.*
*WAKE phase finding: TEST-04a is a mechanism-class failure (sign reversed), not a magnitude miss.*
*This distinction is new as of the 2026-05-09 visitor log (Pass 4); the original 2026-05-05 proposal*
*documented the data correctly but did not formally establish the mechanism-class framing.*
