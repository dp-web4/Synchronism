# Phase-17 (generative axis, QC) — the compatibility design rule, transferred and sharpened (2026-06-24)

**Status:** `[ACTIVE-MRH]` — generative-axis follow-through on the QC coherence lens (dp: "take the
next steps, let's see where the path leads"). Tests whether the repo's reproducible coupling-
coherence result transfers to qubit collective coherence, in a non-circular frustrated-Kuramoto
model. **Result: the design rule transfers and *sharpens* — compatibility (not coupling strength)
gates collective coherence, and in the frustrated case there is a hard COMPATIBILITY FLOOR below
which no coupling strength achieves coherence (`K_c → ∞`), sharper than the homogeneous `p_crit ∝
1/⟨C⟩`.** Generative/applied; **Bucket 0 = 0.**
**Sim:** [`simulations/phase17_qc_compatibility_coherence_design_rule.py`](../simulations/phase17_qc_compatibility_coherence_design_rule.py) · result: `simulations/results/phase17_qc_compatibility_coherence_result.json`
**Author:** CBP-Claude (Opus 4.8).

## The mapping (and the non-trivial bit)

Qubit collective coherence ≈ **phase synchronization** (the resource QEC must protect). Oscillators
= qubits; order parameter `R = |⟨e^{iθ}⟩|` = collective coherence; heterogeneity `ω_i` = noise. The
QC-relevant non-triviality: **incompatible couplings are not merely absent — they FRUSTRATE** (a
destructive-correlation channel), so `J_ij = +1` (compatible, prob `⟨C⟩`) or `−1` (incompatible).
Frustrated Kuramoto: `dθ_i/dt = ω_i + (K/N)Σ_j J_ij sin(θ_j−θ_i)`. This is why "just couple harder"
should fail — frustration can't be out-muscled.

## Results

| test | result |
|---|---|
| **compatibility gates coherence** | sweeping `⟨C⟩` at high `K=8`: `R` rises through a threshold near **`⟨C⟩_crit ≈ 0.65`**; below it, low `R` even at high `K`. ✓ |
| **`p_crit ∝ 1/⟨C⟩`?** | `K_c·⟨C⟩` = 2.5 (`⟨C⟩=1.0`), 2.8 (`0.8`), **5.5 (`0.65`)** — *not* constant; `K_c` **diverges** toward the floor. The clean `1/⟨C⟩` is the homogeneous regime only. ✗ (refined) |
| **agent-zero** | `⟨C⟩=0.5`, `K=20` → `R = 0.078` — coupling strength alone **cannot** rescue a frustrated network. ✓ |
| sanity | `⟨C⟩=1.0` → `K_c≈2.5`, consistent with standard Kuramoto onset for `N(0,1)` frequencies. ✓ |

## The transferred design rule (falsifiable, for real device data next)

**Qubit collective coherence is gated by the *compatibility* of the coupling network, not its
strength or density.** Concretely:
- There is a **compatibility floor** (`⟨C⟩_crit`): below it, collective coherence is unreachable at
  *any* coupling strength (`K_c → ∞`, frustration dominates).
- Above the floor, the coupling needed falls with compatibility (toward the homogeneous `1/⟨C⟩`).
- **Brute-force coupling ("stronger gates / more connectivity") cannot substitute for compatibility**
  — the standard scaling axis is the wrong knob below the floor.

This says something the standard "lower physical error + bigger code" axis underweights: a
qubit-coupling/noise network with too many *incompatible* (destructive-correlation) couplings has a
hard coherence ceiling no QEC overhead removes — you must engineer the *compatibility structure*.
It is consonant with *why* decoherence-free subspaces work (DFS = encoding into the maximally-
compatible / collective subspace), and it generalizes that intuition into a quantitative floor.

## Honesty

- **Not novel physics, and not even novel dynamics:** frustrated Kuramoto having a sync floor is
  known (frustration → spin-glass-like failure to globally synchronize). **The contribution is the
  FRAMING** — *compatibility-of-coupling as the QC coherence design axis* — plus the quantitative
  floor/`1/⟨C⟩` structure, now a concrete hypothesis. **Bucket 0 = 0.**
- **Non-circular:** compatibility is an independent parameter (fraction of `+`/`−` couplings); DFS
  was *not* coded in — the floor *emerged*. The agent-zero (strength-can't-rescue) is the load-
  bearing check and it passed decisively.
- **Convergence is corroboration, not redundancy** (dp's correction): the QC field independently
  moving toward DFS / correlated-noise-awareness, while Synchronism's *independently-developed*
  coherence frame points at the same compatibility axis, is mutual support — and Synchronism adds
  the quantitative floor.

## Next step (the door between vocabulary and demonstrated novelty)

Take a **real device's coupling graph + measured correlated-error data**, define a compatibility
metric `⟨C⟩` on it, and test whether `⟨C⟩` (vs raw connectivity/error-rate) predicts logical-
coherence / QEC performance — and whether a compatibility *floor* shows up. If a real device's
logical performance tracks compatibility where the standard metrics don't, the lens has earned
*demonstrated* generative novelty (per the applied-novelty falsifier). That is reachable with
published device calibration + QEC-experiment data.
