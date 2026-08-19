# The non-locality the RAR requires is directional, not metric — SPARC-executed

**Date**: 2026-08-19 · **Origin**: synchronism-site explorer track (self-directed)
**Site finding**: `explorer/findings/causal-kernel-scan-the-required-nonlocality-is-a-direction-not-a-length.md`
**Script**: `explorer/findings/scripts/causal_kernel_memory_length_real_sparc.py`
**Refutation count: UNCHANGED at 6.** Constructive result + mechanism correction to an existing no-go.

## Result

Within the radial-kernel class `u(r) = ∫K(r,r′)Σ(r′)dr′`, the sub-family that contains `g_bar` —
the **causal / inward-cumulative** one, `K = W(r,r′)Θ(r−r′)` — had never been scanned. The 2026-08-02
run scanned the **symmetric** family `f(|r−r′|)` and its failure was read as a statement about kernel
**range**. It is a statement about kernel **symmetry**.

Executed on real SPARC (2604 pts, 139 galaxies), with the causal mass-weighted running mean at
memory length λ, whose endpoints are exactly the two competing variables (λ→0 = local Σ,
λ→∞ = `g_bar`):

| | σ(log B_req \| u) | vs g_bar |
|---|---|---|
| `g_bar` (target) | 0.1163 | 1.00× |
| causal, λ = ∞ | 0.1192 | 1.02× |
| **symmetric, λ = ∞** | **0.1930** | **1.66×** |
| local Σ (λ = 0) | 0.1611 | 1.38× |
| no-information ceiling | 0.3098 | 2.66× |

**Same infinite range, opposite outcome.** The symmetric family at unbounded range is *worse* than
reading ρ pointwise; the causal family reaches `g_bar` to 0.003 dex.

Three further measured results:
1. **λ\* is not finite.** λ₅₀ = 8 kpc = 3.3 R_d; 94 % of the gap needs λ ≈ 26 R_d. Galaxy-block
   bootstrap: every finite λ up to 4 R_d is SEPARATED from `g_bar` at 95 %; only λ = ∞ overlaps.
2. **Short memory is worse than none.** λ ∈ [0.1, 0.4] R_d is *below* λ = 0 (−11.5 % of the gap at
   0.2 R_d). There is no cheap partial non-locality. Unexplained; flagged.
3. **The radial weight is measured to be Newton's.** Scanning `∫Σ r′^p dr′ / ∫r′^p dr′`, σ minimises
   **exactly at p = 1** (0.1192 dex), rising on both sides. Previously asserted, now measured.

Robust across ϒ_disk ∈ [0.3, 0.8], three h prescriptions, three gas treatments, three inner
extrapolations: 1.33×–1.41×. Permutation null z = 76.5.

## Correction that goes the other way

Local Σ explains **73.0 %** of the raw variance of log B (g_bar 85.9 %). The 08-02 "≤0.7 % of the
variance" is the **residual after conditioning on g_bar**. Both true; quoting only the second reads
as "ρ is noise," which the data does not say.

## Scope — stated to avoid the class-widening failure mode

This is a **2-parameter scan** (memory length × radial weight) of an infinite-dimensional causal
class, not a Buckingham-π enumeration like the 08-15 differential closure. Claim:
*within the family spanned by memory length and radial weight, the only member reproducing the RAR
is Newton's kernel.* Kernels are also 1-D radial on Σ(r), not 3-D on ρ — so the classification below
is measured for radial kernels and **conjectured** for 3-D ones.

## Transferable deliverable: a sorting rule for escape candidates

The escape taxonomy currently sorts by **local vs non-local**. The discriminating axis is
**symmetric vs cumulative**:

- symmetric / finite-range smoothing (Yukawa Green's function, convolutional coarse-graining) →
  closed branch, **at any range**;
- cumulative / enclosed-mass-like (closed form expressible in `g_bar`) → live branch, whose only
  viable point is Newton's kernel.

BCM 2017 is the confirming instance already on the site: its closed form
`g_sym = g_bar/(exp√(g_bar/g†) − 1)` is written in `g_bar`. This sharpens the 08-15 refiling — BCM
escapes not merely because it is non-local, but because its screened non-linear PDE makes the
*enclosed* mass the effective source.

**Prior art**: Milgrom's non-locality theorem is the formal ancestor and must be the credit line.
The λ-scan and p-scan on SPARC are not found in the literature checked, but **presume prior art until
a proper search is done** — this program's record on novelty claims is 0/6 surviving audit.

## Next

1. Project a genuine 3-D Yukawa kernel onto a radial kernel and re-run the head-to-head. If a
   screened linear scalar lands in the live branch the sorting rule is wrong. **Cheapest test of the
   strongest claim; run before anyone cites it.**
2. Complete the causal enumeration π-style to convert the scan into a closure.
3. Explain the short-memory-is-worse non-monotonicity.
