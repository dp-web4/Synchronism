# Phase-16 (velocity-anisotropic phase-transition threshold vs isotropy bounds) (2026-06-24)

**Status:** `[REFUTED]` — adjudicates dp's phase-transition framing (Phase-15 unifying frame)
turned into a *directional* prediction: if pattern-identity thresholds are **substrate-fixed**
(absolute / CMB-frame), then translation at v ≈ 370 km/s relative to the CMB dilates a pattern's
internal frequencies against those fixed thresholds **anisotropically**, shifting phase-transition
points with direction on the sky. **Result: the predicted anisotropy is a fractional-frequency
anisotropy of ≈ (1/2)(v/c)² ≈ 7.6×10⁻⁷ (isotropic dilation) with a quadrupole anisotropic part
≈ (v/c)² ≈ 1.5×10⁻⁶ and a dipole ≈ v/c ≈ 1.2×10⁻³. This maps directly onto the
optical-cavity / Hughes–Drever / clock-comparison isotropy bounds (≈10⁻¹⁷…10⁻²²), where it sits
11–19 orders of magnitude ABOVE the limits. REFUTED in the absolute-threshold case; collapses to
standard SR (no novel content) in the co-moving case. Either fork eliminates it.** This
**consolidates the Phase-13 absolute-time ⊥ Lorentz pattern** (a directly-testable, already-excluded
instance of the same preferred-frame tension). Bucket 0 unchanged (0).
**Sim:** [`simulations/phase16_velocity_anisotropic_threshold.py`](../simulations/phase16_velocity_anisotropic_threshold.py) ·
result: `simulations/results/phase16_velocity_anisotropic_threshold_result.json`
**Author:** CBP-Claude (Opus 4.8), adjudicating dp's Phase-15 phase-transition framing.

> This is a *second* Phase-16 thread, distinct from
> [`phase16-door3-dissipation-vs-oscillation`](2026-06-24-phase16-door3-dissipation-vs-oscillation.md):
> that one treads **door #3** (dissipation/secular). This one treads the **door #2 / preferred-frame**
> consequence of the Phase-15 identity-as-phase-transition picture and adjudicates it against the
> existing Lorentz-invariance bounds the B7/Phase-12/13 notes already catalogued.

## 1. The mechanism and the magnitude

The Phase-15 unifying frame (dp 2026-06-24): a pattern's internal frequencies are observed as
**temperature**; the pattern "holds" only while those frequencies sit within a **quantized band**
(its "fundamental frequency at which it can exist"). The band edges — the **phase-transition
thresholds** — are taken to be **substrate properties (global-/CMB-frame-fixed)**. Translation at
velocity `v` relative to the substrate dilates *all* of the pattern's internal frequencies by the
shared factor (Phase-5 light-clock; Part-A shared-`S` result):

```
S = sqrt(1 - beta^2) = 1/gamma ,    beta = v/c
```

Measured against the **fixed** thresholds, the pattern's frequencies drift by

```
delta f / f = 1 - S = 1 - sqrt(1 - beta^2)  ≈  (1/2) beta^2 .
```

For `v = 370 km/s` (solar-system barycenter vs the CMB rest frame, the Planck dipole):

| quantity | value |
|---|---|
| `beta = v/c` | `1.234×10⁻³` |
| `gamma − 1` | `7.616×10⁻⁷` |
| isotropic dilation `1 − S = (1/2)beta²` | **`7.62×10⁻⁷`** |

**Angular / dipole structure.** Because `v` is a *directed* velocity, the internal-frequency
deviation viewed across the boost carries the standard relativistic-Doppler angular dependence
`f_obs/f₀ = 1/[gamma(1 − beta cosθ)]`, with θ between `v` and the line of sight. Expanding to O(beta²):

| multipole | amplitude | value | character |
|---|---|---|---|
| MONOPOLE | `(1/2)beta²` | `7.6×10⁻⁷` | isotropic time dilation (velocity-magnitude only) |
| **DIPOLE** | `beta` | **`1.2×10⁻³`** | odd, aligned with `v` (the CMB-dipole direction) |
| **QUADRUPOLE** | `beta²` | **`1.5×10⁻⁶`** | leading *anisotropic* (even) modulation |

So the **anisotropic** content — a phase-transition threshold that depends on the lab's *orientation*
relative to the CMB dipole — appears at order `beta²` (quadrupole, ~1.5×10⁻⁶), with an odd dipole
piece at `beta` (~1.2×10⁻³) that an annual-modulation / sidereal analysis isolates. **This is exactly
the multipole structure (sidereal dipole + quadrupole, keyed to the CMB frame) that modern isotropy
experiments are built to detect.**

## 2. Which observable, which bound

The predicted signal is: *a frequency-defined identity threshold whose value depends on velocity and
orientation relative to a preferred (CMB) frame* — i.e. a **fractional-frequency anisotropy tied to
the CMB dipole.** That is the precise observable that the modern Michelson–Morley / Lorentz-invariance
program bounds:

| experiment class | what it bounds | bound on fractional-frequency anisotropy |
|---|---|---|
| Rotating optical cavities (modern MM; Nagel et al. 2015, Herrmann et al. 2009) | orientation-dependent shift of an EM resonance vs a reference | **~10⁻¹⁸** (dim-4 SME `c_μν`) |
| Hughes–Drever / clock-comparison (nuclear & atomic) | orientation/velocity dependence of energy-level splittings | **~10⁻²²** |
| Atomic-clock-comparison networks | sidereal modulation of transition frequencies | **~10⁻¹⁸…10⁻²¹** |

These are precisely the bounds the **B7 / Phase-12 note** already catalogued for **dim-4 LIV**
(species-dependent limiting speed / SME `c_μν`): "bounded at ~10⁻¹⁸–10⁻²² (cavity / Hughes–Drever /
clocks), **not** Planck-suppressed and **reachable**."

**The comparison is brutal.** The predicted anisotropy (quadrupole ~1.5×10⁻⁶; dipole ~1.2×10⁻³) sits

- **~11 orders** above the loosest relevant bound (~10⁻¹⁷),
- **~16 orders** above the tightest (~10⁻²²) for the quadrupole, and
- **~19 orders** above the tightest for the dipole.

A `(v/c)²` ≈ 10⁻⁶ frequency anisotropy is not marginally excluded — it is excluded by **twelve-plus
orders of magnitude**. If real, every atomic clock on Earth would show a `O(10⁻⁶)` sidereal/annual
wobble keyed to the CMB dipole. They show nothing to `10⁻¹⁸`.

### Is a *collective* phase-transition threshold a distinct observable that escapes these bounds?

This was the one honest escape to check (the prompt's "verify whether a collective phase-transition
threshold is bounded the same way"). **It does not escape — for three independent reasons:**

1. **The experiments already measure many-body / collective frequency standards.** An atomic-clock
   transition is a collective electronic state; an optical cavity resonance is a collective EM mode of
   a macroscopic crystal whose *physical length* is set by collective interatomic phase relations; a
   Hughes–Drever signal is a nuclear *collective* level splitting. If a substrate-fixed threshold
   dilated internal frequencies anisotropically, these *are* the systems it would shift. Collectivity
   is not a loophole; it is the regime already bounded.

2. **The threshold is defined *through* internal frequency vs a fixed scale** — which is the literal
   definition of every one of these experiments (a frequency/length referenced to a standard). A
   "phase-transition threshold" expressed in temperature/frequency is not ontologically distinct from
   "a transition frequency"; the experiments bound the velocity/orientation dependence of *any*
   frequency-referenced scale, named or not.

3. **A phase transition is a *more* sensitive amplifier, not a less-bounded observable.** Near a
   critical point, response diverges; a 10⁻⁶ shift in the control parameter (the dilated internal
   frequency) near a sharp band edge produces an *order-one* change in whether the pattern holds. So
   if anything the collective/critical framing makes the signal *easier* to see, not exempt. (Real
   critical systems — superconducting `T_c`, lambda points, lasing thresholds — show no such CMB-keyed
   anisotropy to the precision they are measured.)

There is no version of "it's a collective threshold, so the LIV bounds don't apply." The bounds apply.

## 3. The escape fork

The anisotropy exists **only if the thresholds are absolute (substrate-/CMB-fixed).** The fork:

- **(A) Thresholds CO-MOVE with matter (transform with the pattern's rest frame).** Then there is no
  preferred frame, no sidereal/CMB-keyed modulation; the `(1/2)beta²` dilation is just the ordinary
  SR time dilation every co-moving observer cannot detect on its own clock (Part-A: internal *ratios*
  are invariant, max deviation 2.2×10⁻¹⁶). **Outcome: standard SR, zero novel content.** The
  prediction dissolves into a true statement that is not Synchronism-distinctive.

- **(B) Thresholds are ABSOLUTE (preferred-frame).** Then the anisotropy is real at `O(beta²)` and
  **directly testable** — and **already excluded** by the cavity/clock/Hughes–Drever bounds by
  11–19 orders. **Outcome: refuted.**

This is the **same fork as Phase-13**, made concrete and *cheaply testable on a benchtop*: the
absolute-time / single-global-clock commitment that gives the framework its distinctive content is
precisely the commitment that makes thresholds frame-fixed (fork B) → which is precisely what
Hughes–Drever/cavity/clock isotropy excludes. Keeping the thresholds frame-independent to evade the
bound (fork A) is keeping Lorentz boost invariance — i.e. *giving up* the absolute frame that was the
distinctive content. **There is no fork that is both novel and surviving.**

## 4. Verdict — (a) REFUTED, consolidates Phase-13

**Verdict: (a) REFUTED.** A velocity-/orientation-anisotropic phase-transition threshold at the
`(v/c)²` ≈ 10⁻⁶ level (with a `(v/c)` ≈ 10⁻³ dipole) is excluded by existing isotropy / Lorentz-
invariance experiments (optical-cavity MM ~10⁻¹⁸, Hughes–Drever/clock-comparison ~10⁻²²) by
**11–19 orders of magnitude.** The "collective threshold is a distinct observable" escape fails: the
bounding experiments already use collective frequency standards, the threshold is frequency-referenced
by construction, and criticality amplifies rather than hides the signal. The only fork that survives
the bound (co-moving thresholds) **collapses to standard SR with no novel content** — so there is no
surviving novel reading. **This is Bucket-2 (a tested-and-failed elimination), and it CONSOLIDATES the
Phase-13 pattern** (`absolute time ⊥ a sub-10⁻²² Lorentz-invariant world`): Phase-13 was a loop-level,
interactions-unspecified *obligation*; this is the same tension cashed out as an **already-measured,
tree-level, benchtop refutation** of the absolute-threshold reading. The prior (dp's, "likely (a)")
holds, and the verification asked for is done: the collective threshold is **not** less-bounded.

**The productive elimination:** this pins down *which* commitment dies. It is not the Phase-15 Part-A
clock-universality result (that is fork A — pure SR, fine, and a genuine derivation). It is the leap
from "all internal frequencies dilate by a shared `S`" (true, undetectable internally) to "...against
*absolute* thresholds, so identity transitions are frame-anisotropic" (false, excluded). The shared-`S`
universality that *derives* clock-universality (Phase-15 Part A) is the very thing that **prevents** an
internal observer from seeing the dilation — and the moment you make the thresholds external/absolute
to recover a directional effect, you are in the Hughes–Drever-excluded preferred-frame regime. The two
results (Phase-15 Part A = win; this = refutation) are the two sides of the same shared-`S` coin.

## Honesty / net

- **Bucket 0 = 0**, unchanged. This adds **one Bucket-2 row** (refuted preferred-frame anisotropy,
  consolidating Phase-13). No bet added to Bucket 1 (the surviving fork is SR, not novel).
- **Anti-overclaim:** the magnitude is elementary SR kinematics ((1/2)β²); the *only* Synchronism-
  specific step is the absolute-threshold assumption, and that step is what the data excludes. No
  attempt to rescue it via "collective threshold ≠ frequency standard" — that escape was checked and
  fails on three counts (§2).
- **Productive failure:** a benchtop-already-done refutation of a clean, directional, Synchronism-
  distinctive consequence is *more* valuable than an untestable one — it eliminates the absolute-
  threshold reading of identity-as-phase-transition decisively, and localizes the survivable content
  to the (non-novel) co-moving / SR reading. Door #2's preferred-frame face is closed at tree level
  here, exactly as Phase-13 predicted it would be at loop level.
