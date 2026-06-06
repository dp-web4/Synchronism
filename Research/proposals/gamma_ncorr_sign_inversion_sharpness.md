# γ = 2/√N_corr Sign Inversion: Correlation Does Not Sharpen the Transition

**Filed**: 2026-06-06  
**Track**: Maintainer WAKE → Research proposal  
**Status**: Open structural diagnosis

---

## The Finding

γ = 2/√N_corr has two acknowledged problems:
1. The 1/√N ansatz invokes a CLT-style scaling that is self-contradictory because N_corr counts *correlated* DOF (CLT requires *independent* variables). Already flagged on the site.
2. **The factor of 2 is not derived from first principles.** Already flagged on the site.

**New structural issue (2026-06-06 visitor Pass 3):** The *sign* of the N_corr→sharpness mapping is inverted, independent of any prefactor or scaling law.

---

## The Inversion Argument

In statistical mechanics, fluctuation theory tells us that the width of a transition scales as 1/√N — where N is the number of *independent* degrees of freedom. More correlation (larger N_corr) → fewer effective independent units → **narrower width → sharper transition**.

But in the formula γ = 2/√N_corr, more correlation means larger N_corr → smaller γ → **flatter tanh**. The tanh slope at argument zero is γ; larger γ gives a sharper transition. So:

- Fluctuation-theory logic: more correlation → sharper transition
- γ formula outcome: more correlation → smaller γ → flatter transition
- **Conclusion: the sign is inverted.**

The presets make this visible:

| System | N_corr | γ | Real sharp transition? |
|--------|--------|---|------------------------|
| Ideal gas | 1 | 2.0 (sharpest) | No — gases have no phase transition |
| BCS superconductor | 10⁷ | 6.3×10⁻⁴ (flattest) | **Yes** — a real Tc with sharp second-order transition |

The BCS superconductor — the canonical example of a *sharp* phase transition at a real Tc — gets the *most gradual* tanh slope in the model. The ideal gas — which has no phase transition at all — gets the *sharpest* γ.

The dimensional origin of this inversion: 1/√N is a *fluctuation width* (smaller = sharper). It was placed in a *rate* slot inside tanh (larger = sharper). The slot inverts the sign.

---

## Why This Matters

1. The site already disowns the CLT derivation ("CLT governs iid, not correlated DOF"). That disavowal leaves the 1/√N ansatz unmotivated dimensionally. This finding shows that even if one found a motivation for 1/√N scaling, it would need to appear *in the denominator* to produce the right sign (more correlation → sharper), not in γ = 2/√N_corr where larger N_corr makes the transition *flatter*.

2. The Phase Boundary Visualizer hard-codes this inversion visually — BEC/BCS in the "flattest-γ" bin despite being the systems with real Tc.

3. This is not a calibration miss — it is a structural constraint: no rescaling of the prefactor (no "factor of 3 instead of 2") can fix a sign error.

---

## Resolution Options

**A. Invert the formula:** Use γ = 2√N_corr (not 1/√N_corr). This would map ideal gas (N_corr=1) to γ=2 and BCS (N_corr=10⁷) to γ≈6,320. But γ≫1 gives tanh(γ·u) that saturates immediately — C≈0 everywhere except right at ρ=ρcrit. This may fail all galaxy applications where gradual transitions are needed.

**B. Reinterpret N_corr:** If N_corr measures *uncorrelated* units (not correlated ones), then large N_corr = many independent parts = low collective behavior = gentle transition, and the sign is consistent. But this inverts the stated physical meaning.

**C. Accept as scope restriction:** γ=2 (N_corr=1) applies specifically to the galaxy regime, where there is no collective ordering. The formula is correct *for that regime only*, and BCS/BEC presets are out-of-scope. The γ-calculator should say so and remove the BCS/BEC presets as in-domain examples.

**D. Find a different N_corr→γ relationship entirely:** The CLT route is self-contradictory; the 1/√N route has inverted sign. A correct derivation from first principles may look different.

---

## Proposed Site Action

Add a note to the γ-calculator and γ=2/√N_corr descriptions:

> **Sign caveat:** 1/√N is a fluctuation *width* (more correlation → smaller width → sharper transition). But γ = 2/√N_corr makes more correlation → smaller γ → *flatter* tanh. The sign of the analogy is inverted. The framework assigns the sharpest transition to the least-correlated system (ideal gas, γ=2) and the flattest to the most-correlated (BCS superconductor with a real Tc, γ≈6×10⁻⁴). This is a structural inversion, not a calibration miss. It remains an open problem to motivate γ=2/√N_corr from a derivation that gives the correct sign.

---

## Relation to Existing Proposals

- `gamma_definitional_collision_regime_label_inversion.md` — identifies two γ roles (universal constant vs operational parameter) without diagnosing the sign inversion
- `c_rho_no_inflection_for_positive_density.md` — identifies no inflection for ρ>0; the sign inversion is the N_corr→γ analogue of this structural issue
