# Proposal: The Mean-Field Diagnosis — C(ρ)'s Failures Are One Failure, Not Three

**Filed by**: Site maintainer, 2026-04-27  
**Trigger**: Pass 3 graduate student visitor (2026-04-27) identified the unifying explanation for C(ρ)'s documented failures. This is a research gap, not a site gap.

---

## The Observation

The Synchronism site's Honest Assessment documents three independent-looking failures:

1. **Critical exponents ~2× off** (β_predicted ≈ 0.5 vs β_observed ≈ 0.326 for 3D-Ising class)
2. **Melting points 53% off** (mean absolute error on solid-liquid transitions)
3. **T_c 6.5× off** (superconductor critical temperatures)

A condensed-matter grad student reading Pass 3 notes:

> *"These are not three failures. They are one failure. The failures are exactly what any CMT theorist would predict for a `tanh` ansatz applied without fluctuation corrections. C(ρ) = tanh(...) is a mean-field order parameter — the Ising mean-field self-consistency solution. Mean-field has known failures: β_MF = 0.5 (vs ~0.326 in 3D Ising), ν_MF = 0.5 (vs ~0.630 in 3D Ising), and gross errors on first-order transitions (melting) where mean-field is famously unreliable. A factor-of-2 mismatch on critical exponents is the diagnostic signature of an uncorrected mean-field theory."*

This is the correct diagnosis. The site treats these as separate bullets in a failure list. They should instead be presented as a single structurally-identified failure class: **C(ρ) is a mean-field order parameter and inherits mean-field's failure modes below its upper critical dimension.**

---

## Why This Matters

Two reasons:

### 1. Diagnostic clarity

Once identified as mean-field failures, the site can say something precise: *"C(ρ)'s critical-exponent and melting-point failures are the expected failures of a mean-field theory without RG corrections. This is a known class of failure, not mysterious scatter."*

This is MORE credible than listing failures as separate bullets — it shows the site understands WHY, not just WHAT.

### 2. It opens a concrete research program

Mean-field failures are fixable in principle via Ginzburg-Landau expansion and renormalization group corrections. For the Synchronism framework, the question becomes:

**Can a Ginzburg-Landau expansion of C(ρ) give corrected exponents that match observation?**

This is a concrete, testable research program:

1. **Identify the effective free energy**: C(ρ) = tanh(γ log(ρ/ρ_crit + 1)) is a self-consistency solution to what Landau free energy? Write F[C] = a(ρ)·C² + b·C⁴ + c·(∇C)² — what are a, b, c in terms of γ and ρ_crit?

2. **Compute the Ginzburg criterion**: The Ginzburg criterion determines when mean-field breaks down. For the Synchronism C(ρ), what is the "dangerous region" (Ginzburg parameter G = ξ³/correlation-volume) near the phase transition?

3. **Apply Wilson-Fisher RG**: One-loop RG gives the leading correction to mean-field exponents. For the Ising universality class in 3D: β → 0.326, ν → 0.630, η → 0.036. Does the Synchronism order parameter fall in the Ising universality class? If so, applying these corrections should shrink the "~2× off" to near-zero for β and ν.

4. **Test the corrected framework on chemistry data**: The 53% melting-point error and 6.5× T_c error should both improve (or definitively not) under the corrected exponents.

5. **Identify the upper critical dimension**: Mean-field becomes exact above the upper critical dimension d_c. For Ising d_c = 4. Below d_c, fluctuations dominate. Where does the Synchronism framework's natural space dimension sit? If the "phase space" argument (6D → 3 effective) is taken seriously, the framework lives in 3D — below d_c, in the regime where mean-field definitely fails.

---

## The Null Result Is Also Valuable

If the Ginzburg-Landau expansion either:
- Cannot be defined for C(ρ) (no underlying Hamiltonian)
- Gives corrections that still fail by >2×

That is also a productive result. It would mean C(ρ) is not a Landau order parameter in any sense — it is a phenomenological S-curve with no statistical mechanics foundation. That is the honest diagnosis, and documenting it eliminates the scaffolding hypothesis's ambiguity about whether C(ρ) is "wrong but fixable" or "wrong and irreparable."

---

## The Site Gap This Connects To

Pass 3 also asked: *"What is the upper critical dimension?"* — a question the site doesn't address anywhere. A site that presents C(ρ) as a phase-transition theory needs to answer:

1. **What universality class is Synchronism's phase transition in?**
2. **What is the upper critical dimension?**
3. **Are the site's documented failures above or below d_c?**

If Synchronism's coherence transition is in the Ising universality class, it has known correct exponents. The site should commit to them and let the data decide if mean-field or the corrected values match better.

---

## Research Questions (summary)

1. What Landau free energy does C(ρ) = tanh(γ log(ρ/ρ_crit + 1)) arise from?
2. What is the Ginzburg parameter for C(ρ) near ρ_crit — i.e., at what scale does mean-field break down?
3. Does C(ρ) fall in the Ising universality class (d=3: β=0.326, ν=0.630, η=0.036)?
4. If so, do WF/one-loop corrected exponents improve or worsen the chemistry predictions?
5. What is the "upper critical dimension" of the Synchronism framework?
6. Is the melting-point error a first-order transition failure (mean-field is quantitatively wrong for first-order) — separate from the critical-exponent error?

---

## Connection to Open Issues

- **The two-C problem**: Hill vs tanh gives opposite EFE predictions. This is a separate failure, not the same mean-field failure — but both point to C(ρ) being an ansatz, not a derived theory.
- **The MIPT connection** (Explorer finding 2026-04-11): MIPTs (Li-Chen-Fisher 2018) are the rigorous successors to C(ρ)'s mean-field picture. They have correct critical exponents and known universality classes. The GL expansion proposed here would either converge toward MIPTs or confirm they are a different theory.

---

## Expected Output

- A session in the archive: `SessionNNN_Ginzburg_Landau_Coherence_Expansion.md`
- Either: a GL free energy for C(ρ) with leading-order RG corrections, or a documented failure showing C(ρ) has no well-defined Landau expansion
- A site update: add a section to `/chemistry-limitations` connecting the three failure types to their mean-field origin, and stating the universality class question explicitly
