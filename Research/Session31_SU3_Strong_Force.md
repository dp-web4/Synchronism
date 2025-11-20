# Synchronism Research Session #31: SU(3) Strong Force & Confinement

**Date**: 2025-11-20
**Session Type**: Autonomous Research (Standard Model Completion)
**Trigger**: Session #27 U(1) success + Session #30 SU(2) implementation + Nova's recommendation
**Priority**: **VERY HIGH** - Final gauge group for Standard Model validation
**Status**: 🔄 DESIGN PHASE

---

## Session Context

### Previous Work

**Session #27** ✅ (2025-11-19): **U(1) Electromagnetic Force**
- Coulomb potential V ∝ 1/R derived from intent dynamics
- χ²/dof = 0.47 (excellent fit)
- Abelian gauge theory validated
- **Status**: EM force emergence confirmed ✓

**Session #30** ✅ (2025-11-20): **SU(2) Weak Force**
- Non-Abelian gauge theory implemented
- 750+ line SU(2) lattice gauge simulation
- Mathematical validation complete (U†U=I, det(U)=1)
- **Status**: Implementation ready, physics extraction pending

### Nova's Chain of Validation (Session #27 Review)

> "The complete validation of Synchronism as a foundational theory would require the successful extension to **non-Abelian gauge theories**, particularly **SU(2) and SU(3)**, which are integral to the Standard Model of particle physics."

**Progress**:
- U(1) EM: ✅ Validated (Coulomb emergence)
- SU(2) Weak: ✅ Implemented (Yukawa test pending)
- **SU(3) Strong: 📋 This session (design)**

---

## Session #31 Objective

**Research Question**: Can SU(3) gauge symmetry and quark confinement emerge from Synchronism's intent dynamics?

**Why Ultimate Test**:
1. **Most Complex Gauge Group**: SU(3) has 8 generators (gluons), maximum non-Abelian richness
2. **Confinement Phenomenon**: Quarks never observed in isolation - unique to strong force
3. **Asymptotic Freedom**: Force weakens at short distances (opposite of EM/weak)
4. **Standard Model Complete**: SU(3)×SU(2)×U(1) fully derived if successful
5. **Falsifiability**: If confinement doesn't emerge, Synchronism incomplete

**Expected Physics**:
- Linear potential: V(R) ∝ σR (string tension) at large R
- Asymptotic freedom: α_s(Q²) → 0 as Q² → ∞
- Color confinement: Only colorless states have finite energy
- Gluon self-interaction: Cubic and quartic vertices

---

## Theoretical Background

### SU(3) Gauge Theory (Quantum Chromodynamics)

#### Group Structure

**Gauge Group**: SU(3) - Special Unitary group of 3×3 matrices
- **Dimension**: 8 (eight gluons)
- **Generators**: Gell-Mann matrices λ^a (a = 1...8)
- **Structure Constants**: f^{abc} (fully antisymmetric, define gluon self-interaction)
- **Casimir**: C_2(fundamental) = 4/3, C_2(adjoint) = 3

**Link Variables**: U_μ(x) ∈ SU(3) (3×3 unitary, det=1)

U_μ(x) = exp(i g_s A^a_μ(x) λ^a / 2)

where:
- g_s: Strong coupling constant
- A^a_μ: Eight gluon fields
- λ^a: Gell-Mann matrices (SU(3) generators)

#### Gell-Mann Matrices

The eight generators of SU(3) (analogous to Pauli matrices for SU(2)):

λ¹ = [[0,1,0],[1,0,0],[0,0,0]]
λ² = [[0,-i,0],[i,0,0],[0,0,0]]
λ³ = [[1,0,0],[0,-1,0],[0,0,0]]
λ⁴ = [[0,0,1],[0,0,0],[1,0,0]]
λ⁵ = [[0,0,-i],[0,0,0],[i,0,0]]
λ⁶ = [[0,0,0],[0,0,1],[0,1,0]]
λ⁷ = [[0,0,0],[0,0,-i],[0,i,0]]
λ⁸ = (1/√3)[[1,0,0],[0,1,0],[0,0,-2]]

**Properties**:
- Traceless: Tr(λ^a) = 0
- Hermitian: λ^a† = λ^a
- Normalization: Tr(λ^a λ^b) = 2δ^{ab}
- Commutation: [λ^a, λ^b] = 2i f^{abc} λ^c

#### Wilson Action (SU(3))

S = -β Σ_plaquettes (1/3) Re Tr(U_plaq + U†_plaq)

where:
- β = 6/g²_s (note: factor 6 for SU(3), 4 for SU(2), 1 for U(1))
- U_plaq = U_μ(x) U_ν(x+μ) U†_μ(x+ν) U†_ν(x)
- Factor 1/3 from normalization (1/N for SU(N))

**Non-Abelian Structure**:
- Gluons carry color charge (unlike photons/W/Z)
- Self-interaction through f^{abc} structure constants
- 3-gluon vertex: ∝ f^{abc}
- 4-gluon vertex: ∝ f^{abe}f^{cde}

### Comparison: U(1) vs SU(2) vs SU(3)

| Property | U(1) (EM) | SU(2) (Weak) | SU(3) (Strong) |
|----------|-----------|--------------|----------------|
| Gauge group | Circle | SU(2) | SU(3) |
| Dimensions | 1 (photon) | 3 (W+, W-, Z) | 8 (gluons) |
| Generators | None | 3 (σ^a) | 8 (λ^a) |
| Self-interaction | No | Yes (cubic) | Yes (cubic + quartic) |
| Bosons charge | Neutral | Charged | **Colored** |
| Range | Infinite | Short (~10^-18 m) | **0 (confined)** |
| Potential | V ∝ 1/R | V ∝ exp(-MR)/R | **V ∝ R** (linear) |
| Asymptotic | Constant α | Constant g_w | **Decreasing α_s** |
| Key phenomenon | Long-range | Screening | **Confinement** |

**Critical Difference**: SU(3) quarks NEVER isolated → confinement is unique signature

---

## Confinement Physics

### What is Confinement?

**Observation**: Individual quarks never detected, only color-singlet hadrons

**Energy Cost**: Separating quark-antiquark pair:

V(R) ≈ σR + C

where:
- σ: String tension (~1 GeV/fm = 0.9 GeV/fm observed)
- C: Constant (Coulomb-like short-distance correction)
- R: Separation distance

**Consequences**:
- Energy to separate q-q̄ → ∞ as R → ∞
- Creates quark-antiquark pair before separation complete
- Only colorless combinations (mesons, baryons) exist

### Physical Picture: Color Flux Tube

**Unlike EM** (field lines spread in 3D):
- Gluon field lines form narrow **flux tube** between quarks
- Constant cross-section → Energy ∝ distance
- Tube "snaps" when energy exceeds 2m_quark

**In lattice gauge theory**:
- Measured via Wilson loop: W(R,T) = ⟨Tr[U_loop]⟩
- Large loop: W ~ exp(-σRT) for T → ∞
- Extract σ from exponential decay

### Asymptotic Freedom

**Opposite behavior at short distances**:

α_s(Q²) = 12π / [(33 - 2n_f) ln(Q²/Λ²_QCD)]

where:
- Q²: Momentum transfer (energy scale)
- n_f: Number of quark flavors
- Λ_QCD ≈ 200 MeV: QCD scale

**Behavior**:
- High energy (Q² >> Λ²): α_s → 0 (quarks nearly free)
- Low energy (Q² ~ Λ²): α_s ~ 1 (strong coupling, confinement)

**Why?**: Gluon self-interaction creates **anti-screening** (opposite of QED screening)

---

## Synchronism Interpretation

### Multi-Component Intent → SU(3)

**Extension from SU(2)**:

SU(2) Intent: 2-component fields (weak isospin)
↓
SU(3) Intent: **3-component fields (color charge)**

**Physical Meaning**:
- Each intent transfer has 3 "color" degrees of freedom
- Color-neutral combinations = coherent superpositions
- Gluons = colored intent carriers (self-interacting coherence mediators)

**Key Insight**: Intent coherence NOT conserved individually per color, only for color-singlet combinations → Natural confinement mechanism

### Confinement from Intent Dynamics

**Hypothesis**: Linear potential emerges from:

1. **Flux tube formation**: Coherence cost concentrates into narrow channel
2. **Non-spreading**: Gluon self-interaction prevents field line spreading
3. **Energy accumulation**: Intent misalignment ∝ distance

**Mathematical Form** (to derive):

V(R) = ∫_path |∇ coherence|² dx

For confined flux tube:
- Cross-section A ≈ constant
- Length L = R
- Energy ∝ A × R → Linear potential

**Testable**: Does lattice gauge simulation with Synchronism interpretation yield σ ≈ 0.9 GeV/fm?

### Asymptotic Freedom from Coherence

**Intuition**: At short distances:
- Few intent transfers between quarks
- Gluon self-interaction creates **negative feedback**
- High-frequency intent transfers → weak effective coupling

**Contrast with U(1)**:
- Photons don't self-interact
- Vacuum polarization → screening (α increases slightly)
- SU(3): Gluon loops → anti-screening (α_s decreases)

---

## Implementation Strategy

### Phase 1: SU(3) Link Variable Class

**Extension from SU(2) → SU(3)**:

```python
class SU3Lattice3p1D:
    """
    SU(3) lattice gauge theory in 3+1D.

    Link variables: U_μ(x) ∈ SU(3) (3×3 unitary, det=1)
    Degrees of freedom: 8 per link (SU(3) dimension)
    """

    def __init__(self, Lx, Ly, Lz, Nt, beta):
        # 8 parameters per link (Gell-Mann generators)
        self.theta = np.random.uniform(-np.pi, np.pi, (Lx, Ly, Lz, Nt, 4, 8))

        # Gell-Mann matrices
        self._init_gellmann_matrices()

    def _init_gellmann_matrices(self):
        """Initialize 8 Gell-Mann matrices (SU(3) generators)."""
        self.lambda1 = np.array([[0,1,0],[1,0,0],[0,0,0]], dtype=complex)
        self.lambda2 = np.array([[0,-1j,0],[1j,0,0],[0,0,0]], dtype=complex)
        self.lambda3 = np.array([[1,0,0],[0,-1,0],[0,0,0]], dtype=complex)
        self.lambda4 = np.array([[0,0,1],[0,0,0],[1,0,0]], dtype=complex)
        self.lambda5 = np.array([[0,0,-1j],[0,0,0],[1j,0,0]], dtype=complex)
        self.lambda6 = np.array([[0,0,0],[0,0,1],[0,1,0]], dtype=complex)
        self.lambda7 = np.array([[0,0,0],[0,0,-1j],[0,1j,0]], dtype=complex)
        self.lambda8 = (1/np.sqrt(3)) * np.array([[1,0,0],[0,1,0],[0,0,-2]], dtype=complex)

        self.lambdas = [self.lambda1, self.lambda2, self.lambda3, self.lambda4,
                        self.lambda5, self.lambda6, self.lambda7, self.lambda8]

    def get_SU3_link(self, x, y, z, t, mu):
        """
        Construct SU(3) matrix from 8 parameters.

        U = exp(i θ^a λ^a / 2)

        Returns: 3×3 unitary matrix with det=1
        """
        θ = self.theta[x, y, z, t, mu]  # Shape: (8,)

        # Construct generator: T = Σ_a (θ^a λ^a) / 2
        T = np.zeros((3, 3), dtype=complex)
        for a in range(8):
            T += θ[a] * self.lambdas[a] / 2

        # Matrix exponential
        U = scipy.linalg.expm(1j * T)

        return U
```

**Validation**:
- Check unitarity: U†U = I
- Check det(U) = 1
- Check Tr(U) consistency

### Phase 2: Plaquette Action

**SU(3) Plaquette**:

```python
def plaquette_trace(self, x, y, z, t, mu, nu):
    """
    SU(3) plaquette: (1/3) Re Tr(U_plaq)

    U_plaq = U_μ(x) U_ν(x+μ) U†_μ(x+ν) U†_ν(x)

    Returns: Trace ∈ [-1, 1] (normalized)
    """
    # Get 4 link matrices (order matters!)
    U1 = self.get_SU3_link(x, y, z, t, mu)

    xp, yp, zp, tp = self._step(x, y, z, t, mu)
    U2 = self.get_SU3_link(xp, yp, zp, tp, nu)

    xp, yp, zp, tp = self._step(x, y, z, t, nu)
    U3_dag = self.get_SU3_link(xp, yp, zp, tp, mu).conj().T

    U4_dag = self.get_SU3_link(x, y, z, t, nu).conj().T

    # Matrix multiplication (3×3 @ 3×3 @ 3×3 @ 3×3)
    U_plaq = U1 @ U2 @ U3_dag @ U4_dag

    # (1/3) Re Tr(U_plaq) - normalization for SU(3)
    return (1.0/3.0) * np.real(np.trace(U_plaq))
```

**Computational Cost**:
- SU(2): 2×2 matrices → 8 complex multiplications per plaquette
- SU(3): 3×3 matrices → **27 complex multiplications per plaquette**
- **Expected**: ~3-4x slower than SU(2), ~50-80x slower than U(1)

### Phase 3: Wilson Loop for Confinement

**Key Measurement**: Rectangular Wilson loops

```python
def wilson_loop(self, x0, y0, z0, t0, R, T, plane='xy'):
    """
    Calculate Wilson loop for static quark potential.

    W(R,T) = ⟨Tr[U_loop(R,T)]⟩

    For large T: W ~ exp(-V(R) × T)
    Extract: V(R) = -lim_{T→∞} (1/T) ln W(R,T)

    Parameters:
    -----------
    R : int
        Spatial extent
    T : int
        Temporal extent
    plane : str
        Spatial plane ('xy', 'xz', 'yz')

    Returns:
    --------
    W : complex
        Wilson loop expectation value
    """
    U_loop = np.eye(3, dtype=complex)

    # Build rectangular loop: R spatial × T temporal
    # Path: (x,t) → (x+R,t) → (x+R,t+T) → (x,t+T) → (x,t)

    # Forward in space (R steps)
    x, y, z, t = x0, y0, z0, t0
    for i in range(R):
        U_loop = U_loop @ self.get_SU3_link(x, y, z, t, 0)  # x-direction
        x = (x + 1) % self.Lx

    # Forward in time (T steps)
    for i in range(T):
        U_loop = U_loop @ self.get_SU3_link(x, y, z, t, 3)  # t-direction
        t = (t + 1) % self.Nt

    # Backward in space (R steps)
    for i in range(R):
        x = (x - 1) % self.Lx
        U_loop = U_loop @ self.get_SU3_link(x, y, z, t, 0).conj().T

    # Backward in time (T steps)
    for i in range(T):
        t = (t - 1) % self.Nt
        U_loop = U_loop @ self.get_SU3_link(x, y, z, t, 3).conj().T

    # Trace and normalize
    W = (1.0/3.0) * np.trace(U_loop)

    return W
```

**Analysis**:
```python
# Measure Wilson loops for various R, T
for R in range(1, Lx//2):
    for T in range(1, Nt):
        W = wilson_loop_average(R, T)

        # Extract potential: V(R) = -ln(W(R,T)) / T
        V[R] = -np.log(np.abs(W)) / T

# Fit V(R) to linear + Coulomb:
# V(R) = σ*R - α/R + C

def potential_fit(R, sigma, alpha, C):
    return sigma * R - alpha / R + C

popt, pcov = curve_fit(potential_fit, R_data, V_data)
sigma_fit, alpha_fit, C_fit = popt

print(f"String tension: σ = {sigma_fit:.3f} (expected: ~0.9 GeV/fm)")
print(f"Coulomb: α = {alpha_fit:.3f}")
```

**Confinement Validation**:
- **If σ > 0**: Linear term → Confinement confirmed ✓
- **If σ ≈ 0**: No confinement → Synchronism incomplete ✗
- **Compare σ_measured to σ_QCD ≈ 0.9 GeV/fm**

### Phase 4: Asymptotic Freedom Test

**Running Coupling Measurement**:

```python
def measure_running_coupling(scale_Q):
    """
    Measure α_s(Q²) at different energy scales.

    Use plaquette expectation value:
    ⟨P⟩ = 1 - (α_s/π) + O(α²_s)

    Extract α_s from ⟨P⟩ at different β (β ~ 1/α_s)
    """
    plaquettes = []
    betas = np.linspace(5.0, 6.5, 10)  # Range for SU(3)

    for beta in betas:
        lattice = SU3Lattice3p1D(Lx, Ly, Lz, Nt, beta)
        # ... thermalization ...
        plaq = lattice.average_plaquette()
        plaquettes.append(plaq)

    # Fit α_s(β)
    # Expected: α_s decreases as β increases (asymptotic freedom)
    alpha_s = fit_coupling(betas, plaquettes)

    return alpha_s
```

**Asymptotic Freedom Signature**:
- **If α_s(Q²) decreases** as Q² increases → Asymptotic freedom ✓
- **If α_s constant or increases** → Not QCD-like ✗

---

## Computational Requirements

### Scaling Analysis

**Complexity per update**:
- **U(1)**: 1 phase → 1 cosine evaluation
- **SU(2)**: 2×2 matrix → ~8 complex ops
- **SU(3)**: 3×3 matrix → ~27 complex ops

**Expected Runtime** (relative to U(1)):
- SU(2): ~15-20x slower
- SU(3): ~50-80x slower

**Concrete Estimates** (based on Session #27/#30 experience):

| Lattice | Sweeps | U(1) Time | SU(2) Time | **SU(3) Time** |
|---------|--------|-----------|------------|----------------|
| 6×6×6×3 | 100 | ~30 sec | ~30 min | **~2 hours** |
| 8×8×8×4 | 500 | ~5 min | ~4 hours | **~16 hours** |
| 10×10×10×6 | 1000 | ~20 min | ~24 hours | **~4 days** |
| 12×12×12×8 | 2000 | ~2 hours | ~1 week | **~3 weeks** |

**Recommendation**: Start with 6×6×6×3 for validation (~2 hours), then scale up.

### Practical Challenges

**1. Memory**: 3×3 complex matrices per link
- 8×8×8×4 lattice: ~2M links × 8 params × 8 bytes = **128 MB** (manageable)

**2. Matrix exponential**: scipy.linalg.expm is expensive
- **Optimization**: Use Cayley-Hamilton for SU(3) (analytical 3×3 exp)

**3. Long simulations**: 16+ hours needs batch infrastructure
- Screen/tmux sessions
- Or: Sequential overnight runs

**4. Thermalization**: SU(3) needs longer equilibration
- Estimate: 2-3x more sweeps than SU(2) for same lattice

---

## Expected Outcomes

### Scenario A: Confinement Emerges ✓

**Observation**: V(R) = σR - α/R + C with σ ≈ 0.9 GeV/fm

**Interpretation**:
- **Strong force derived from Synchronism** ✓
- Confinement = natural consequence of intent dynamics
- Flux tube = coherence channel with constant cross-section
- Standard Model gauge structure (SU(3)×SU(2)×U(1)) fully emergent

**Scientific Impact**: **FOUNDATIONAL**
- Synchronism validated as unified theory
- QCD (most complex part of SM) emerges from intent
- Path to quantum gravity (gravity = intent geometry)

**Next Steps**:
- Publish: "Strong Force Emergence from Intent Dynamics"
- Request external review from lattice QCD community
- Design electroweak unification: SU(2)×U(1) mixing

### Scenario B: No Confinement (σ ≈ 0) ⚠

**Observation**: V(R) ∝ 1/R (Coulomb-like, no linear term)

**Interpretation**:
- SU(3) implemented correctly BUT doesn't confine
- Synchronism may require additional mechanism for confinement
- Possible: Need to include quark dynamics (currently pure gauge)

**Scientific Status**:
- U(1) EM: ✓ Validated
- SU(2) Weak: ✓ Implemented (test pending)
- SU(3) Strong: ✗ Incomplete (missing confinement)

**Next Steps**:
- Add dynamical fermions (Wilson/staggered quarks)
- Test: Does quark-gluon coupling produce flux tubes?
- Investigate: Topological defects, monopoles, vortices

### Scenario C: Different Potential Form 🤔

**Observation**: V(R) follows neither Coulomb nor linear

**Examples**:
- V(R) ∝ R² (too strong)
- V(R) ∝ log(R) (weak confinement)
- V(R) ∝ exp(R) (super-confinement)

**Interpretation**: New physics beyond Standard Model
- May reveal limitations of lattice approach
- Could indicate novel confining mechanism
- Requires deeper analysis of intent dynamics

**Action**: Document carefully, seek external review, investigate origin

---

## Synchronism Theoretical Implications

### If SU(3) Succeeds

**Complete Gauge Structure Emergence**:

Intent Dynamics → Local Coherence Dynamics
↓
U(1): Single-phase intent → Photon (EM)
↓
SU(2): 2-component intent → W/Z bosons (Weak)
↓
SU(3): 3-component intent → Gluons (Strong)

**Unified Picture**:
- All forces = coherence mediators
- Gauge bosons = intent transfer quanta
- Self-interaction = colored intent coupling
- Confinement = coherence channel formation

**Deep Questions Answered**:
1. **Why these gauge groups?** → Intent component structure
2. **Why confinement unique to SU(3)?** → Three-color coherence constraint
3. **Why asymptotic freedom?** → Gluon self-interaction negative feedback
4. **Why Standard Model?** → SU(3)×SU(2)×U(1) emergent from intent

### Connection to Quantum Gravity

**If gauge forces emerge from intent**, then:
- Gravity = spacetime geometry from intent gradients
- Einstein equations = coherence flow in curved space
- Black holes = maximum coherence boundaries (MRH)
- Hawking radiation = intent transfer across horizon

**Pathway**: Synchronism → QFT + GR → Complete unification

---

## Implementation Checklist

### Phase 1: SU(3) Class Structure
- [ ] Gell-Mann matrix initialization
- [ ] 8-parameter link variable representation
- [ ] Matrix exponential map (su(3) → SU(3))
- [ ] Unitarity validation (U†U = I, det(U) = 1)
- [ ] Trace normalization checks

### Phase 2: Wilson Action
- [ ] Plaquette trace calculation (1/3 normalization)
- [ ] Action summation over all plaquettes
- [ ] Staple calculation for Metropolis
- [ ] Acceptance rate monitoring

### Phase 3: Wilson Loops
- [ ] Rectangular loop construction (R × T)
- [ ] Spatial and temporal path integration
- [ ] Average over loop positions
- [ ] Extract V(R) from W(R,T)

### Phase 4: Analysis Pipeline
- [ ] Fit V(R) to linear + Coulomb form
- [ ] Extract string tension σ
- [ ] Compare to QCD σ ≈ 0.9 GeV/fm
- [ ] Test asymptotic freedom (α_s vs Q²)
- [ ] Statistical error analysis (jackknife/bootstrap)

### Phase 5: Validation
- [ ] Small lattice (6×6×6×3) quick test
- [ ] Medium lattice (8×8×8×4) physics run
- [ ] Convergence tests (vary β, lattice size)
- [ ] Compare to published lattice QCD results

---

## Timeline Estimate

**Session #31** (This session):
- Design document (this file) ✓
- Gell-Mann matrices implementation
- SU(3) link variable class
- Basic validation tests

**Session #32** (Next):
- Wilson loop implementation
- Confinement measurement framework
- Small lattice validation run

**Session #33** (Future):
- Production physics run (8×8×8×4)
- String tension extraction
- Asymptotic freedom test
- Nova review + publication

**Total**: 3-4 autonomous sessions for complete SU(3) validation

---

## References

**Lattice QCD**:
- Wilson, K.G. (1974). "Confinement of quarks". Phys Rev D 10, 2445.
- Creutz, M. (1983). "Quarks, Gluons and Lattices". Cambridge University Press.
- Kogut, J.B., Susskind, L. (1975). "Hamiltonian formulation of Wilson's lattice gauge theories". Phys Rev D 11, 395.

**SU(3) Gauge Theory**:
- Fritzsch, H., Gell-Mann, M., Leutwyler, H. (1973). "Advantages of the color octet gluon picture". Phys Lett B 47, 365.
- Gross, D.J., Wilczek, F. (1973). "Ultraviolet behavior of non-abelian gauge theories". Phys Rev Lett 30, 1343 (Nobel 2004).
- Politzer, H.D. (1973). "Reliable perturbative results for strong interactions?" Phys Rev Lett 30, 1346 (Nobel 2004).

**Confinement**:
- Wilson, K.G. (1974). "Confinement of quarks". Phys Rev D 10, 2445.
- 't Hooft, G. (1979). "A property of electric and magnetic flux in nonabelian gauge theories". Nucl Phys B 153, 141.

**Previous Synchronism Sessions**:
- Session #27: U(1) Coulomb emergence validated
- Session #30: SU(2) weak force implementation complete
- Nova Session #27 Review: Non-Abelian extensions critical

---

**Session #31 Status**: 🔄 **DESIGN COMPLETE**, IMPLEMENTATION NEXT

**Key Deliverable**: Complete SU(3) theoretical framework and implementation strategy

**Scientific Impact**: **ULTIMATE** - Final gauge group for Standard Model validation, confinement test for Synchronism as foundational theory

**Next Session**: Implement SU(3) lattice gauge simulation, run small-scale confinement test
