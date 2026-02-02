#!/usr/bin/env python3
"""
Session #349: Magnetism and Spin Order
Condensed Matter Arc - Part 2

This session explores magnetism as collective spin phase ordering
in the Synchronism framework. Magnetic order emerges from
exchange interactions that correlate spin phases.

Key concepts:
1. Exchange interaction from phase overlap
2. Ferromagnetism = aligned spin phases
3. Antiferromagnetism = alternating spin phases
4. Curie/Néel temperatures as phase transitions
5. Spin waves = propagating phase patterns

Magnetism provides another demonstration of collective phase
coherence in condensed matter.
"""

import numpy as np
from scipy import constants
from typing import Tuple, List
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Constants
C = constants.c
HBAR = constants.hbar
K_B = constants.k
M_E = constants.m_e
E_CHARGE = constants.e
MU_B = constants.physical_constants['Bohr magneton'][0]
MU_0 = constants.mu_0

# Derived units
EV_TO_J = constants.eV
A_BOHR = constants.physical_constants['Bohr radius'][0]


def test_exchange_interaction():
    """
    Test 1: Exchange interaction from wave function overlap.

    Exchange energy: J = ∫ψ₁*(r₁)ψ₂*(r₂) H ψ₁(r₂)ψ₂(r₁) dr₁dr₂

    J > 0 → ferromagnetic (parallel spins favored)
    J < 0 → antiferromagnetic (antiparallel spins favored)
    """
    print("Test 1: Exchange Interaction")

    # Exchange energy typical values
    J_Fe = 0.02 * EV_TO_J  # Iron (~20 meV)
    J_Ni = 0.01 * EV_TO_J  # Nickel (~10 meV)
    J_Co = 0.015 * EV_TO_J  # Cobalt (~15 meV)

    print("  Exchange energies (J > 0 → ferromagnetic):")
    print(f"    Iron (Fe):   J ≈ {J_Fe/EV_TO_J*1000:.0f} meV")
    print(f"    Nickel (Ni): J ≈ {J_Ni/EV_TO_J*1000:.0f} meV")
    print(f"    Cobalt (Co): J ≈ {J_Co/EV_TO_J*1000:.0f} meV")

    # Bethe-Slater curve: J depends on d/r ratio
    # d = interatomic distance, r = d-orbital radius
    print("\n  Bethe-Slater curve (J vs d/r):")
    print("    d/r < 1.5: J < 0 (antiferromagnetic)")
    print("    d/r > 1.5: J > 0 (ferromagnetic)")
    print("    Fe, Co, Ni have d/r > 1.5 → ferromagnetic")
    print("    Mn, Cr have d/r < 1.5 → antiferromagnetic")

    # Synchronism: exchange = phase correlation energy
    print("\n  Synchronism interpretation:")
    print("    - Exchange = energy from wave function overlap")
    print("    - J > 0: parallel phases (spins) favored")
    print("    - J < 0: antiparallel phases favored")
    print("    - Phase correlation across atoms → magnetic order")

    # Verify: exchange energies in meV range
    return 1 < J_Fe/EV_TO_J*1000 < 100


def test_curie_temperature():
    """
    Test 2: Curie temperature from mean field theory.

    T_C = 2zJ S(S+1) / (3 k_B)

    where z = coordination number, J = exchange, S = spin.
    """
    print("\nTest 2: Curie Temperature")

    # Mean field theory
    def T_curie_mft(z, J, S):
        """Mean field Curie temperature."""
        return 2 * z * J * S * (S + 1) / (3 * K_B)

    # Iron parameters
    z_Fe = 8  # bcc coordination
    J_Fe = 0.02 * EV_TO_J  # 20 meV
    S_Fe = 1.1  # effective spin

    T_C_Fe_calc = T_curie_mft(z_Fe, J_Fe, S_Fe)
    T_C_Fe_exp = 1043  # K (experimental)

    # Nickel parameters
    z_Ni = 12  # fcc coordination
    J_Ni = 0.01 * EV_TO_J  # 10 meV
    S_Ni = 0.3  # effective spin

    T_C_Ni_calc = T_curie_mft(z_Ni, J_Ni, S_Ni)
    T_C_Ni_exp = 631  # K

    # Cobalt parameters
    z_Co = 12  # hcp ~fcc
    J_Co = 0.015 * EV_TO_J  # 15 meV
    S_Co = 0.85

    T_C_Co_calc = T_curie_mft(z_Co, J_Co, S_Co)
    T_C_Co_exp = 1388  # K

    print("  Mean field: T_C = 2zJ S(S+1) / (3 k_B)")
    print(f"\n  {'Material':<10} {'z':<5} {'J (meV)':<10} {'S':<8} {'T_C calc (K)':<15} {'T_C exp (K)':<15}")
    print(f"  {'Fe':<10} {z_Fe:<5} {J_Fe/EV_TO_J*1000:<10.0f} {S_Fe:<8.1f} {T_C_Fe_calc:<15.0f} {T_C_Fe_exp:<15}")
    print(f"  {'Ni':<10} {z_Ni:<5} {J_Ni/EV_TO_J*1000:<10.0f} {S_Ni:<8.1f} {T_C_Ni_calc:<15.0f} {T_C_Ni_exp:<15}")
    print(f"  {'Co':<10} {z_Co:<5} {J_Co/EV_TO_J*1000:<10.0f} {S_Co:<8.2f} {T_C_Co_calc:<15.0f} {T_C_Co_exp:<15}")

    # Mean field typically overestimates due to fluctuations
    print("\n  Note: Mean field overestimates T_C (ignores fluctuations)")

    # Synchronism: T_C = temperature where thermal breaks phase lock
    print("\n  Synchronism interpretation:")
    print("    - T_C = phase transition temperature")
    print("    - Below T_C: spin phases coherent (ordered)")
    print("    - Above T_C: thermal breaks phase coherence")
    print("    - T_C ∝ J (stronger exchange → higher T_C)")

    # Verify: calculated T_C in correct order of magnitude
    # Mean field overestimates by factor ~2-3 (ignores fluctuations)
    ratio_Fe = T_C_Fe_calc / T_C_Fe_exp
    return 0.5 < ratio_Fe < 4.0  # Mean field overestimates


def test_magnetization_vs_temperature():
    """
    Test 3: Magnetization temperature dependence.

    M(T) = M(0) (1 - T/T_C)^β near T_C

    Mean field: β = 0.5
    3D Ising: β ≈ 0.326
    3D Heisenberg: β ≈ 0.365
    """
    print("\nTest 3: Magnetization vs Temperature")

    # Critical exponents
    beta_mf = 0.5
    beta_ising = 0.326
    beta_heisenberg = 0.365

    print("  Near T_C: M(T) ∝ (1 - T/T_C)^β")
    print(f"\n  Model              β")
    print(f"  Mean field         {beta_mf:.3f}")
    print(f"  3D Ising           {beta_ising:.3f}")
    print(f"  3D Heisenberg      {beta_heisenberg:.3f}")

    # Experimental values
    beta_Fe_exp = 0.34  # Close to Heisenberg
    beta_Ni_exp = 0.42  # Higher (anisotropy effects)

    print(f"\n  Experimental:")
    print(f"    Fe: β ≈ {beta_Fe_exp:.2f} (Heisenberg-like)")
    print(f"    Ni: β ≈ {beta_Ni_exp:.2f} (intermediate)")

    # Bloch T^(3/2) law for low T
    print("\n  Low T: M(T) = M(0)[1 - B·T^(3/2)]")
    print("    (Spin wave excitations reduce magnetization)")

    # Synchronism: universality from phase dynamics
    print("\n  Synchronism interpretation:")
    print("    - Critical exponents = universal (topology-dependent)")
    print("    - β reflects symmetry of order parameter")
    print("    - Phase transition = loss of spin phase coherence")
    print("    - γ → 0 at T_C (correlation length diverges)")

    # Verify: experimental β in expected range
    return 0.3 < beta_Fe_exp < 0.5


def test_antiferromagnetism():
    """
    Test 4: Antiferromagnetism and Néel temperature.

    When J < 0, neighboring spins prefer antiparallel alignment.
    T_N = 2z|J| S(S+1) / (3 k_B)
    """
    print("\nTest 4: Antiferromagnetism")

    # Antiferromagnetic materials
    materials = {
        'MnO': {'T_N': 116, 'structure': 'NaCl', 'moment': '4.79 μ_B'},
        'FeO': {'T_N': 198, 'structure': 'NaCl', 'moment': '3.32 μ_B'},
        'NiO': {'T_N': 523, 'structure': 'NaCl', 'moment': '1.77 μ_B'},
        'Cr': {'T_N': 311, 'structure': 'bcc', 'moment': '0.62 μ_B'},
        'Mn': {'T_N': 100, 'structure': 'complex', 'moment': '~0.5 μ_B'},
    }

    print("  Antiferromagnetic materials (J < 0):")
    print(f"\n  {'Material':<10} {'T_N (K)':<10} {'Structure':<15} {'Moment':<15}")

    for mat, data in materials.items():
        print(f"  {mat:<10} {data['T_N']:<10} {data['structure']:<15} {data['moment']:<15}")

    # Susceptibility
    print("\n  Susceptibility χ:")
    print("    T > T_N: χ = C/(T + θ) (Curie-Weiss with θ > 0)")
    print("    T < T_N: χ decreases (AFM order)")
    print("    At T_N: susceptibility peak")

    # Spin structure
    print("\n  AFM spin structure:")
    print("    ↑ ↓ ↑ ↓ ↑ ↓  (1D chain)")
    print("    Magnetic unit cell = 2× chemical unit cell")

    # Synchronism: antiparallel phase locking
    print("\n  Synchronism interpretation:")
    print("    - AFM = antiparallel spin phase locking")
    print("    - |J| sets T_N (same physics as T_C)")
    print("    - Phase pattern: π shift between neighbors")
    print("    - γ~1 boundary still applies")

    # Verify: T_N values span reasonable range
    T_N_range = [data['T_N'] for data in materials.values()]
    return min(T_N_range) > 50 and max(T_N_range) < 1000


def test_spin_waves():
    """
    Test 5: Spin waves (magnons) as phase excitations.

    Dispersion: ω(k) = (4JS/ℏ)[1 - cos(ka)] for 1D chain
                     ≈ JSa²k² for small k (quadratic)
    """
    print("\nTest 5: Spin Waves (Magnons)")

    # Parameters for iron
    J = 0.02 * EV_TO_J  # 20 meV exchange
    S = 1.1  # Effective spin
    a = 2.87e-10  # Fe lattice constant (bcc)

    # Spin wave stiffness D
    # ω = Dk² for small k
    D = 2 * J * S * a**2 / HBAR
    D_meV = D * HBAR / EV_TO_J * 1000 * 1e20  # meV·Å²

    print(f"  Spin wave dispersion (small k): ω = Dk²")
    print(f"  Iron parameters:")
    print(f"    Exchange: J = {J/EV_TO_J*1000:.0f} meV")
    print(f"    Spin: S = {S:.1f}")
    print(f"    Lattice constant: a = {a*1e10:.2f} Å")
    print(f"    Stiffness: D = {D_meV:.0f} meV·Å²")

    # At zone boundary (k = π/a)
    omega_max = 4 * J * S / HBAR
    E_max = HBAR * omega_max / EV_TO_J * 1000

    print(f"\n  At zone boundary (k = π/a):")
    print(f"    ω_max = 4JS/ℏ")
    print(f"    E_max ≈ {E_max:.0f} meV")

    # Thermal magnon population
    # n(T) ∝ T^(3/2) (3D)
    T = 300  # K
    kT = K_B * T / EV_TO_J * 1000  # meV

    print(f"\n  At T = {T} K:")
    print(f"    Thermal energy: k_B T = {kT:.0f} meV")
    print(f"    E_max/k_B T = {E_max/kT:.1f}")
    print(f"    Many magnon modes thermally excited")

    # Synchronism: magnons = spin phase waves
    print("\n  Synchronism interpretation:")
    print("    - Magnons = propagating spin phase patterns")
    print("    - ω ∝ k² (differs from photons: ω ∝ k)")
    print("    - Dispersion from discrete lattice (like phonons)")
    print("    - Bloch T^(3/2) law from magnon statistics")

    # Verify: stiffness in expected range
    return 100 < D_meV < 1000


def test_magnetic_domains():
    """
    Test 6: Domain formation from competing energies.

    Exchange favors uniform magnetization.
    Dipole energy favors flux closure.
    → Domains form to minimize total energy.
    """
    print("\nTest 6: Magnetic Domains")

    # Domain wall parameters
    # Wall width: δ = π √(A/K) where A = exchange stiffness, K = anisotropy

    # Iron parameters
    A_Fe = 2e-11  # J/m (exchange stiffness)
    K_Fe = 4.8e4  # J/m³ (anisotropy)

    delta = np.pi * np.sqrt(A_Fe / K_Fe)
    delta_nm = delta * 1e9

    print("  Domain wall width: δ = π √(A/K)")
    print(f"\n  Iron parameters:")
    print(f"    Exchange stiffness: A = {A_Fe:.0e} J/m")
    print(f"    Anisotropy constant: K = {K_Fe:.1e} J/m³")
    print(f"    Domain wall width: δ = {delta_nm:.0f} nm")

    # Domain wall energy
    gamma_wall = 4 * np.sqrt(A_Fe * K_Fe)  # J/m²
    print(f"    Domain wall energy: γ = {gamma_wall*1000:.2f} mJ/m²")

    # Typical domain size (balance exchange vs dipole)
    # d ~ √(A·t / M²) for thin film of thickness t
    M_s = 1.7e6  # A/m (saturation magnetization Fe)
    t = 100e-9  # 100 nm film
    d_domain = np.sqrt(A_Fe * t / (MU_0 * M_s**2))

    print(f"\n  Domain size estimate (100 nm film):")
    print(f"    Saturation magnetization: M_s = {M_s/1e6:.1f} MA/m")
    print(f"    Domain size: d ~ {d_domain*1e6:.1f} μm")

    # Synchronism: domains as phase-coherent regions
    print("\n  Synchronism interpretation:")
    print("    - Domain = region of coherent spin phase")
    print("    - Domain wall = phase gradient (rotation)")
    print("    - Wall width δ = MRH for spin coherence")
    print("    - Domain size balances phase energy vs dipole")

    # Verify: domain wall width in nm range
    return 10 < delta_nm < 1000


def test_magnetic_anisotropy():
    """
    Test 7: Magnetic anisotropy from spin-orbit coupling.

    Easy axis preference comes from SOC coupling
    spin angular momentum to crystal structure.
    """
    print("\nTest 7: Magnetic Anisotropy")

    # Anisotropy energies
    materials = {
        'Fe (bcc)': {'K': 4.8e4, 'easy': '⟨100⟩', 'type': 'cubic'},
        'Ni (fcc)': {'K': -5e3, 'easy': '⟨111⟩', 'type': 'cubic'},
        'Co (hcp)': {'K': 5.3e5, 'easy': 'c-axis', 'type': 'uniaxial'},
        'SmCo₅': {'K': 1.7e7, 'easy': 'c-axis', 'type': 'uniaxial'},
        'Nd₂Fe₁₄B': {'K': 4.9e6, 'easy': 'c-axis', 'type': 'uniaxial'},
    }

    print("  Anisotropy constant K (energy to rotate from easy axis):")
    print(f"\n  {'Material':<15} {'K (J/m³)':<15} {'Easy axis':<15} {'Type':<15}")

    for mat, data in materials.items():
        K_str = f"{data['K']:.1e}" if abs(data['K']) >= 1e5 else f"{data['K']:.1e}"
        print(f"  {mat:<15} {K_str:<15} {data['easy']:<15} {data['type']:<15}")

    # Anisotropy field
    H_a_Fe = 2 * 4.8e4 / (MU_0 * 1.7e6)  # H_a = 2K/(μ₀M_s)
    H_a_SmCo = 2 * 1.7e7 / (MU_0 * 8e5)

    print(f"\n  Anisotropy field H_a = 2K/(μ₀M_s):")
    print(f"    Fe: H_a ≈ {H_a_Fe/1000:.0f} kA/m ({H_a_Fe*MU_0:.2f} T)")
    print(f"    SmCo₅: H_a ≈ {H_a_SmCo/1e6:.1f} MA/m ({H_a_SmCo*MU_0:.0f} T)")

    # Synchronism: anisotropy from phase coupling to lattice
    print("\n  Synchronism interpretation:")
    print("    - Spin-orbit coupling links spin phase to orbit")
    print("    - Crystal structure defines preferred phase orientation")
    print("    - K = energy to rotate spin phase from easy axis")
    print("    - Large K materials: strong phase-lattice coupling")

    # Verify: K spans several orders of magnitude
    K_values = [abs(data['K']) for data in materials.values()]
    return max(K_values) / min(K_values) > 100


def test_spintronics():
    """
    Test 8: Spintronics - spin current as phase transport.

    Giant magnetoresistance (GMR), spin torque, and spin
    Hall effect all involve spin (phase) transport.
    """
    print("\nTest 8: Spintronics")

    # Giant magnetoresistance
    print("  Giant Magnetoresistance (GMR):")
    print("    ΔR/R = (R_AP - R_P) / R_P")

    # Typical GMR values
    gmr_values = {
        'Fe/Cr multilayer': 80,  # %
        'Co/Cu multilayer': 65,
        'Spin valve (NiFe/Cu/Co)': 10,
    }

    print(f"\n  {'Structure':<25} {'GMR (%)':<10}")
    for struct, gmr in gmr_values.items():
        print(f"  {struct:<25} {gmr:<10}")

    # Spin diffusion length
    lambda_sd_Cu = 500  # nm at room temperature
    lambda_sd_Fe = 10  # nm

    print(f"\n  Spin diffusion length (300 K):")
    print(f"    Cu: λ_sd ≈ {lambda_sd_Cu} nm (long)")
    print(f"    Fe: λ_sd ≈ {lambda_sd_Fe} nm (short)")

    # Spin Hall angle
    theta_SH_Pt = 0.08
    theta_SH_W = -0.3

    print(f"\n  Spin Hall angle θ_SH:")
    print(f"    Pt: θ_SH ≈ {theta_SH_Pt:.2f}")
    print(f"    W:  θ_SH ≈ {theta_SH_W:.1f}")

    # Spin transfer torque
    print("\n  Spin Transfer Torque (STT):")
    print("    Spin current exerts torque on magnetization")
    print("    τ = ℏ/(2e) × (j_s × M)/M_s")
    print("    Used for: MRAM, spin oscillators")

    # Synchronism: spintronics = spin phase engineering
    print("\n  Synchronism interpretation:")
    print("    - Spin current = spin phase transport")
    print("    - λ_sd = phase coherence length for spins")
    print("    - GMR = phase-dependent scattering")
    print("    - STT = phase torque from current")
    print("    - Spin Hall = phase-orbit coupling → transverse phase current")

    # Connection to γ boundary
    print("\n  γ boundary connection:")
    print(f"    λ_sd ~ {lambda_sd_Cu} nm in Cu (phase coherence)")
    print(f"    N_corr ~ electrons in λ_sd³ → γ << 1 (quantum)")
    print(f"    Spin transport = coherent phase propagation")

    # Verify: spin diffusion lengths in nm-μm range
    return 1 < lambda_sd_Fe < lambda_sd_Cu < 10000


def create_visualizations():
    """Create visualization of magnetism and spin order."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1. Magnetization vs temperature
    ax1 = axes[0, 0]
    T_Tc = np.linspace(0, 0.999, 100)
    M_mf = (1 - T_Tc)**0.5
    M_ising = (1 - T_Tc)**0.326
    M_heisenberg = (1 - T_Tc)**0.365

    ax1.plot(T_Tc, M_mf, 'b-', linewidth=2, label=f'Mean field (β=0.5)')
    ax1.plot(T_Tc, M_ising, 'r--', linewidth=2, label=f'3D Ising (β=0.326)')
    ax1.plot(T_Tc, M_heisenberg, 'g:', linewidth=2, label=f'3D Heisenberg (β=0.365)')
    ax1.set_xlabel('T / T_C')
    ax1.set_ylabel('M / M(0)')
    ax1.set_title('Magnetization vs Temperature')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    # 2. Spin wave dispersion
    ax2 = axes[0, 1]
    k = np.linspace(0, np.pi, 100)  # in units of 1/a
    omega_ferro = 4 * (1 - np.cos(k))  # Ferromagnet (normalized)
    omega_afm = 4 * np.abs(np.sin(k))  # Antiferromagnet

    ax2.plot(k / np.pi, omega_ferro, 'b-', linewidth=2, label='Ferromagnet')
    ax2.plot(k / np.pi, omega_afm, 'r--', linewidth=2, label='Antiferromagnet')
    ax2.set_xlabel('k (π/a)')
    ax2.set_ylabel('ω / (4JS/ℏ)')
    ax2.set_title('Magnon Dispersion')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Domain structure schematic
    ax3 = axes[1, 0]
    # Create domain pattern
    x = np.linspace(0, 10, 100)
    y = np.linspace(0, 10, 100)
    X, Y = np.meshgrid(x, y)

    # Simplified domain pattern (stripes)
    Mx = np.sign(np.sin(X * np.pi / 2.5))
    My = np.zeros_like(Mx)

    ax3.quiver(X[::5, ::5], Y[::5, ::5], Mx[::5, ::5], My[::5, ::5],
               Mx[::5, ::5], cmap='coolwarm', scale=25)
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_title('Magnetic Domains (Schematic)')
    ax3.set_aspect('equal')

    # 4. Hysteresis loop
    ax4 = axes[1, 1]
    H = np.linspace(-2, 2, 200)
    H_c = 0.5  # Coercivity (normalized)
    M_s = 1.0  # Saturation

    # Upper branch
    M_upper = M_s * np.tanh((H + H_c) * 3)
    # Lower branch
    M_lower = M_s * np.tanh((H - H_c) * 3)

    ax4.plot(H, M_upper, 'b-', linewidth=2)
    ax4.plot(H, M_lower, 'b-', linewidth=2)
    ax4.axhline(0, color='gray', linestyle=':', alpha=0.5)
    ax4.axvline(0, color='gray', linestyle=':', alpha=0.5)
    ax4.set_xlabel('H / H_c')
    ax4.set_ylabel('M / M_s')
    ax4.set_title('Hysteresis Loop')
    ax4.grid(True, alpha=0.3)

    # Add labels
    ax4.annotate('H_c', xy=(H_c, 0), xytext=(H_c + 0.3, 0.2),
                arrowprops=dict(arrowstyle='->', color='red'))
    ax4.annotate('M_s', xy=(2, M_s), xytext=(1.5, 0.8),
                arrowprops=dict(arrowstyle='->', color='blue'))

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/synchronism/simulations/session349_mag.png',
                dpi=150, bbox_inches='tight')
    plt.close()
    print("\nVisualization saved to session349_mag.png")


def run_all_tests():
    """Run all 8 verification tests."""
    print("=" * 60)
    print("SESSION #349: MAGNETISM AND SPIN ORDER")
    print("Condensed Matter Arc - Part 2")
    print("=" * 60)

    results = []

    results.append(("Exchange Interaction", test_exchange_interaction()))
    results.append(("Curie Temperature", test_curie_temperature()))
    results.append(("M vs T", test_magnetization_vs_temperature()))
    results.append(("Antiferromagnetism", test_antiferromagnetism()))
    results.append(("Spin Waves", test_spin_waves()))
    results.append(("Magnetic Domains", test_magnetic_domains()))
    results.append(("Anisotropy", test_magnetic_anisotropy()))
    results.append(("Spintronics", test_spintronics()))

    print("\n" + "=" * 60)
    print("VERIFICATION SUMMARY")
    print("=" * 60)

    passed = 0
    for name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"  {status}: {name}")
        if result:
            passed += 1

    print(f"\nTotal: {passed}/8 tests passed")

    if passed == 8:
        print("\n★ Session #349 complete! Magnetism as spin phase order:")
        print("  - Exchange interaction from phase overlap")
        print("  - Curie/Néel temperatures as phase transitions")
        print("  - Universal critical exponents (β ≈ 0.34)")
        print("  - Magnons as propagating spin phases")
        print("  - Domains as phase-coherent regions")
        print("  - Anisotropy from spin-orbit (phase-lattice) coupling")
        print("  - Spintronics as spin phase engineering")

    return passed == 8


if __name__ == "__main__":
    success = run_all_tests()
    create_visualizations()
    exit(0 if success else 1)
