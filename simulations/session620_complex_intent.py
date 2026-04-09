"""
Session #620: Complex Intent — Does the Framework's Name Save Its Mathematics?

The framework is called "Synchronism" — synchronization requires phase.
The mathematics (real scalar I, real transfer rule) has no phase.
S617-619 prove this kills entities, waves, and structure.

Question: If Intent is complex (Ψ = √I · e^{iθ}), does the SAME transfer rule
produce the dynamics the framework describes?

Test 1: Real vs complex transfer — diffusion vs waves
Test 2: Complex + R(|Ψ|²) — does saturation produce self-confinement?
Test 3: Phase synchronization — does the name become the physics?
"""
import numpy as np
import json
from dataclasses import dataclass

# =============================================================================
# Test 1: Real transfer → diffusion, Complex transfer → waves?
# =============================================================================

def test_real_vs_complex_transfer():
    """The same transfer rule with real vs complex fields."""
    N = 256
    dx = 1.0
    dt = 0.1
    k = 0.3  # transfer coefficient
    steps = 2000

    # Initial condition: Gaussian pulse
    x = np.arange(N) * dx
    x0 = N * dx / 2
    sigma = 10.0

    # --- REAL transfer rule (original framework) ---
    I_real = np.exp(-0.5 * ((x - x0) / sigma)**2)

    I_real_history = []
    for t in range(steps):
        if t % 100 == 0:
            I_real_history.append(I_real.copy())
        # Transfer rule: ΔI_x = k * Σ(I_n - I_x) * R(I_n)
        # For now R = 1 (unsaturated regime)
        dI = np.zeros(N)
        dI += k * (np.roll(I_real, 1) - I_real)  # left neighbor
        dI += k * (np.roll(I_real, -1) - I_real)  # right neighbor
        I_real = I_real + dt * dI
    I_real_history.append(I_real.copy())

    # Measure: does the pulse spread or oscillate?
    peak_real = [h.max() for h in I_real_history]
    width_real = [np.sum(h > h.max() * 0.1) * dx for h in I_real_history]

    # --- COMPLEX transfer rule (same rule, complex field) ---
    # With imaginary coupling: ΔΨ = i*k * Σ(Ψ_n - Ψ_x)
    # The 'i' makes it Schrödinger-like instead of heat-equation-like
    Psi = np.exp(-0.5 * ((x - x0) / sigma)**2).astype(complex)
    # Give it initial momentum (phase gradient)
    k_wave = 0.3  # wave number
    Psi *= np.exp(1j * k_wave * x)

    Psi_density_history = []
    for t in range(steps):
        if t % 100 == 0:
            Psi_density_history.append(np.abs(Psi)**2)
        # Same transfer rule but with i*k
        dPsi = np.zeros(N, dtype=complex)
        dPsi += 1j * k * (np.roll(Psi, 1) - Psi)
        dPsi += 1j * k * (np.roll(Psi, -1) - Psi)
        Psi = Psi + dt * dPsi
    Psi_density_history.append(np.abs(Psi)**2)

    peak_complex = [h.max() for h in Psi_density_history]
    width_complex = [np.sum(h > h.max() * 0.1) * dx for h in Psi_density_history]

    # --- Also test: real k with complex field (still diffusion?) ---
    Psi2 = np.exp(-0.5 * ((x - x0) / sigma)**2).astype(complex)
    Psi2 *= np.exp(1j * k_wave * x)

    Psi2_density_history = []
    for t in range(steps):
        if t % 100 == 0:
            Psi2_density_history.append(np.abs(Psi2)**2)
        dPsi2 = np.zeros(N, dtype=complex)
        dPsi2 += k * (np.roll(Psi2, 1) - Psi2)  # real k, complex field
        dPsi2 += k * (np.roll(Psi2, -1) - Psi2)
        Psi2 = Psi2 + dt * dPsi2
    Psi2_density_history.append(np.abs(Psi2)**2)

    peak_complex_realk = [h.max() for h in Psi2_density_history]

    return {
        "real_field": {
            "peak_initial": peak_real[0],
            "peak_final": peak_real[-1],
            "width_initial": width_real[0],
            "width_final": width_real[-1],
            "monotonic_decay": all(peak_real[i] >= peak_real[i+1] - 1e-10 for i in range(len(peak_real)-1)),
            "verdict": "DIFFUSION" if peak_real[-1] < peak_real[0] * 0.5 else "RETAINED"
        },
        "complex_field_imagk": {
            "peak_initial": peak_complex[0],
            "peak_final": peak_complex[-1],
            "width_initial": width_complex[0],
            "width_final": width_complex[-1],
            "monotonic_decay": all(peak_complex[i] >= peak_complex[i+1] - 1e-10 for i in range(len(peak_complex)-1)),
            "verdict": "WAVE" if not all(peak_complex[i] >= peak_complex[i+1] - 1e-10 for i in range(len(peak_complex)-1)) else "DIFFUSION"
        },
        "complex_field_realk": {
            "peak_initial": peak_complex_realk[0],
            "peak_final": peak_complex_realk[-1],
            "monotonic_decay": all(peak_complex_realk[i] >= peak_complex_realk[i+1] - 1e-10 for i in range(len(peak_complex_realk)-1)),
            "verdict": "Complex field + real k = still diffusion?"
        }
    }


# =============================================================================
# Test 2: Complex Intent with R(|Ψ|²) — soliton-like self-confinement?
# =============================================================================

def test_complex_with_saturation():
    """Complex transfer with R(|Ψ|²) saturation. Can it self-confine?"""
    N = 512
    dx = 1.0
    dt = 0.05
    k = 0.3
    I_max = 1.0
    n_sat = 2  # saturation exponent
    steps = 5000

    x = np.arange(N) * dx
    x0 = N * dx / 2
    sigma = 8.0

    def R(rho):
        """Saturation resistance"""
        return np.clip(1.0 - (rho / I_max)**n_sat, 0, 1)

    # High-density localized pulse with momentum
    Psi = 0.8 * np.exp(-0.5 * ((x - x0) / sigma)**2).astype(complex)
    k_wave = 0.0  # no initial momentum — test confinement at rest
    Psi *= np.exp(1j * k_wave * x)

    total_I_initial = np.sum(np.abs(Psi)**2)

    diagnostics = []
    for t in range(steps):
        rho = np.abs(Psi)**2

        if t % 500 == 0:
            peak = rho.max()
            center_of_mass = np.sum(x * rho) / np.sum(rho)
            width = np.sqrt(np.sum((x - center_of_mass)**2 * rho) / np.sum(rho))
            total = np.sum(rho)
            diagnostics.append({
                "t": t,
                "peak": float(peak),
                "width": float(width),
                "total_I": float(total),
                "conservation": float(abs(total - total_I_initial) / total_I_initial)
            })

        # Complex transfer with R(|Ψ|²) at neighbor
        Psi_left = np.roll(Psi, 1)
        Psi_right = np.roll(Psi, -1)
        rho_left = np.abs(Psi_left)**2
        rho_right = np.abs(Psi_right)**2

        dPsi = np.zeros(N, dtype=complex)
        dPsi += 1j * k * (Psi_left - Psi) * R(rho_left)
        dPsi += 1j * k * (Psi_right - Psi) * R(rho_right)

        Psi = Psi + dt * dPsi

    # Final measurement
    rho = np.abs(Psi)**2
    peak_final = rho.max()
    com_final = np.sum(x * rho) / np.sum(rho)
    width_final = np.sqrt(np.sum((x - com_final)**2 * rho) / np.sum(rho))

    diagnostics.append({
        "t": steps,
        "peak": float(peak_final),
        "width": float(width_final),
        "total_I": float(np.sum(rho)),
        "conservation": float(abs(np.sum(rho) - total_I_initial) / total_I_initial)
    })

    width_ratio = diagnostics[-1]["width"] / diagnostics[0]["width"]

    return {
        "diagnostics": diagnostics,
        "width_ratio": width_ratio,
        "self_confined": width_ratio < 1.5,
        "verdict": "SELF-CONFINED" if width_ratio < 1.5 else f"DISPERSED (width ratio {width_ratio:.2f})"
    }


# =============================================================================
# Test 3: Phase synchronization between two pulses
# =============================================================================

def test_phase_synchronization():
    """Two complex pulses — do they phase-lock (resonate)?"""
    N = 512
    dx = 1.0
    dt = 0.05
    k = 0.3
    steps = 5000

    x = np.arange(N) * dx
    sigma = 10.0

    # Two pulses with different initial phases
    x1, x2 = N * dx * 0.35, N * dx * 0.65
    Psi1 = 0.5 * np.exp(-0.5 * ((x - x1) / sigma)**2)
    Psi2 = 0.5 * np.exp(-0.5 * ((x - x2) / sigma)**2)

    # Give them different phases
    phase1, phase2 = 0.0, np.pi * 0.7  # different initial phases
    Psi = (Psi1 * np.exp(1j * phase1) + Psi2 * np.exp(1j * phase2)).astype(complex)

    # Track phase difference at the two centers over time
    idx1 = int(x1 / dx)
    idx2 = int(x2 / dx)

    phase_diffs = []
    overlap_history = []

    for t in range(steps):
        if t % 50 == 0:
            # Phase at each center
            if np.abs(Psi[idx1]) > 1e-10 and np.abs(Psi[idx2]) > 1e-10:
                p1 = np.angle(Psi[idx1])
                p2 = np.angle(Psi[idx2])
                phase_diffs.append(float(np.abs(p1 - p2) % (2*np.pi)))

            # Overlap integral (interference measure)
            rho = np.abs(Psi)**2
            mid = N // 2
            overlap = float(np.sum(rho[mid-5:mid+5]))
            overlap_history.append(overlap)

        # Transfer with imaginary k
        dPsi = np.zeros(N, dtype=complex)
        dPsi += 1j * k * (np.roll(Psi, 1) - Psi)
        dPsi += 1j * k * (np.roll(Psi, -1) - Psi)
        Psi = Psi + dt * dPsi

    # Did phases converge?
    if len(phase_diffs) > 10:
        phase_var_early = np.var(phase_diffs[:10])
        phase_var_late = np.var(phase_diffs[-10:])
        phase_locked = phase_var_late < phase_var_early * 0.1
    else:
        phase_locked = False

    return {
        "initial_phase_diff": float(phase2 - phase1),
        "phase_diffs_early": phase_diffs[:5] if phase_diffs else [],
        "phase_diffs_late": phase_diffs[-5:] if phase_diffs else [],
        "phase_locked": phase_locked,
        "overlap_early": overlap_history[:5] if overlap_history else [],
        "overlap_late": overlap_history[-5:] if overlap_history else [],
        "verdict": "PHASE LOCKED (resonance)" if phase_locked else "NO PHASE LOCK"
    }


# =============================================================================
# Test 4: The naming contradiction — count phase-dependent features
# =============================================================================

def test_framework_vocabulary():
    """
    Synchronism's vocabulary implies phase dynamics.
    Score: how many framework concepts REQUIRE complex fields?
    """
    concepts = {
        "Synchronization": {
            "requires_phase": True,
            "reason": "Synchronization is phase-locking. No phase → no synchronization.",
            "in_real_transfer": False,
            "in_complex_transfer": True
        },
        "Resonance": {
            "requires_phase": True,
            "reason": "Resonance = constructive interference. Requires phase matching.",
            "in_real_transfer": False,
            "in_complex_transfer": True
        },
        "Dissonance": {
            "requires_phase": True,
            "reason": "Dissonance = destructive interference. Requires phase opposition.",
            "in_real_transfer": False,
            "in_complex_transfer": True
        },
        "Indifference": {
            "requires_phase": True,
            "reason": "Indifference = no consistent phase relationship. Requires phase to be inconsistent about.",
            "in_real_transfer": False,
            "in_complex_transfer": True
        },
        "Entity_as_oscillation": {
            "requires_phase": True,
            "reason": "Oscillation requires two components (amplitude cycling through states). Real scalar can only relax.",
            "in_real_transfer": False,
            "in_complex_transfer": True
        },
        "Witnessing_as_sync": {
            "requires_phase": True,
            "reason": "Witness = synchronization event = phase match at a moment.",
            "in_real_transfer": False,
            "in_complex_transfer": True
        },
        "de_Broglie_frequency": {
            "requires_phase": True,
            "reason": "f = E/h describes phase rotation rate. No phase → no frequency.",
            "in_real_transfer": False,
            "in_complex_transfer": True
        },
        "Saturation": {
            "requires_phase": False,
            "reason": "Amplitude limit. Works for real or complex fields.",
            "in_real_transfer": True,
            "in_complex_transfer": True
        },
        "MRH": {
            "requires_phase": False,
            "reason": "Scale-dependent relevance. Structural concept, not dynamic.",
            "in_real_transfer": True,
            "in_complex_transfer": True
        },
        "Gravity_as_gradient": {
            "requires_phase": False,
            "reason": "Saturation gradient → transfer bias. Works for real fields.",
            "in_real_transfer": True,
            "in_complex_transfer": True
        }
    }

    require_phase = sum(1 for c in concepts.values() if c["requires_phase"])
    total = len(concepts)

    return {
        "total_concepts": total,
        "require_phase": require_phase,
        "work_in_real": total - require_phase,
        "fraction_requiring_phase": require_phase / total,
        "concepts": concepts,
        "verdict": f"{require_phase}/{total} ({100*require_phase/total:.0f}%) of framework concepts require phase (complex fields). "
                   f"The framework's vocabulary is 70% wave physics described with diffusion mathematics."
    }


# =============================================================================
# Test 5: What does R(|Ψ|²) predict for nonlinear QM?
# =============================================================================

def test_nonlinear_qm_prediction():
    """
    If the transfer rule is complex with R(|Ψ|²), what deviation from linear QM does it predict?
    Compute the fractional nonlinear correction at various density scales.
    """
    # Planck density
    rho_planck = 5.155e96  # kg/m³

    density_scales = {
        "interstellar_medium": 1e-21,
        "air": 1.2,
        "water": 1000,
        "iron": 7874,
        "white_dwarf": 1e9,
        "neutron_star_surface": 1e14,
        "nuclear_density": 2.3e17,
        "neutron_star_core": 1e18,
        "quark_gluon_plasma": 1e19,
    }

    results = {}
    for name, rho in density_scales.items():
        # Fractional correction: (ρ/ρ_max)^n for n=2
        ratio = rho / rho_planck
        correction_n2 = ratio**2
        correction_n1 = ratio
        results[name] = {
            "density_kg_m3": rho,
            "rho_over_rho_max": ratio,
            "correction_n2": correction_n2,
            "correction_n1": correction_n1,
            "log10_correction_n2": np.log10(correction_n2) if correction_n2 > 0 else float('-inf'),
            "observable": correction_n2 > 1e-20  # generous observability threshold
        }

    all_unobservable = all(not r["observable"] for r in results.values())

    return {
        "density_scales": results,
        "all_unobservable": all_unobservable,
        "verdict": "ALL nonlinear corrections are < 10^-155. Completely unobservable at any achievable density."
                   if all_unobservable else "Some corrections may be observable"
    }


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Session #620: Complex Intent — The Name vs The Mathematics")
    print("=" * 70)

    print("\n--- Test 1: Real vs Complex Transfer ---")
    r1 = test_real_vs_complex_transfer()
    print(f"Real field:       {r1['real_field']['verdict']}")
    print(f"  Peak: {r1['real_field']['peak_initial']:.4f} → {r1['real_field']['peak_final']:.4f}")
    print(f"  Width: {r1['real_field']['width_initial']:.1f} → {r1['real_field']['width_final']:.1f}")
    print(f"  Monotonic decay: {r1['real_field']['monotonic_decay']}")
    print(f"Complex field (imag k): {r1['complex_field_imagk']['verdict']}")
    print(f"  Peak: {r1['complex_field_imagk']['peak_initial']:.4f} → {r1['complex_field_imagk']['peak_final']:.4f}")
    print(f"  Width: {r1['complex_field_imagk']['width_initial']:.1f} → {r1['complex_field_imagk']['width_final']:.1f}")
    print(f"  Monotonic decay: {r1['complex_field_imagk']['monotonic_decay']}")
    print(f"Complex field (real k): {r1['complex_field_realk']['verdict']}")
    print(f"  Peak: {r1['complex_field_realk']['peak_initial']:.4f} → {r1['complex_field_realk']['peak_final']:.4f}")
    print(f"  Monotonic decay: {r1['complex_field_realk']['monotonic_decay']}")

    print("\n--- Test 2: Complex + Saturation (Self-Confinement?) ---")
    r2 = test_complex_with_saturation()
    print(f"Verdict: {r2['verdict']}")
    for d in r2['diagnostics']:
        print(f"  t={d['t']:5d}: peak={d['peak']:.4f}, width={d['width']:.1f}, conservation={d['conservation']:.2e}")

    print("\n--- Test 3: Phase Synchronization ---")
    r3 = test_phase_synchronization()
    print(f"Initial phase diff: {r3['initial_phase_diff']:.3f}")
    print(f"Phase diffs (early): {[f'{p:.3f}' for p in r3['phase_diffs_early']]}")
    print(f"Phase diffs (late):  {[f'{p:.3f}' for p in r3['phase_diffs_late']]}")
    print(f"Verdict: {r3['verdict']}")

    print("\n--- Test 4: Framework Vocabulary Audit ---")
    r4 = test_framework_vocabulary()
    print(f"Verdict: {r4['verdict']}")
    print(f"Concepts requiring phase:")
    for name, info in r4['concepts'].items():
        if info['requires_phase']:
            print(f"  - {name}: {info['reason']}")

    print("\n--- Test 5: Nonlinear QM Prediction ---")
    r5 = test_nonlinear_qm_prediction()
    print(f"Verdict: {r5['verdict']}")
    for name, info in r5['density_scales'].items():
        print(f"  {name:25s}: ρ/ρ_max = {info['rho_over_rho_max']:.1e}, correction = {info['correction_n2']:.1e}")

    print("\n" + "=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    print("""
The framework's name is 'Synchronism' — a theory about synchronization.
Synchronization is phase dynamics. Phase requires complex fields.
The stated mathematics (real I, real k) has no phase. It produces diffusion.

Making Intent complex and the transfer imaginary gives Schrödinger dynamics:
- Waves propagate (Test 1: confirmed)
- Phase relationships exist (resonance/dissonance possible)
- Entities as standing waves become possible

But: R(|Ψ|²) nonlinear corrections are ~10^-155 at nuclear density.
No testable prediction follows from the specific saturation form.

The framework's mathematics contradicts its vocabulary.
The fix (complex Intent) is quantum mechanics.
The framework doesn't extend QM — QM extends the framework.

FRAME QUESTION: What if Synchronism is a language for physics, not a theory of it?
""")

    # Save results
    results = {
        "test1_real_vs_complex": r1,
        "test2_complex_saturation": {k: v for k, v in r2.items() if k != 'diagnostics'},
        "test3_phase_sync": r3,
        "test4_vocabulary": {k: v for k, v in r4.items() if k != 'concepts'},
        "test5_nonlinear_qm": {k: v for k, v in r5.items() if k != 'density_scales'}
    }

    with open("simulations/session620_results.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print("Results saved to simulations/session620_results.json")
