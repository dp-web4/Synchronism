"""
Session 618: Waveguide Hypothesis Test

QUESTION: Does density-dependent viscosity μ(ρ) = D·[1-(ρ/ρ_max)^n] in 2-DOF
compressible Navier-Stokes create a natural waveguide?

The idea: If μ → 0 at high ρ, then dense regions are nearly inviscid (superfluid-like)
while sparse regions are highly viscous. A wave launched inside a dense region should
stay confined because the surrounding medium damps it.

This would mean the SPECIFIC form of R(I) produces self-confining structures in 2-DOF,
which generic constant-viscosity N-S does not.

SETUP:
- 1D compressible N-S with density-dependent viscosity
- Dense core (ρ ≈ 0.9 ρ_max) embedded in sparse background (ρ ≈ 0.1 ρ_max)
- Small velocity perturbation in the center
- Compare: R(I) viscosity vs constant viscosity (matched mean)
- Measure: how much wave energy remains in the dense core over time

EQUATIONS (1D compressible):
  ∂ρ/∂t + ∂(ρv)/∂x = 0
  ∂(ρv)/∂t + ∂(ρv² + P)/∂x = ∂/∂x[μ(ρ)·∂v/∂x]

  P = K·ρ^γ  (polytropic EOS, γ=1 isothermal for simplicity)
  μ(ρ) = D·[1-(ρ/ρ_max)^n]  (Synchronism's R(I))

CONTROLS:
  A) μ(ρ) = D·[1-(ρ/ρ_max)^n]  — Synchronism viscosity
  B) μ = D_eff (constant, matched to mean of A)  — standard N-S
  C) μ = 0 (inviscid) — baseline
"""

import numpy as np
import json
from pathlib import Path

# ── Parameters ──
Nx = 1000       # grid points
Lx = 10.0       # domain length
dx = Lx / Nx
x = np.linspace(0.5*dx, Lx - 0.5*dx, Nx)

rho_max = 1.0
n_sat = 2       # saturation exponent
D_visc = 0.1    # base viscosity coefficient
K_eos = 1.0     # EOS coefficient (P = K*rho for gamma=1)
gamma_eos = 1.0 # isothermal

# Time stepping — CFL = c_s * dt/dx < 0.5; c_s = sqrt(K) = 1
dt = 0.4 * dx  # CFL = 0.4, safe for Lax-Friedrichs
T_final = 3.0
N_steps = int(T_final / dt)
save_every = max(1, N_steps // 200)  # save ~200 snapshots

# ── Initial conditions ──
# Dense core from x=3 to x=7, sparse background
rho_core = 0.9 * rho_max
rho_bg = 0.1 * rho_max
core_mask = (x > 3.0) & (x < 7.0)
transition_width = 0.3  # smooth transition

def init_density():
    """Smooth density profile: high in core, low outside."""
    rho = np.full(Nx, rho_bg)
    # Smooth tanh transitions
    rho += (rho_core - rho_bg) * 0.5 * (
        np.tanh((x - 3.0) / transition_width) - np.tanh((x - 7.0) / transition_width)
    )
    return rho

def init_velocity():
    """Small Gaussian velocity perturbation in the center."""
    v = 0.01 * np.exp(-((x - 5.0)**2) / 0.1)
    return v

# ── Viscosity models ──
def mu_synchronism(rho):
    """Density-dependent: μ = D·[1 - (ρ/ρ_max)^n]"""
    return D_visc * (1.0 - np.clip(rho / rho_max, 0, 1)**n_sat)

def mu_constant(rho):
    """Constant viscosity matched to spatial mean of Synchronism model."""
    rho0 = init_density()
    mu_mean = np.mean(mu_synchronism(rho0))
    return np.full_like(rho, mu_mean)

def mu_inviscid(rho):
    """Zero viscosity."""
    return np.zeros_like(rho)

# ── Solver: Lax-Friedrichs + central diffusion ──
def pressure(rho):
    return K_eos * rho**gamma_eos

def sound_speed(rho):
    return np.sqrt(gamma_eos * K_eos * np.maximum(rho, 1e-10)**(gamma_eos - 1))

def step(rho, mom, mu_func, dt_step):
    """One timestep of 1D compressible N-S with density-dependent viscosity.

    Conservative variables: rho, mom = rho*v
    """
    v = mom / np.maximum(rho, 1e-10)
    P = pressure(rho)
    mu = mu_func(rho)

    # Fluxes
    flux_rho = mom  # rho * v
    flux_mom = mom * v + P  # rho*v^2 + P

    # Lax-Friedrichs flux splitting
    # Max wave speed for stability
    c = sound_speed(rho)
    alpha = np.max(np.abs(v) + c) * 1.1  # slight over-estimate for safety

    # Left/right states (periodic boundaries)
    rho_L = rho
    rho_R = np.roll(rho, -1)
    mom_L = mom
    mom_R = np.roll(mom, -1)

    fL_rho = flux_rho
    fR_rho = np.roll(flux_rho, -1)
    fL_mom = flux_mom
    fR_mom = np.roll(flux_mom, -1)

    # Numerical flux at i+1/2
    nf_rho = 0.5 * (fL_rho + fR_rho) - 0.5 * alpha * (rho_R - rho_L)
    nf_mom = 0.5 * (fL_mom + fR_mom) - 0.5 * alpha * (mom_R - mom_L)

    # Flux differences: F_{i+1/2} - F_{i-1/2}
    drho = (nf_rho - np.roll(nf_rho, 1)) / dx
    dmom = (nf_mom - np.roll(nf_mom, 1)) / dx

    # Viscous term: ∂/∂x[μ(ρ)·∂v/∂x]
    # Central differences
    dv_dx = (np.roll(v, -1) - np.roll(v, 1)) / (2 * dx)
    mu_dv = mu * dv_dx
    visc_term = (np.roll(mu_dv, -1) - np.roll(mu_dv, 1)) / (2 * dx)

    # Update
    rho_new = rho - dt_step * drho
    mom_new = mom - dt_step * dmom + dt_step * visc_term

    # Floor
    rho_new = np.maximum(rho_new, 1e-10)

    return rho_new, mom_new

def measure_core_energy(rho, mom, core_mask):
    """Kinetic energy in the core region."""
    v = mom / np.maximum(rho, 1e-10)
    KE = 0.5 * rho * v**2
    return np.sum(KE[core_mask]) * dx

def measure_total_energy(rho, mom):
    """Total kinetic + internal energy."""
    v = mom / np.maximum(rho, 1e-10)
    KE = 0.5 * rho * v**2
    IE = K_eos * rho**gamma_eos / (gamma_eos - 1) if gamma_eos != 1 else K_eos * rho * np.log(np.maximum(rho, 1e-10))
    return np.sum(KE + IE) * dx

def run_simulation(mu_func, label):
    """Run one simulation and return diagnostics."""
    rho = init_density()
    v = init_velocity()
    mom = rho * v

    core = (x > 3.0) & (x < 7.0)

    results = {
        'label': label,
        'times': [],
        'core_KE': [],
        'total_KE': [],
        'total_mass': [],
        'max_v': [],
        'v_snapshots': [],
        'rho_snapshots': [],
        'snapshot_times': [],
    }

    KE0_core = measure_core_energy(rho, mom, core)
    KE0_total = np.sum(0.5 * rho * v**2) * dx

    print(f"\n{'='*60}")
    print(f"Running: {label}")
    print(f"  Initial core KE: {KE0_core:.6e}")
    print(f"  Initial total KE: {KE0_total:.6e}")
    print(f"  Steps: {N_steps}, dt={dt:.2e}, dx={dx:.4f}")

    for step_i in range(N_steps + 1):
        t = step_i * dt

        if step_i % save_every == 0:
            v_curr = mom / np.maximum(rho, 1e-10)
            core_KE = measure_core_energy(rho, mom, core)
            total_KE = np.sum(0.5 * rho * v_curr**2) * dx

            results['times'].append(t)
            results['core_KE'].append(core_KE)
            results['total_KE'].append(total_KE)
            results['total_mass'].append(np.sum(rho) * dx)
            results['max_v'].append(np.max(np.abs(v_curr)))

        # Save a few full snapshots
        if step_i in [0, N_steps//10, N_steps//4, N_steps//2, N_steps]:
            v_curr = mom / np.maximum(rho, 1e-10)
            results['v_snapshots'].append(v_curr.tolist())
            results['rho_snapshots'].append(rho.tolist())
            results['snapshot_times'].append(t)

        if step_i < N_steps:
            rho, mom = step(rho, mom, mu_func, dt)

            # Check for blowup
            if np.any(np.isnan(rho)) or np.any(np.isnan(mom)) or np.max(np.abs(mom)) > 1e10:
                print(f"  BLOWUP at step {step_i}, t={t:.4f}")
                break

    # Final diagnostics
    v_final = mom / np.maximum(rho, 1e-10)
    core_KE_final = measure_core_energy(rho, mom, core)
    total_KE_final = np.sum(0.5 * rho * v_final**2) * dx

    print(f"  Final core KE: {core_KE_final:.6e}")
    print(f"  Final total KE: {total_KE_final:.6e}")
    if KE0_core > 0:
        print(f"  Core KE retention: {core_KE_final/KE0_core*100:.2f}%")
    if KE0_total > 0:
        print(f"  Total KE retention: {total_KE_final/KE0_total*100:.2f}%")
    print(f"  Mass conservation: {results['total_mass'][-1]/results['total_mass'][0]*100:.6f}%")

    return results

# ── Run all three ──
if __name__ == '__main__':
    print("Session 618: Waveguide Hypothesis Test")
    print(f"Grid: {Nx} points, domain: [0, {Lx}]")
    print(f"Core: ρ={rho_core}, Background: ρ={rho_bg}")
    print(f"Viscosity: D={D_visc}, n={n_sat}")
    print(f"EOS: P = {K_eos}·ρ^{gamma_eos}")

    # Compute and report mean viscosity for matching
    rho0 = init_density()
    mu_sync = mu_synchronism(rho0)
    print(f"\nViscosity profile:")
    print(f"  Core (ρ=0.9): μ = {D_visc * (1 - 0.9**n_sat):.4f}")
    print(f"  Background (ρ=0.1): μ = {D_visc * (1 - 0.1**n_sat):.4f}")
    print(f"  Spatial mean: μ = {np.mean(mu_sync):.4f}")
    print(f"  Ratio (bg/core): {(1 - 0.1**n_sat)/(1 - 0.9**n_sat):.2f}x")

    results_sync = run_simulation(mu_synchronism, "Synchronism μ(ρ)")
    results_const = run_simulation(mu_constant, "Constant μ (matched mean)")
    results_inv = run_simulation(mu_inviscid, "Inviscid (μ=0)")

    # ── Compare ──
    print("\n" + "="*60)
    print("COMPARISON: Core KE Retention Over Time")
    print("="*60)

    # Compare at several time points
    for frac in [0.1, 0.25, 0.5, 1.0]:
        idx = min(int(frac * len(results_sync['times'])), len(results_sync['times'])-1)
        t = results_sync['times'][idx]

        ke_s = results_sync['core_KE'][idx]
        ke_c = results_const['core_KE'][idx]
        ke_i = results_inv['core_KE'][idx]

        ke0 = results_sync['core_KE'][0]

        print(f"\n  t = {t:.2f}:")
        if ke0 > 0:
            print(f"    Synchronism:  {ke_s/ke0*100:.2f}%")
            print(f"    Constant:     {ke_c/ke0*100:.2f}%")
            print(f"    Inviscid:     {ke_i/ke0*100:.2f}%")

    # ── The key test ──
    print("\n" + "="*60)
    print("KEY TEST: Does Synchronism viscosity confine waves better?")
    print("="*60)

    # Final core KE ratios
    ke0 = results_sync['core_KE'][0]
    if ke0 > 0 and len(results_sync['core_KE']) > 1:
        sync_final = results_sync['core_KE'][-1] / ke0
        const_final = results_const['core_KE'][-1] / ke0
        inv_final = results_inv['core_KE'][-1] / ke0

        print(f"\n  Final core KE retention:")
        print(f"    Synchronism μ(ρ): {sync_final*100:.2f}%")
        print(f"    Constant μ:       {const_final*100:.2f}%")
        print(f"    Inviscid:         {inv_final*100:.2f}%")

        if sync_final > const_final * 1.1:
            print(f"\n  ✓ WAVEGUIDE EFFECT DETECTED")
            print(f"    Synchronism retains {sync_final/const_final:.2f}x more core energy")
            print(f"    This is a SPECIFIC prediction of density-dependent viscosity")
        elif sync_final < const_final * 0.9:
            print(f"\n  ✗ ANTI-WAVEGUIDE: Synchronism loses MORE energy from core")
            print(f"    Constant viscosity retains {const_final/sync_final:.2f}x more")
        else:
            print(f"\n  ≈ NO SIGNIFICANT DIFFERENCE between viscosity models")
            print(f"    Ratio: {sync_final/const_final:.3f}")

    # ── Also test: does the density profile itself evolve differently? ──
    print("\n" + "="*60)
    print("DENSITY EVOLUTION: Does the core structure persist?")
    print("="*60)

    for label, res in [("Synchronism", results_sync), ("Constant", results_const), ("Inviscid", results_inv)]:
        if len(res['rho_snapshots']) >= 2:
            rho_init = np.array(res['rho_snapshots'][0])
            rho_final = np.array(res['rho_snapshots'][-1])
            core_init = np.mean(rho_init[(np.array(x) > 4) & (np.array(x) < 6)])
            core_final = np.mean(rho_final[(np.array(x) > 4) & (np.array(x) < 6)])
            print(f"  {label}: core ρ {core_init:.4f} → {core_final:.4f} ({core_final/core_init*100:.1f}%)")

    # ── Supplementary: test with stronger contrast ──
    print("\n" + "="*60)
    print("SUPPLEMENTARY: Higher saturation (ρ_core = 0.99)")
    print("="*60)

    # Save original and modify
    rho_core_orig = rho_core

    # Redefine init for higher contrast
    def init_density_high():
        rho = np.full(Nx, 0.05 * rho_max)
        rho += (0.99 * rho_max - 0.05 * rho_max) * 0.5 * (
            np.tanh((x - 3.0) / transition_width) - np.tanh((x - 7.0) / transition_width)
        )
        return rho

    # Override init_density temporarily
    import types
    original_init = init_density

    # Monkey-patch for high-contrast test
    globals()['init_density'] = init_density_high

    results_sync_hi = run_simulation(mu_synchronism, "Synchronism μ(ρ), ρ_core=0.99")
    results_const_hi = run_simulation(mu_constant, "Constant μ, ρ_core=0.99")

    ke0_hi = results_sync_hi['core_KE'][0]
    if ke0_hi > 0:
        sync_hi = results_sync_hi['core_KE'][-1] / ke0_hi
        const_hi = results_const_hi['core_KE'][-1] / ke0_hi
        print(f"\n  High-contrast final core KE retention:")
        print(f"    Synchronism: {sync_hi*100:.2f}%")
        print(f"    Constant:    {const_hi*100:.2f}%")
        if sync_hi > 0 and const_hi > 0:
            print(f"    Ratio: {sync_hi/const_hi:.3f}")

    # Restore
    globals()['init_density'] = original_init

    print("\n" + "="*60)
    print("SESSION 618 COMPLETE")
    print("="*60)
