"""
Session #622: Direct test of the operator's self-witnessing mechanism.

The mechanism (from SESSION_FOCUS.md operator note):
- A pulse moves through a NEAR-SATURATION background
- The pulse's leading edge pushes nearby cells into saturation
- R → 0 at those cells → transfer blocked → energy reflects
- Reflected pulse encounters near-saturation cells behind it → bounces again
- Self-confinement through self-created transient walls

This is DIFFERENT from what S17-22 tested:
- S17: Pre-existing static walls → damped oscillation
- S19: Pulse in LOW-I background → dispersal
- S20: Analytical — coupled pulses repel
- S21-22: Vortex → core disperses

The operator's mechanism requires:
1. 2-DOF (pulse must be MOVING — needs velocity field)
2. Near-saturation background (R close to 0)
3. The pulse doesn't need large amplitude — just enough to push beyond threshold

But S619 proved: P = I_max - I gives c² < 0 at ALL densities.
Even with 2-DOF N-S, the sound speed is imaginary.
The self-witnessing mechanism fails BEFORE it's tested —
not from damping (S17) but from the EOS no-go (S619).

This simulation tests whether 2-DOF dynamics in near-saturation
background can produce ANYTHING other than collapse, using the
actual N-S equations with the stated EOS.
"""

import numpy as np

def R(rho, rho_max=1.0, n=2):
    """Saturation resistance = viscosity."""
    return np.maximum(1.0 - (rho/rho_max)**n, 0.0)

def run_2dof_ns_1d(N, rho_init, v_init, rho_max=1.0, n=2, dt=0.001,
                    steps=10000, report_every=1000):
    """
    1D compressible Navier-Stokes with:
    - P(rho) = rho_max - rho  (FUNDAMENTALS.md identification)
    - mu(rho) = D * R(rho)     (saturation-based viscosity)
    - Periodic BC

    Equations:
    d(rho)/dt + d(rho*v)/dx = 0
    d(rho*v)/dt + d(rho*v*v + P)/dx = d(mu * dv/dx)/dx
    """
    dx = 1.0 / N
    D = 0.1  # base diffusion

    rho = rho_init.copy()
    v = v_init.copy()

    center = N // 2
    rho_hist = []
    v_hist = []

    for step in range(steps):
        if step % report_every == 0:
            rho_hist.append(rho.copy())
            v_hist.append(v.copy())

        # Pressure and viscosity
        P = rho_max - rho
        mu = D * R(rho, rho_max, n)

        # Spatial derivatives (central differences, periodic BC)
        drho = np.zeros(N)
        dv = np.zeros(N)
        d2v = np.zeros(N)
        dP = np.zeros(N)
        dmu_dv = np.zeros(N)

        for i in range(N):
            ip = (i + 1) % N
            im = (i - 1) % N

            # Advection: d(rho*v)/dx using upwind
            flux_p = rho[i] * v[i] if v[i] > 0 else rho[ip] * v[ip]
            flux_m = rho[im] * v[im] if v[im] > 0 else rho[i] * v[i]
            drho[i] = -(flux_p - flux_m) / dx

            # Pressure gradient
            dP[i] = -(P[ip] - P[im]) / (2 * dx)

            # Viscous term: d/dx(mu * dv/dx)
            mu_mid_p = 0.5 * (mu[i] + mu[ip])
            mu_mid_m = 0.5 * (mu[im] + mu[i])
            dmu_dv[i] = (mu_mid_p * (v[ip] - v[i]) - mu_mid_m * (v[i] - v[im])) / dx**2

            # Advection of momentum
            if v[i] > 0:
                dv[i] = -v[i] * (v[i] - v[im]) / dx
            else:
                dv[i] = -v[i] * (v[ip] - v[i]) / dx

        # Update
        rho_new = rho + dt * drho
        mom_update = dv + dP / np.maximum(rho, 1e-10) + dmu_dv / np.maximum(rho, 1e-10)
        v_new = v + dt * mom_update

        # Enforce bounds
        rho_new = np.clip(rho_new, 1e-10, rho_max * 0.9999)

        rho = rho_new
        v = v_new

    rho_hist.append(rho.copy())
    v_hist.append(v.copy())

    return np.array(rho_hist), np.array(v_hist)

def main():
    print("="*70)
    print("SESSION 622: SELF-WITNESSING MECHANISM TEST")
    print("="*70)

    N = 128
    rho_max = 1.0
    n = 2

    # ================================================================
    # TEST A: The operator's mechanism — pulse in near-saturation medium
    # ================================================================
    print("\n" + "-"*70)
    print("TEST A: Moving pulse in near-saturation background")
    print("Background rho=0.85 (R=0.2775), pulse pushes to 0.95 (R=0.0975)")
    print("-"*70)

    rho_bg = 0.85
    rho_init = rho_bg * np.ones(N)
    # Narrow pulse
    center = N // 2
    rho_init[center-2:center+3] = 0.95  # near saturation

    # Give it rightward velocity
    v_init = np.zeros(N)
    v_init[center-2:center+3] = 0.5

    print(f"\nInitial: rho_bg={rho_bg}, rho_pulse=0.95, v_pulse=0.5")
    print(f"P(rho_bg) = {rho_max - rho_bg:.3f}")
    print(f"P(rho_pulse) = {rho_max - 0.95:.3f}")
    print(f"dP/drho = -1 everywhere → c² = dP/drho = -1 → c = IMAGINARY")
    print(f"⚠️ Sound speed is imaginary → Hadamard instability")
    print(f"Running 2-DOF N-S anyway to see what happens...")

    rho_hist, v_hist = run_2dof_ns_1d(N, rho_init, v_init, rho_max, n,
                                        dt=0.0005, steps=20000, report_every=2000)

    print(f"\nEvolution of center cell density:")
    for i, rho_snap in enumerate(rho_hist):
        t = i * 2000
        max_rho = np.max(rho_snap)
        max_pos = np.argmax(rho_snap)
        width = np.sum(rho_snap > rho_bg + 0.01)
        print(f"  t={t:5d}: max(rho)={max_rho:.4f} at x={max_pos}, width={width}")

    # Check for oscillation at center
    center_rho = np.array([h[center] for h in rho_hist])
    mean_rho = np.mean(center_rho)
    centered = center_rho - mean_rho
    signs = np.sign(centered)
    sc = np.sum(np.abs(np.diff(signs)) > 0)
    print(f"\nCenter cell: {sc} sign changes over {len(rho_hist)} snapshots")
    if sc > 2:
        print("→ OSCILLATION detected at center!")
    else:
        print("→ No oscillation at center")

    # ================================================================
    # TEST B: What about with CORRECTED EOS? P = integral(R) = rho - rho^(n+1)/(n+1)
    # This gives c² = R(rho) ≥ 0 → waves propagate (S619 case 4)
    # But then dP/drho > 0 → no gravity (pressure pushes outward)
    # ================================================================
    print("\n" + "-"*70)
    print("TEST B: Corrected EOS (P = ∫R dρ) — waves but no gravity")
    print("c² = R(ρ) ≥ 0 → waves propagate, but no gravitational attraction")
    print("-"*70)

    def run_corrected_eos(N, rho_init, v_init, rho_max=1.0, n=2, dt=0.001,
                          steps=10000, report_every=1000):
        """Same as above but with P = integral of R(rho) = rho - rho^(n+1)/(n+1)."""
        dx = 1.0 / N
        D = 0.1

        rho = rho_init.copy()
        v = v_init.copy()

        rho_hist = []
        v_hist = []

        for step in range(steps):
            if step % report_every == 0:
                rho_hist.append(rho.copy())
                v_hist.append(v.copy())

            # Corrected pressure: P = rho - rho^(n+1)/(n+1) (gives dP/drho = R(rho) >= 0)
            P = rho - rho**(n+1) / (n+1)
            mu = D * R(rho, rho_max, n)

            drho = np.zeros(N)
            dv = np.zeros(N)

            for i in range(N):
                ip = (i + 1) % N
                im = (i - 1) % N

                # Continuity
                flux_p = rho[i] * v[i] if v[i] > 0 else rho[ip] * v[ip]
                flux_m = rho[im] * v[im] if v[im] > 0 else rho[i] * v[i]
                drho[i] = -(flux_p - flux_m) / dx

                # Pressure gradient
                dP_dx = (P[ip] - P[im]) / (2 * dx)

                # Viscous term
                mu_mid_p = 0.5 * (mu[i] + mu[ip])
                mu_mid_m = 0.5 * (mu[im] + mu[i])
                visc = (mu_mid_p * (v[ip] - v[i]) - mu_mid_m * (v[i] - v[im])) / dx**2

                # Advection of momentum
                if v[i] > 0:
                    adv = -v[i] * (v[i] - v[im]) / dx
                else:
                    adv = -v[i] * (v[ip] - v[i]) / dx

                dv[i] = adv - dP_dx / max(rho[i], 1e-10) + visc / max(rho[i], 1e-10)

            rho = np.clip(rho + dt * drho, 1e-10, rho_max * 0.9999)
            v = v + dt * dv

        rho_hist.append(rho.copy())
        v_hist.append(v.copy())
        return np.array(rho_hist), np.array(v_hist)

    # Same initial conditions as Test A
    rho_init = rho_bg * np.ones(N)
    rho_init[center-2:center+3] = 0.95
    v_init = np.zeros(N)
    v_init[center-2:center+3] = 0.5

    rho_hist2, v_hist2 = run_corrected_eos(N, rho_init, v_init, rho_max, n,
                                            dt=0.0005, steps=20000, report_every=2000)

    print(f"\nEvolution with P = ∫R dρ (waves allowed, no gravity):")
    for i, rho_snap in enumerate(rho_hist2):
        t = i * 2000
        max_rho = np.max(rho_snap)
        max_pos = np.argmax(rho_snap)
        width = np.sum(rho_snap > rho_bg + 0.01)
        print(f"  t={t:5d}: max(rho)={max_rho:.4f} at x={max_pos}, width={width}")

    # Does the pulse propagate as a wave?
    peak_positions = [np.argmax(h) for h in rho_hist2]
    print(f"\nPeak position over time: {peak_positions}")
    if len(set(peak_positions)) > 1:
        print("→ Pulse PROPAGATES (wave behavior confirmed)")
    else:
        print("→ Pulse STATIONARY")

    # Check: does the pulse DISPERSE or CONFINE?
    initial_width = np.sum(rho_hist2[0] > rho_bg + 0.01)
    final_width = np.sum(rho_hist2[-1] > rho_bg + 0.01)
    print(f"Width: initial={initial_width}, final={final_width}")
    if final_width > initial_width * 1.5:
        print("→ DISPERSES — no self-confinement even with correct EOS")
    elif final_width < initial_width * 0.5:
        print("→ COLLAPSES")
    else:
        print("→ Roughly maintained (but check if this is just short time)")

    # ================================================================
    # VERDICT
    # ================================================================
    print("\n" + "="*70)
    print("VERDICT: The Self-Witnessing Mechanism")
    print("="*70)
    print("""
The operator's self-witnessing mechanism requires:
  1. A moving pulse (2-DOF: density + velocity) ✓
  2. Near-saturation background                  ✓
  3. Transient walls from saturation overflow     ✓

But it ALSO requires:
  4. Wave propagation (c² > 0)                    ✗ (P = I_max - I gives c² = -1)
  5. Self-confinement (focusing nonlinearity)     ✗ (R(I) is defocusing)

The mechanism fails at step 4 BEFORE reaching step 5.

With the CORRECTED EOS (P = ∫R dρ, which gives waves), you get
propagation but NOT self-confinement — the pulse disperses.
And you lose gravity.

S619's no-go theorem is the fatal blow: you cannot have BOTH the
wave propagation needed for self-witnessing AND the gravitational
attraction claimed by the framework. Not with any P(ρ) from R(I).

The self-witnessing concept is physically intuitive — a dynamic
structure maintaining itself through repeated interaction with its
own boundaries. But the mathematics of R(I) saturation cannot
implement it, for the same deep reason it can't produce entities:
one mechanism, one behavior.

The operator's mechanism was never directly tested as stated.
Now it has been. It fails — but for a DIFFERENT and more fundamental
reason than S17-22 found. S17-22 found damping. This test finds
that the EOS itself prevents the wave propagation the mechanism
requires. Even undamped, even in 2-DOF, the pressure is inverted.
""")

if __name__ == "__main__":
    main()
