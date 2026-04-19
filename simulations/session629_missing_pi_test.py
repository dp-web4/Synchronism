"""
Session 629: The Missing π — Does Synchronism Have Universal Constants?

The operator's framing: "Intent is a computational abstraction like π."
This defends Intent against the units critique (it's dimensionless, like π).
But it implies a stronger structural claim: Intent should, like π, bring
specific universal constants to the theory's predictions.

π has these properties:
  - Specific numerical value (3.14159...)
  - Emerges from geometry (not posited)
  - Constrains formulas (can't be replaced with arbitrary number)
  - Universal across all Euclidean contexts

For Synchronism, the natural candidate is k_crit — the critical coupling
at which checkerboard instability (Class 2 oscillation) first emerges.
If k_crit is a universal constant of Intent dynamics, it should be:
  - Independent of dimensionality (1D, 2D, 3D)
  - Independent of n (the exponent in R(I) = [1-ρ^n])
  - Independent of lattice type (regular, random, coordination number)

If k_crit varies with these, it's a parameter, not a constant.
This is the π analogy test.

Method: Measure k_crit across dimensions and n values.
"""
import numpy as np

def find_k_crit(dim, n, L=128, T=400, n_trials=3, k_range=(0.1, 0.7), k_steps=25):
    """Find critical k at which uniform state first becomes unstable.

    Monotonic R(I) = [1-(I/I_max)^n], I_max=1, initial condition near uniform.
    k_crit = smallest k at which perturbation grows (Lyapunov > 0).
    """
    if dim == 1:
        shape = (L,)
    elif dim == 2:
        shape = (L // 4, L // 4)  # keep total volume tractable
    elif dim == 3:
        shape = (16, 16, 16)
    else:
        raise ValueError(dim)

    ks = np.linspace(*k_range, k_steps)
    k_crit_trials = []

    for trial in range(n_trials):
        rng = np.random.default_rng(42 + trial)
        for k in ks:
            # Uniform + tiny perturbation
            I = 0.5 * np.ones(shape) + 1e-4 * rng.standard_normal(shape)
            I0_norm = np.linalg.norm(I - 0.5)

            # Evolve
            for t in range(T):
                R = np.clip(1 - I**n, 0, 1)
                # Sum of neighbor differences weighted by R at neighbor
                laplacian = np.zeros_like(I)
                for ax in range(len(shape)):
                    I_plus = np.roll(I, -1, axis=ax)
                    I_minus = np.roll(I, 1, axis=ax)
                    R_plus = np.roll(R, -1, axis=ax)
                    R_minus = np.roll(R, 1, axis=ax)
                    laplacian += (I_plus - I) * R_plus + (I_minus - I) * R_minus
                I = I + k * laplacian
                I = np.clip(I, 0, 1)

                if not np.isfinite(I).all():
                    break

            I_final_norm = np.linalg.norm(I - np.mean(I))

            # Crit if perturbation grew significantly
            if I_final_norm > 10 * I0_norm:
                k_crit_trials.append(k)
                break
        else:
            k_crit_trials.append(np.nan)

    return np.nanmean(k_crit_trials), np.nanstd(k_crit_trials)


if __name__ == "__main__":
    print("=" * 70)
    print("SESSION 629: Does k_crit behave like π?")
    print("=" * 70)
    print()
    print("If k_crit is a universal Synchronism constant, it should be")
    print("invariant across dimensionality and exponent n.")
    print()
    print(f"{'dim':>4} {'n':>4} {'k_crit':>10} {'±':>10}")
    print("-" * 32)

    results = []
    for dim in [1, 2, 3]:
        for n in [1, 2, 3, 4]:
            k_crit, k_std = find_k_crit(dim, n)
            results.append((dim, n, k_crit, k_std))
            print(f"{dim:>4} {n:>4} {k_crit:>10.4f} {k_std:>10.4f}")

    print()
    print("=" * 70)
    print("ANALYSIS")
    print("=" * 70)

    # Check variation across dim at fixed n
    for n in [1, 2, 3, 4]:
        kc_by_dim = [r[2] for r in results if r[1] == n]
        if all(np.isfinite(kc_by_dim)):
            spread = max(kc_by_dim) - min(kc_by_dim)
            print(f"n={n}: k_crit across dim ranges {min(kc_by_dim):.3f}–{max(kc_by_dim):.3f} (spread {spread:.3f})")

    print()
    for dim in [1, 2, 3]:
        kc_by_n = [r[2] for r in results if r[0] == dim]
        if all(np.isfinite(kc_by_n)):
            spread = max(kc_by_n) - min(kc_by_n)
            print(f"dim={dim}: k_crit across n ranges {min(kc_by_n):.3f}–{max(kc_by_n):.3f} (spread {spread:.3f})")

    print()
    all_kc = [r[2] for r in results if np.isfinite(r[2])]
    print(f"Overall: k_crit range = {min(all_kc):.3f} to {max(all_kc):.3f}")
    print(f"Range/mean = {(max(all_kc)-min(all_kc))/np.mean(all_kc):.1%}")
    print()
    print("VERDICT:")
    if (max(all_kc) - min(all_kc)) / np.mean(all_kc) < 0.05:
        print("  k_crit appears ~universal. Candidate for a Synchronism constant.")
    else:
        print("  k_crit varies with dim/n. It is a parameter, not a constant.")
        print("  The π-analogy fails: no structural invariant emerges.")
