#!/usr/bin/env python3
"""
Phase-14 (gravity, RADIATIVE sector) — the inflow model has no spin-2 gravitational waves,
because the substrate is flat-space + a scalar (irrotational) flow. This BREAKS the
"GR-in-coordinates" verdict in the one sector the arc never tested: radiation.

Context
-------
Phase-3c/9 established the inflow/swimmer gravity reproduces GR EXACTLY — but only for
STATIC, spherical configurations (Schwarzschild: light bending, precession). The verdict
"Synchronism gravity = Gullstrand-Painlevé = GR in coordinates" was getting tidy. A
coordinate change reproduces GR by definition; the place to stress it is a
diffeomorphism-INVARIANT observable the static tests don't probe: GRAVITATIONAL RADIATION.

The framework's defining commitments (confirmed in PREDICTIONS.md / FUNDAMENTALS.md):
  - space is a FLAT absolute Planck grid (Phase-3c: "absolute-time, flat-space substrate");
  - the substrate is a 1-DOF SCALAR, IRROTATIONAL (curl-free) field (S617/S665/S666).
The eikonal inflow defines the acoustic/Gullstrand-Painlevé metric
    ds² = -(c²-u²)dt² - 2 u·dx dt + δ_ij dx^i dx^j,
whose SPATIAL part is FLAT: g_ij = δ_ij. All "gravity" sits in g_00 (scalar/Newtonian) and
g_0i (vector/gravitomagnetic). GR's gravitational waves are the TRANSVERSE-TRACELESS
(spin-2) part of the SPATIAL metric, h_ij^TT. A permanently-flat spatial metric, or any
metric built from a single scalar field, has h_ij^TT ≡ 0.

This script demonstrates:
  (A) The transverse-traceless (spin-2) content of ANY scalar-sourced spatial metric
      (h_ij = δ_ij ψ + ∂_i∂_j χ) is identically zero — verified numerically on a 3D
      Fourier grid — while a genuine tensor source has TT content ≈ 1. So a scalar/
      irrotational substrate carries NO spin-2 radiative degrees of freedom.
  (B) The acoustic/GP metric the arc uses has g_ij = δ_ij (flat spatial) ⇒ h_ij^TT = 0,
      independently of (A): no graviton sector at all.
  (C) Consequence: the inflow model predicts ZERO quadrupole (spin-2) gravitational-wave
      luminosity. Binary-pulsar timing (Hulse-Taylor PSR B1913+16: observed orbital decay
      matches the GR quadrupole formula to 0.2%) and LIGO/Virgo spin-2 polarization tests
      then REFUTE it: the framework's flat-space + scalar-flow gravity cannot produce the
      observed tensor radiation. This is a structural refutation tied to a CORE commitment,
      not a tunable parameter — the analog-gravity limitation (acoustic metrics reproduce
      kinematics, not the Einstein/radiative sector) made concrete for Synchronism.

numpy only. Headless. Writes results JSON.
Author: CBP-Claude (Opus 4.8), autonomous.
"""
import json, os
import numpy as np


def tt_fraction(h_ij_k, kx, ky, kz):
    """Given h_ij in Fourier space (arrays over a 3D k-grid, 6 independent components),
       return ||h^TT||² / ||h||² using the transverse-traceless projector
       Λ_ij,kl = P_ik P_jl − ½ P_ij P_kl,  P_ij = δ_ij − k_i k_j/|k|²."""
    shape = kx.shape
    k2 = kx**2 + ky**2 + kz**2
    k2[k2 == 0] = np.inf   # kill k=0 (no projector / no radiation)
    K = [kx, ky, kz]
    # projector P_ij
    P = [[ (1.0 if i==j else 0.0) - K[i]*K[j]/k2 for j in range(3)] for i in range(3)]
    # h as 3x3 symmetric
    h = h_ij_k
    # TT projection: hTT_ij = (P_ik P_jl - 0.5 P_ij P_kl) h_kl
    hTT = [[np.zeros(shape, dtype=complex) for _ in range(3)] for _ in range(3)]
    # first term P_ik P_jl h_kl
    for i in range(3):
        for j in range(3):
            s = np.zeros(shape, dtype=complex)
            for k in range(3):
                for l in range(3):
                    s += P[i][k]*P[j][l]*h[k][l]
            hTT[i][j] = s
    # trace term 0.5 P_ij (P_kl h_kl)
    Ph = np.zeros(shape, dtype=complex)
    for k in range(3):
        for l in range(3):
            Ph += P[k][l]*h[k][l]
    for i in range(3):
        for j in range(3):
            hTT[i][j] -= 0.5*P[i][j]*Ph
    # norms (sum over grid and components)
    def norm2(t):
        return sum(np.sum(np.abs(t[i][j])**2) for i in range(3) for j in range(3))
    return float(norm2(hTT)/norm2(h))


def main():
    n = 24
    L = 2*np.pi
    k1 = np.fft.fftfreq(n, d=L/n)*2*np.pi
    kx, ky, kz = np.meshgrid(k1, k1, k1, indexing='ij')
    rng = np.random.default_rng(0)

    print("="*76)
    print("Phase-14: does the flat-space + scalar-flow substrate carry spin-2 grav. waves?")
    print("="*76)

    # a generic scalar field in Fourier space
    phi = (rng.standard_normal((n,n,n)) + 1j*rng.standard_normal((n,n,n)))
    phi[0,0,0] = 0

    print("\n(A) Transverse-traceless (spin-2) content of candidate spatial metrics h_ij:")

    # scalar-trace metric  h_ij = δ_ij φ
    h_trace = [[ (phi if i==j else np.zeros_like(phi)) for j in range(3)] for i in range(3)]
    f_trace = tt_fraction(h_trace, kx, ky, kz)
    print(f"    scalar-trace      h_ij = δ_ij·φ          : ||h^TT||²/||h||² = {f_trace:.3e}")

    # scalar-longitudinal metric h_ij = ∂_i∂_j φ  (Fourier: -k_i k_j φ)
    K=[kx,ky,kz]
    h_long = [[ -K[i]*K[j]*phi for j in range(3)] for i in range(3)]
    f_long = tt_fraction(h_long, kx, ky, kz)
    print(f"    scalar-longitudinal h_ij = ∂_i∂_j φ       : ||h^TT||²/||h||² = {f_long:.3e}")

    # general scalar metric a·δ_ij φ + b·∂_i∂_j χ
    chi = (rng.standard_normal((n,n,n)) + 1j*rng.standard_normal((n,n,n))); chi[0,0,0]=0
    h_gen = [[ 0.7*(phi if i==j else 0.0) - 1.3*K[i]*K[j]*chi for j in range(3)] for i in range(3)]
    f_gen = tt_fraction(h_gen, kx, ky, kz)
    print(f"    general scalar    aδ_ijφ + b∂_i∂_jχ       : ||h^TT||²/||h||² = {f_gen:.3e}")

    # control: a GENUINE tensor source (random symmetric tensor) keeps O(1) TT content
    H = [[ (rng.standard_normal((n,n,n))+1j*rng.standard_normal((n,n,n))) for j in range(3)] for i in range(3)]
    for i in range(3):
        for j in range(3):
            H[i][j] = 0.5*(H[i][j]+H[j][i])  # symmetric
            H[i][j][0,0,0]=0
    f_tensor = tt_fraction(H, kx, ky, kz)
    print(f"    control: generic tensor source            : ||h^TT||²/||h||² = {f_tensor:.3e}")
    print("    => any SCALAR-sourced spatial metric has ZERO spin-2 (TT) content; only a genuine")
    print("       tensor source radiates spin-2. The framework's 1-DOF scalar/irrotational")
    print("       substrate has NO graviton sector.")

    print("\n(B) The acoustic / Gullstrand-Painlevé metric the arc uses:")
    print("    ds² = -(c²-u²)dt² - 2 u·dx dt + δ_ij dx^i dx^j")
    print("    spatial part g_ij = δ_ij (FLAT) ⇒ h_ij = 0 ⇒ h_ij^TT = 0, trivially. Gravity lives")
    print("    only in g_00 (scalar/Newtonian) and g_0i (vector/gravitomagnetic). No tensor sector.")

    print("\n(C) Radiative consequence and refutation:")
    # GR quadrupole luminosity for a fiducial binary (order-of-magnitude), vs inflow = 0.
    G=6.674e-11; c=2.998e8; Msun=1.989e30
    m1=m2=1.4*Msun; a=1.95e9; e=0.617   # ~Hulse-Taylor-like
    mu=m1*m2/(m1+m2); M=m1+m2
    # Peters (1964) orbit-averaged GW luminosity
    f_e = (1 + (73/24)*e**2 + (37/96)*e**4)/(1-e**2)**3.5
    L_GR = (32.0/5.0)*(G**4/c**5)*(m1*m2)**2*M/a**5 * f_e
    print(f"    fiducial Hulse-Taylor-like binary: GR quadrupole (spin-2) luminosity L_GR ≈ {L_GR:.3e} W")
    print(f"    inflow / flat-space-scalar model: TENSOR (spin-2) luminosity = 0 — no TT sector (A,B).")
    print(f"    It is NOT non-radiating: the scalar substrate radiates at QUADRUPOLE order, but with")
    print(f"    SCALAR (spin-0, 'breathing') polarization and a different numerical coefficient — plus")
    print(f"    generic DIPOLE radiation if a star's intent-charge ≠ its mass (strongly excluded).")
    print(f"    The refutation is therefore WRONG POLARIZATION, not 'no decay':")
    print(f"     • OBSERVED (PSR B1913+16): orbital-decay / GR-SPIN-2-quadrupole = 1.0013 ± 0.0021")
    print(f"       (Weisberg & Huang 2016): the observed rate matches the TENSOR coefficient to ~0.2%.")
    print(f"       Scalar quadrupole has a different coefficient ⇒ generically misses this match;")
    print(f"       any dipole channel misses it grossly.")
    print(f"     • LIGO/Virgo multi-detector polarization tests (e.g. GW170817) favor PURE TENSOR")
    print(f"       (+,×) and disfavor scalar/breathing GW modes.")
    print(f"    => flat-space + scalar-flow gravity predicts the WRONG GW polarization (spin-0, not")
    print(f"       spin-2); refuted by polarization tests and the precise binary-pulsar tensor match.")

    print("\nConclusion")
    print("-"*76)
    print("'Synchronism gravity = GR in coordinates' holds ONLY in the static sector (Phase-9).")
    print("In the RADIATIVE sector it FAILS structurally: flat absolute space + a scalar")
    print("irrotational flow ⇒ no transverse-traceless spatial curvature ⇒ no spin-2 gravitational")
    print("waves ⇒ refuted by binary-pulsar timing and direct GW polarization. The break is tied to")
    print("the framework's CORE commitments (flat Planck grid + scalar substrate), not a parameter —")
    print("the same pattern as Phase-13 (absolute time → dim-4 LIV): the defining ontology forces")
    print("an already-refuted prediction in the dynamical/radiative sector, while only reproducing")
    print("known physics (as a coordinate change) in the static sector. Bucket 0 unchanged (0); this")
    print("is a refutation candidate (Bucket 2) for the inflow gravity's radiative completion.")

    out=dict(
        tt_fraction_scalar_trace=f_trace, tt_fraction_scalar_long=f_long,
        tt_fraction_general_scalar=f_gen, tt_fraction_tensor_control=f_tensor,
        L_GR_fiducial_W=L_GR, L_inflow_tensor_W=0.0,
        hulse_taylor_obs_over_GR="1.0013 ± 0.0021 (Weisberg & Huang 2016)",
        radiation_note="scalar substrate radiates at quadrupole order but with SCALAR (spin-0, "
                   "breathing) polarization, not GR's tensor (spin-2); refutation is wrong "
                   "polarization + wrong coefficient, not 'no decay'.",
        conclusion="flat-space + scalar irrotational substrate has zero spin-2 (TT) radiative DOF; "
                   "acoustic/GP metric has flat spatial sections (g_ij=delta_ij); GW polarization is "
                   "scalar not tensor; refuted by binary-pulsar tensor-quadrupole match (0.2%) + "
                   "LIGO spin-2 polarization tests. 'GR in coordinates' holds only in the static "
                   "sector; radiative sector is structurally non-GR and refuted, tied to core "
                   "commitments (flat Planck grid + scalar substrate).")
    os.makedirs("simulations/results",exist_ok=True)
    with open("simulations/results/phase14_no_tensor_gravitational_waves_result.json","w") as f:
        json.dump(out,f,indent=2)
    print("\nsaved -> simulations/results/phase14_no_tensor_gravitational_waves_result.json")


if __name__=="__main__":
    main()
