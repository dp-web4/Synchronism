#!/usr/bin/env python3
"""
Phase-12 (door #2, discreteness/LIV) — stress-test the "symmetry-protected untestable"
conclusion: close the named leak, then find the channel the leak-argument omits.

Context
-------
Phase-9 reduced the Bucket-0 frontier to two doors; Phase-11 closed door #1 (galactic =
MOND cage). That left door #2 (discreteness/LIV) as the sole channel. A same-day proposal
(`Research/proposals/liv_frontier_symmetry_protected_untestable.md`) argues door #2 is ALSO
untestable via a three-lock argument on the TIME-OF-FLIGHT channel:
  (1) even-k lattice dispersion ω²=m²+2c²(1−cos ka)/a² ⇒ leading LIV is QUADRATIC (n=2),
      no linear (n=1) term. n=1 would already be REFUTED (GRB bounds > E_Planck).
  (2) n=2 prediction ≈ E_Pl sits ~10⁷ below the n=2 bound ⇒ unreachable.
  (3) n=2 LIV is generic (LQG/causal-sets/DSR) ⇒ non-confirming even if reached.
It flags ONE leak worth a session: the substrate has a PREFERRED FRAME (universal clock);
does that produce a reachable sidereal-modulation or effective-n=1 signal?

This script:
  (A) CLOSES THE LEAK. Boost the even-k lattice dispersion to a lab frame moving at β
      relative to the substrate. Show the lab-frame LIV correction stays ∝ a² (n=2) — a
      Lorentz boost cannot change the power of the lattice spacing a — so the preferred
      frame gives an ANISOTROPIC n=2 modulation (∝ β·(E/E_Pl)²), NOT an effective n=1.
      Strictly weaker than the isotropic n=2 already deemed unreachable. Leak sealed.
  (B) Tree-level pattern universality: on one substrate with one (pattern-unaware) cutoff,
      all low-energy excitations share the same c; differences are O((ka)²) = n=2. So the
      tree-level dim-4 "species-dependent c" channel is also Planck-suppressed (safe).
  (C) The channel the leak-argument OMITS. The even-k symmetry forbids ODD-in-k (n=1,
      dimension-5) terms — but a dimension-4 LIV operator (different limiting speed,
      ω²=m²+(1+c_LIV)k²) is ALSO EVEN IN k, so the symmetry does NOT forbid it. Dim-4 LIV
      is bounded at ~1e-22 (Hughes-Drever/cavity/clock), NOT Planck-suppressed. Tree-level
      universality (B) keeps it ~(m/E_Pl)² — but RADIATIVE corrections in an interacting
      discrete theory generically percolate Planck-scale LIV UP into dim-4 at O(α/π)
      unless a symmetry protects (Collins-Perez-Sudarsky-Gambini-Pullin, PRL 93, 191301,
      2004). The even-k/cubic symmetry does NOT protect dim-4 species-LIV. So door #2 has a
      REACHABLE, refutation-capable sub-channel the three-lock omits.

numpy/scipy. Headless. Writes results JSON.
Author: CBP-Claude (Opus 4.8), autonomous.
"""
import json, os
import numpy as np
from scipy.optimize import brentq


# ----------------------------------------------------------------------
# Lattice dispersion in the SUBSTRATE rest frame (massless, 1D), c=1, a=1:
#   ω_s(k_s) = (2/a) sin(|k_s| a / 2)   (standard nearest-neighbour band)
#            = |k_s| (1 − (k_s a)²/24 + ...)   for small k_s
# ----------------------------------------------------------------------
def omega_sub(ks, a=1.0):
    return (2.0/a)*np.sin(np.abs(ks)*a/2.0)

# substrate-frame group (ray) velocity of a massless mode: u_s = dω_s/dk_s = cos(k_s a/2)
def u_sub(ks, a=1.0):
    return np.cos(ks*a/2.0)


# ----------------------------------------------------------------------
# (A) Preferred-frame leak. A wave packet with substrate-frame ray velocity u_s
#     has lab-frame ray velocity (boost β along the axis, emergent-Lorentz observer):
#        u_lab = (u_s − β)/(1 − u_s β)      [relativistic velocity addition of the worldline]
#     Forward photon: k_s>0, u_s=+cos(k_s a/2). Backward: k_s<0, u_s=−cos(k_s a/2).
#     Deficit c−|u_lab| measures the LIV. Fit its power in k_s and its β-asymmetry.
# ----------------------------------------------------------------------
def u_lab(u_s, beta):
    return (u_s - beta)/(1.0 - u_s*beta)

def part_A():
    print("(A) Boost the even-k lattice dispersion to a frame moving at β (preferred-frame leak)")
    a=1.0
    ks = np.array([2e-3, 4e-3, 8e-3, 1.6e-2, 3.2e-2])   # small vs BZ (π/a), EFT regime
    beta = 1.23e-3   # Earth vs CMB (~370 km/s)
    print(f"    β = {beta:.3e} (Earth vs CMB);  c=a=1; lattice (Planck) spacing a=1")
    print(f"    {'k/(π/a)':>10} {'δv_fwd=c−|v|(+k)':>18} {'δv_bwd=c−|v|(−k)':>18} {'(fwd−bwd)/mean':>16}")
    rowsA=[]
    for k in ks:
        us_f = +u_sub(k, a)            # forward photon
        us_b = -u_sub(k, a)            # backward photon
        vf = u_lab(us_f, beta)         # lab ray velocity (signed)
        vb = u_lab(us_b, beta)
        df, db = 1.0-abs(vf), 1.0-abs(vb)   # speed deficits (c=1)
        asym = (df-db)/(0.5*(df+db))
        print(f"    {k/np.pi:>10.4e} {df:>18.6e} {db:>18.6e} {asym:>16.4e}")
        rowsA.append(dict(k=float(k), dv_fwd=float(df), dv_bwd=float(db), asym=float(asym)))
    # power law of δv vs k (isotropic part): fit log-log
    dv_iso=np.array([0.5*(r['dv_fwd']+r['dv_bwd']) for r in rowsA])
    p_iso=np.polyfit(np.log(ks),np.log(dv_iso),1)[0]
    # anisotropy magnitude vs β: asym/β should be O(1) and ~const (n=2 modulated by β)
    asym_over_beta=np.mean([r['asym'] for r in rowsA])/beta
    print(f"    => δv power-law in k:  slope = {p_iso:.3f}  (==2 ⇒ n=2 / quadratic LIV)")
    print(f"    => forward/backward asymmetry / β = {asym_over_beta:.3f} (O(1) ⇒ modulation is")
    print(f"       n=2 × O(β), NOT an effective n=1 term: a boost cannot lower the power of a).")
    print(f"       LEAK SEALED: preferred frame ⇒ anisotropic n=2 (∝ β·(E/E_Pl)²), weaker than")
    print(f"       the isotropic n=2 already ruled unreachable.")
    return dict(rows=rowsA, dv_power_in_k=p_iso, asym_over_beta=asym_over_beta, beta=beta)


# ----------------------------------------------------------------------
# (B) Tree-level pattern universality: group velocity vs k for the shared
#     substrate. All low-k excitations → c; the spread is O((ka)²) = n=2.
#     A "species" with characteristic scale k_char has mean speed deficit
#     ~ (k_char a)². For real particles k_char a ~ m/E_Pl ⇒ (m/E_Pl)².
# ----------------------------------------------------------------------
def part_B():
    print("\n(B) Tree-level species universality (one substrate, one pattern-unaware cutoff)")
    a=1.0
    for kchar,label in [(1e-3,"k_char=1e-3 (E/E_Pl)"),(1e-2,"k_char=1e-2"),(1e-1,"k_char=1e-1")]:
        vg=u_sub(kchar,a)              # group velocity = cos(k a/2)
        deficit=1.0-vg
        print(f"    {label:>26}:  c−v_g = {deficit:.3e}  ≈ (k a)²/8 = {(kchar*a)**2/8:.3e}  (n=2)")
    print("    => all species share the same c at low energy; differences are (m/E_Pl)²-suppressed.")
    print("       So the tree-level dim-4 'species-dependent c' channel is ALSO Planck-suppressed.")
    return dict(note="tree-level dim-4 species-LIV ~ (m/E_Pl)^2, safe")


# ----------------------------------------------------------------------
# (C) The operator-dimension hierarchy: what even-k forbids vs what it doesn't.
# ----------------------------------------------------------------------
def part_C():
    print("\n(C) Operator-dimension hierarchy — what the even-k symmetry actually protects")
    rows=[
        # dim, parity in k, forbidden by even-k?, current bound (in E_Pl units or fractional), reachable?
        ("dim-5 (n=1, linear LIV, Myers-Pospelov)", "ODD in k",  "YES (forbidden)",
         "E_QG,1 > 5.9 E_Pl (LHAASO GRB221009A)", "would be REFUTED if present"),
        ("dim-6 (n=2, quadratic LIV, time-of-flight)", "EVEN in k", "no (allowed)",
         "E_QG,2 > 6e-8 E_Pl", "~1e7 below reach (structurally untestable)"),
        ("dim-4 (species-dependent c / c-anisotropy, SME c_μν)", "EVEN in k", "NO — NOT forbidden",
         "|c_LIV| < ~1e-22 (Hughes-Drever / cavity / clocks)", "REACHABLE — refutation-capable"),
    ]
    print(f"    {'operator':>52} | {'k-parity':>9} | {'even-k forbids?':>16} | reach")
    for d,par,forb,bound,reach in rows:
        print(f"    {d:>52} | {par:>9} | {forb:>16} | {bound} -> {reach}")
    print("\n    KEY POINT: the even-k symmetry forbids only the ODD-in-k (n=1, dim-5) term.")
    print("    A dim-4 LIV operator (different limiting speed) is EVEN in k — NOT forbidden — and")
    print("    is bounded at ~1e-22, i.e. NOT Planck-suppressed and REACHABLE. Tree-level")
    print("    universality (B) keeps it ~(m/E_Pl)², but radiative corrections in an interacting")
    print("    discrete theory generically percolate Planck LIV up to dim-4 at O(α/π) unless a")
    print("    symmetry protects (Collins-Perez-Sudarsky-Gambini-Pullin, PRL 93, 191301, 2004).")
    print("    The even-k/cubic symmetry does NOT protect dim-4 species-LIV. So door #2 is NOT")
    print("    'structurally untestable' — it has a reachable sub-channel the three-lock omits.")
    return rows


def main():
    print("="*78)
    print("Phase-12: stress-test 'door #2 (LIV) is symmetry-protected untestable'")
    print("="*78,"\n")
    A=part_A()
    B=part_B()
    C=part_C()
    print("\nConclusion")
    print("-"*78)
    print("1. The proposal's named leak is SEALED: the preferred frame gives an anisotropic n=2")
    print("   modulation (∝ β), never an effective n=1. The time-of-flight channel is uniformly")
    print("   n=2 and structurally unreachable — the three-lock holds on its own terms.")
    print("2. But the three-lock is scoped to TIME-OF-FLIGHT (dim≥5). It omits the dim-4 channel")
    print("   (species-dependent limiting speeds / c-anisotropy), which is EVEN in k (so NOT")
    print("   forbidden by the even-k symmetry it relies on) and bounded at ~1e-22 (REACHABLE).")
    print("3. Door #2 is therefore NOT 'structurally untestable'. It is untested at the one place")
    print("   it could be REFUTED — contingent on whether the framework's (unspecified) interacting")
    print("   UV completion radiatively percolates LIV into dim-4. 'Untestable, hence safe' is the")
    print("   neat conclusion that absorbed the challenge; the honest status is 'falsifiable at")
    print("   dim-4, pending a UV completion the framework still owes'. Bucket 0 unchanged (0).")
    out=dict(part_A_leak_sealed=A, part_B_universality=B,
             part_C_hierarchy=[dict(operator=r[0],k_parity=r[1],evenk_forbids=r[2],
                                    bound=r[3],reach=r[4]) for r in C],
             conclusion="leak sealed (n=2×β); but dim-4 species-LIV is even-k-allowed, "
                        "reachable (~1e-22), radiatively generic -> door2 falsifiable at dim-4, "
                        "not structurally untestable; pending UV completion")
    os.makedirs("simulations/results",exist_ok=True)
    with open("simulations/results/phase12_liv_door2_dim4_channel_result.json","w") as f:
        json.dump(out,f,indent=2)
    print("\nsaved -> simulations/results/phase12_liv_door2_dim4_channel_result.json")


if __name__=="__main__":
    main()
