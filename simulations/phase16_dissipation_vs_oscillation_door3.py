#!/usr/bin/env python3
"""
Phase-16 (door #3, secular/dissipative sector) — the dissipation that would make a novel
prediction is in the substrate that cannot host the matter to test it.

Context
-------
The recalibration (dp, 2026-06-24) reopened door #3 — the secular / time-domain /
non-equilibrium sector — as the place least examined by the arc and most distinctive for
the framework's defining commitments (absolute time + a *dissipative* substrate on a global
tick, S666). PREDICTIONS.md lists "dissipation / arrow-of-time" among door-#3 candidates.

The natural door-#3 novel prediction from a DISSIPATIVE substrate is **intrinsic
decoherence**: even a perfectly isolated pattern loses coherence to the substrate's
dissipation, at a rate set by the substrate diffusion constant. Standard QM is unitary (no
intrinsic decoherence), so this is genuinely distinctive — exactly the kind of thing door #3
is supposed to harbor. This realizes the session's standing tension #4 (the oscillation
basis and the dissipative C(rho)/diffusion account may be *competing*, not harmonious).

But the framework already established (S617/S665/S666; Phase-1; FUNDAMENTALS §3):
  - the substrate's own rule is 1-DOF scalar DIFFUSION — irrotational, DISSIPATIVE — and it
    **cannot host particle dynamics** (it disperses/dissipates at 0% pass rate);
  - stable matter required switching to a *different*, CONSERVATIVE substrate: a 2nd-order
    wave equation with a focusing nonlinearity (a discrete breather / Klein-Gordon soliton).
These are "two distinct substrates connected only by narrative" (PREDICTIONS Bucket 2).

This script demonstrates the resulting obstruction concretely:
  (A) DISSIPATIVE substrate (diffusion ∂I/∂t = D ∂²I/∂x²): carries the door-#3 signature —
      intrinsic decoherence, mode k decaying at rate Γ(k)=2Dk² (verified) — but CANNOT hold
      a localized oscillating pattern (no matter: a bump spreads to uniform).
  (B) CONSERVATIVE substrate (wave ∂²I/∂t² = c²∂²I/∂x², the matter substrate's linear core):
      hosts persistent patterns and conserves energy to machine precision — NO intrinsic
      decoherence (the door-#3 signature is absent).
  (C) The exclusivity: a single linear substrate is either dissipative (signature, no matter)
      or conservative (matter, no signature). The door-#3 novel prediction lives in the
      substrate that cannot host the matter that would exhibit it.
  (D) Magnitude, if the framework committed to dissipative-fundamental: Γ=2Dk² with a
      Planck-scale diffusion D~c·ℓ_Pl, for lab superpositions, vs CSL / matter-wave bounds —
      a real but unpinned door-#3 bet, blocked by the matter obstruction.

numpy only. Headless. Writes results JSON.
Author: CBP-Claude (Opus 4.8), autonomous — treads the recalibration's door #3 + tension #4.
"""
import json, os
import numpy as np


def main():
    N = 512
    L = 2*np.pi
    x = np.linspace(0, L, N, endpoint=False)
    dx = x[1]-x[0]
    k = 2*np.pi*np.fft.fftfreq(N, d=dx)     # wavenumbers
    D = 0.1                                   # diffusion constant (toy units)
    c = 1.0                                   # wave speed

    print("="*78)
    print("Phase-16 (door #3): dissipation gives the signature but not the matter")
    print("="*78)

    # ---- (A) DISSIPATIVE substrate: intrinsic decoherence, but no matter ----
    print("\n(A) DISSIPATIVE substrate  ∂I/∂t = D ∂²I/∂x²   (D=%.2f) — the door-#3 signature" % D)
    # single-mode decay: Î(k,t) = Î(k,0) e^{-D k² t}; intensity mode decays at 2Dk²
    print("    mode-resolved decay of ∫I² (intrinsic decoherence rate Γ(k)):")
    print(f"    {'k':>6} {'measured Γ':>12} {'2·D·k²':>12} {'ratio':>8}")
    rowsA=[]
    for kk in [1.0, 2.0, 3.0, 4.0]:
        I0 = np.cos(kk*x)                      # pure mode k
        # evolve a short time, measure d/dt ln(∫I²)
        def intensity_at(t):
            Ihat = np.fft.fft(I0)*np.exp(-D*k**2*t)
            I = np.real(np.fft.ifft(Ihat))
            return np.sum(I**2)*dx
        t0, dt = 0.0, 1e-3
        gamma = -(np.log(intensity_at(t0+dt)) - np.log(intensity_at(t0)))/dt
        pred = 2*D*kk**2
        print(f"    {kk:>6.1f} {gamma:>12.5f} {pred:>12.5f} {gamma/pred:>8.4f}")
        rowsA.append(dict(k=kk, gamma=float(gamma), pred_2Dk2=float(pred)))
    # matter test: a localized bump spreads to uniform (no persistent pattern)
    bump = np.exp(-((x-np.pi)**2)/(2*0.2**2))
    def peak_after(t):
        Ih = np.fft.fft(bump)*np.exp(-D*k**2*t)
        return np.max(np.real(np.fft.ifft(Ih)))
    p0, p1, p2 = peak_after(0), peak_after(1.0), peak_after(10.0)
    print(f"    localized bump peak: t=0:{p0:.3f}  t=1:{p1:.3f}  t=10:{p2:.3f}  -> spreads to uniform")
    print("    => dissipative substrate HAS intrinsic decoherence (Γ=2Dk²) but CANNOT hold matter.")

    # ---- (B) CONSERVATIVE substrate: matter, but no decoherence ----
    print("\n(B) CONSERVATIVE substrate  ∂²I/∂t² = c²∂²I/∂x²  (matter substrate's linear core)")
    # state (I, Idot); spectral exact: Î(t)=Î0 cos(ckt) + (Îdot0/(ck)) sin(ckt)
    I0 = np.exp(-((x-np.pi)**2)/(2*0.4**2)) - np.mean(np.exp(-((x-np.pi)**2)/(2*0.4**2)))
    Idot0 = np.zeros_like(I0)
    Ih0 = np.fft.fft(I0); Idh0 = np.fft.fft(Idot0)
    ck = c*np.abs(k); ck[0]=1.0
    def energy_at(t):
        Iht = Ih0*np.cos(ck*t) + np.where(np.abs(k)>0, Idh0/ck*np.sin(ck*t), 0)
        Idht = -Ih0*ck*np.sin(ck*t) + Idh0*np.cos(ck*t)
        I = np.real(np.fft.ifft(Iht)); Idot = np.real(np.fft.ifft(Idht))
        Ix = np.real(np.fft.ifft(1j*k*Iht))
        return 0.5*np.sum(Idot**2 + c**2*Ix**2)*dx
    E = [energy_at(t) for t in [0,1,5,20,100]]
    print(f"    energy E(t) at t=0,1,5,20,100: " + " ".join(f"{e:.6f}" for e in E))
    print(f"    relative drift max|E-E0|/E0 = {max(abs(e-E[0]) for e in E)/E[0]:.2e}  (conserved)")
    print("    => conservative substrate HOLDS matter / conserves energy but has NO intrinsic")
    print("       decoherence: the door-#3 signature is absent.")

    # ---- (C) the exclusivity ----
    print("\n(C) Exclusivity: a single LINEAR substrate is dissipative XOR conservative.")
    print("    dissipative ⇒ intrinsic decoherence (door-#3 signature) ⇒ no stable matter;")
    print("    conservative ⇒ stable matter ⇒ no intrinsic decoherence.")
    print("    The framework's matter (Phase-1 wave+focusing) is conservative; its distinctive")
    print("    dissipative substrate (S666) cannot host that matter. So door #3's novel signature")
    print("    lives in the substrate that cannot make the matter that would exhibit it —")
    print("    the 'two substrates connected only by narrative' problem, now the SPECIFIC")
    print("    obstruction blocking door #3. (This realizes tension #4: the oscillation/persistent-")
    print("    pattern account and the dissipative/diffusion account are COMPETING, not harmonious.)")

    # ---- (D) magnitude IF dissipative-fundamental ----
    print("\n(D) Magnitude IF the framework committed to dissipative-fundamental (Planck diffusion):")
    c_si=2.998e8; lP=1.616e-35
    D_planck = c_si*lP                          # ~ c·ℓ_Pl, a natural Planck diffusion constant
    print(f"    D ~ c·ℓ_Pl = {D_planck:.2e} m²/s")
    for Lsys,label in [(1e-10,"atomic 1 Å"),(1e-9,"molecule 1 nm"),(1e-6,"meso 1 µm")]:
        kk = 1.0/Lsys
        Gamma = 2*D_planck*kk**2                 # intrinsic decoherence rate 1/s
        print(f"    superposition scale {label:>14}:  Γ=2D/L² = {Gamma:.2e} /s   (τ={1/Gamma:.2e} s)")
    print("    For comparison, CSL/matter-wave intrinsic-decoherence bounds sit near ~1e-8..1e-16 /s")
    print("    (model-dependent). So a Planck-diffusion door-#3 bet is small but NOT obviously")
    print("    Planck-*killed* at nm scale — it is a real, unpinned bet. BUT it is BLOCKED: this")
    print("    rate applies to a dissipative substrate that holds no matter, so there is no")
    print("    isolated material superposition for it to act on within the framework as it stands.")

    print("\nConclusion")
    print("-"*78)
    print("Door #3's most natural novel prediction — intrinsic decoherence ∝ D·k² from the")
    print("dissipative substrate — is genuinely distinctive (standard QM is unitary) and lands in")
    print("the least-examined, most-snapshot-robust sector. But it is obstructed from WITHIN: the")
    print("dissipative substrate that produces the signature cannot host stable matter, and the")
    print("conservative wave+focusing substrate that hosts matter produces no signature. Door #3")
    print("therefore cannot deliver a novel prediction until the framework UNIFIES the two")
    print("substrates (or commits to dissipative-fundamental with matter as quasi-stable solitons")
    print("that resist diffusion — Phase-1 showed focusing, a CONSERVATIVE mechanism, is what makes")
    print("patterns, so this is open work). Not a structural closure (respecting the recalibration):")
    print("a SPECIFIC, nameable obstruction + the constructive condition that would turn door #3")
    print("into a real Bucket-1 bet. Bucket 0 unchanged (0).")

    out=dict(
        A_decoherence_rates=rowsA, A_bump_peaks=dict(t0=float(p0),t1=float(p1),t10=float(p2)),
        B_energy=[float(e) for e in E],
        B_energy_drift=float(max(abs(e-E[0]) for e in E)/E[0]),
        D_planck_diffusion_m2_s=float(D_planck),
        conclusion="dissipative substrate gives intrinsic-decoherence signature (Gamma=2Dk^2) but "
                   "no matter; conservative wave+focusing substrate gives matter but no signature; "
                   "door-#3 novel prediction is blocked by the two-substrate split (realizes "
                   "tension #4: oscillation vs dissipation are competing not harmonious). Specific "
                   "obstruction + constructive condition, not a structural closure. Bucket 0 = 0.")
    os.makedirs("simulations/results",exist_ok=True)
    with open("simulations/results/phase16_dissipation_vs_oscillation_door3_result.json","w") as f:
        json.dump(out,f,indent=2)
    print("\nsaved -> simulations/results/phase16_dissipation_vs_oscillation_door3_result.json")


if __name__=="__main__":
    main()
