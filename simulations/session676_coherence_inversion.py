#!/usr/bin/env python3
"""
Session 676: Verify the "coherence" naming inversion against the actual equation.

Proposal coherence_classicality_naming_and_test03_test05_double_filing.md (Problem 1)
+ explorer adjudication claim: C = tanh(gamma * ln(rho/rho_crit + 1)) with
gamma = 2/sqrt(N_corr) assigns LOW C to the most collective / most
quantum-phase-coherent systems (BEC, BCS superconductors), because gamma DECREASES
with N_corr. So "coherence" C is anti-correlated with (a) quantum phase coherence
(ODLRO, maximized in BEC/BCS) and (b) synchronization (a marching band and a BEC
are both maximally collective). Verify by COMPUTATION, not assertion (S669/S675).

C depends on density rho AND collectivity N_corr. To isolate the collectivity axis
(the naming claim), hold rho/rho_crit fixed and vary N_corr across a ladder of
systems ordered by ACTUAL macroscopic quantum coherence / synchronization.
"""
import numpy as np

def C(gamma, rho_ratio):
    return np.tanh(gamma * np.log(rho_ratio + 1.0))

def gamma_of(Ncorr):
    return 2.0 / np.sqrt(Ncorr)

# Ladder ordered by INCREASING macroscopic quantum coherence / collectivity.
# (N_corr = number of correlated DOF; the framework's own collectivity measure.)
ladder = [
    ("lone electron",        1,      "minimal collectivity; single-particle"),
    ("small molecule",       1e1,    "few correlated DOF"),
    ("nanoparticle",         1e3,    "mesoscopic"),
    ("BCS superconductor",   1e8,    "Cooper-pair condensate (ODLRO, maximally quantum-coherent)"),
    ("BEC",                  1e23,   "macroscopic wavefunction (ODLRO, maximally quantum-coherent)"),
]

for rho_ratio in [10.0, 1.0, 100.0]:
    print("=" * 76)
    print(f"C = tanh(gamma*ln(rho/rho_crit+1)), gamma=2/sqrt(N_corr), rho/rho_crit={rho_ratio:g}")
    print("=" * 76)
    print(f"   {'system':<20} {'N_corr':>8} {'gamma':>10} {'C':>10}  {'note'}")
    Cs = []
    for name, Nc, note in ladder:
        g = gamma_of(Nc)
        c = C(g, rho_ratio)
        Cs.append(c)
        print(f"   {name:<20} {Nc:>8.0e} {g:>10.2e} {c:>10.4f}  {note}")
    # quantum coherence increases DOWN the ladder; check C's direction
    monotone_down = all(Cs[i] >= Cs[i+1] for i in range(len(Cs)-1))
    print(f"   C decreases monotonically as quantum coherence/collectivity increases: "
          f"{monotone_down}")
    print(f"   lone electron C={Cs[0]:.4f}  -->  BEC C={Cs[-1]:.2e}")

print()
print("=" * 76)
print("Why it is structural (not a tuning artifact)")
print("=" * 76)
print("   gamma = 2/sqrt(N_corr) is strictly DECREASING in N_corr.")
print("   C = tanh(gamma * ln(rho/rho_crit+1)) is strictly INCREASING in gamma")
print("   (for rho>rho_crit, ln(.)>0, tanh increasing). Therefore, at any fixed")
print("   density, C is strictly DECREASING in N_corr. More collective => lower C.")
print("   Maximally quantum-coherent macroscopic systems (BEC/BCS) have the largest")
print("   N_corr, hence the smallest C. The inversion is forced by the equation.")

print()
print("=" * 76)
print("Three things the name 'coherence' / 'Synchronism' invoke, vs what C does")
print("=" * 76)
print("   (1) Quantum phase coherence (ODLRO): MAX in BEC/BCS -> but C ~ 0 there.")
print("       => C is ANTI-correlated with quantum coherence.")
print("   (2) Synchronization (lockstep): a marching band and a BEC are both")
print("       maximally collective/synchronized -> both large N_corr -> both low C.")
print("       => C is ANTI-correlated with synchronization (the framework's namesake).")
print("   (3) The framework's DEFINING metaphor (landing page/glossary): 'electrons")
print("       in a superconductor, in lockstep = high coherence' -- but a")
print("       superconductor (large N_corr -> tiny gamma) is pinned at C~0 by the")
print("       equation. The exemplar used to DEFINE high coherence is a low-C system.")
print()
print("   What C actually is (per the equation): a density-driven saturation index,")
print("   high for DENSE & weakly-correlated matter, low for SPARSE *or* strongly-")
print("   correlated matter. It carries no hbar / T / action / decoherence rate")
print("   (S613: no decoherence parameter) -> zero quantum-vs-classical content.")
print("   Renaming to 'classicality'/'decoherence-fraction' would assert a FALSEHOOD")
print("   (that a lone electron is 'maximally classical', a BEC 'maximally quantum').")
print("   The honest fix is a SCOPE statement, not a vocabulary swap.")
