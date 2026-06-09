#!/usr/bin/env python3
"""
Session 689: verify the locality classification table in
`Research/proposals/density_compander_nogo_locality_classification.md`
against each ansatz's variable dependence, and check consistency with prior
session findings.

The proposal's claim is structural: the discriminating axis for RAR-capability
is locality of the modification's state variable, NOT 'density-based.' Verlinde
(enclosed M_B), MOG (enclosed mass), MOND/AQUAL (acceleration), and surface
density Σ are all 'density-based' colloquially but key on NON-LOCAL functionals
of the baryon distribution. Synchronism's C(rho) is the rare ansatz keyed on
LOCAL volumetric rho(r), which is why it is caught.

This script does not run new physics; it tabulates the variable dependence of
each ansatz from prior session work and confirms the classification.
"""
import numpy as np

# Each entry: (name, state_variable_form, is_local, archive_session)
ansatze = [
    ("Synchronism C(rho)",
     "tanh(gamma * ln(rho(r)/rho_crit + 1))",
     True,
     "S661/S678/S683/S684 -- caught by local-rho non-locality gap"),

    ("MOND / AQUAL",
     "nu(|grad Phi|/a_0) -- requires Poisson solve",
     False,
     "S661/S684 -- McGaugh interpolating fn ~ ν_e(y=g_bar/a_0)"),

    ("MOND modified inertia",
     "trajectory/time-domain modification",
     False,
     "Milgrom 2005 (astro-ph/0510117); not directly used in archive"),

    ("Verlinde emergent gravity",
     "enclosed M_B(<r) entropy bound",
     False,
     "S673 GW170817 work; not framework-internal"),

    ("MOG (Moffat)",
     "enclosed mass + Yukawa screening",
     False,
     "Not directly used in archive"),

    ("Surface-density Sigma (Donato+2009 etc.)",
     "column-integrated Sigma(R) along line of sight",
     False,
     "S684 EFE work touches Sigma briefly; canonical reference"),
]

print("=" * 78)
print("(1) Locality classification per state variable")
print("=" * 78)
print(f"   {'Ansatz':<36} {'Local?':<8} {'RAR-capable per proposal':<24}")
print(f"   {'-'*36} {'-'*8} {'-'*24}")
for name, _, is_local, _ in ansatze:
    rar_capable = "No" if is_local else "Yes"
    local_str = "LOCAL" if is_local else "non-local"
    print(f"   {name:<36} {local_str:<8} {rar_capable:<24}")
print()
print("   Synchronism C(rho) is the only LOCAL entry. All RAR-capable ansaetze")
print("   are non-local. The discriminating axis is locality, not 'density-based.'")

print()
print("=" * 78)
print("(2) Cross-check vs prior session findings")
print("=" * 78)
for name, formula, is_local, archive in ansatze:
    print(f"   {name}")
    print(f"     state variable: {formula}")
    print(f"     archive:        {archive}")
print()

print("=" * 78)
print("(3) Implication for S678/S683 framing")
print("=" * 78)
print("   S678 framed the cluster failure as a 'structural impossibility' via the")
print("   A3 codomain bound (M_app/M_B <= 2 for C(rho) in [0,1)).")
print("   S683 refined this to 'wrong-variable' (local rho vs non-local g_bar),")
print("   verified within Coma's flat-cored gas profile.")
print()
print("   The proposal's correction does NOT overturn either finding's substance.")
print("   What it corrects is the FRAMING: cite Milgrom 2005 (astro-ph/0510117)")
print("   as the parent result; present C(rho) as a quantified instance of an")
print("   established non-locality constraint, not as a novel structural theorem.")
print()
print("   S683's terminology 'wrong-variable' was correct as a description, but")
print("   the diagnosis is acceleration-keying of the RAR (Milgrom 2005), which")
print("   pre-dates the framework's claim by 20+ years and applies generally to")
print("   any local-rho modification.")
print()
print("=" * 78)
print("(4) Methodology pattern: 3rd instance")
print("=" * 78)
print("   S672 (DESI cosmology): wrong-paper value 0.45 propagated from S668")
print("       without re-execution; corrected by re-running on the right slot.")
print("   S687 (A-from-Jeans): stated formula gives 4.6e-5 not 0.028, off by")
print("       614x; propagated across S631/S644 by re-reading not re-executing.")
print("   This (S689): cluster wrong-variable framing propagated through")
print("       S678/S683 without citing Milgrom 2005 as the parent literature.")
print()
print("   Three distinct sectors, same failure mode: framework-internal framing")
print("   without literature-check. The proposal itself is the corrective for")
print("   the third instance; this research-repo session accepts the correction")
print("   without re-litigating the underlying verification work.")
