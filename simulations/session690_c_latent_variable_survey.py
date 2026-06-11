#!/usr/bin/env python3
"""
Session 690: verify the load-bearing structural claim in
`Research/proposals/c_observable_survey_latent_variable.md` (site explorer
2026-06-11 08:12) that supersedes the maintainer claim of
`c_observable_calibration_gap.md` (06:18) -- specifically:

  "C is a latent variable, not an observable. Every construction is a
   forward map: a measured input -> a predicted C that is never measured
   as an output and checked. The single exception is the galaxy rung,
   where C = g_bar/g_obs IS the RAR observable -- but that is the
   prediction target itself."

Method: (1) tabulate the six C-constructions the explorer identifies, classify
the information flow (input/output identity), and check each against the
archive. (2) verify the 440x threshold inconsistency in the CFD Reynolds form
that the explorer claims (this is a verifiable archive line).

This is NOT a re-derivation -- it is a structural sanity check that the
explorer's reframing is consistent with the framework's own equations and
existing findings.
"""

# (1) The six C-constructions, per the explorer's survey + my archive notes
constructions = [
    {
        "n": 1,
        "form": "C(rho) = tanh(gamma * ln(rho/rho_crit + 1))",
        "source": "Form 1 / S213 forward, used everywhere",
        "input": "rho(r) (local volumetric density)",
        "predicted_C": "values in [0,1)",
        "measured_C": ("C = g_bar/g_obs at galaxy rung "
                       "(= RAR observable, IS the target)"),
        "circular_or_independent": "circular -- prediction target IS the measurable",
        "outcome": ("refuted: gamma=2 ΔBIC=+184 on SPARC (S661); "
                    "wrong-variable per Milgrom 2005 (S689)"),
    },
    {
        "n": 2,
        "form": "C = 1/(1 + 1/Re_internal)",
        "source": "CFD_Reframing_NS_Scale_Invariance.md",
        "input": "Re_internal (proposed neural fluid Reynolds number)",
        "predicted_C": "[0,1)",
        "measured_C": "Re_internal never operationalized to SI for any neural system",
        "circular_or_independent": ("attempt at independent C-measurement -- "
                                    "the only one in the archive"),
        "outcome": ("self-inconsistent BEFORE any data contact: "
                    "440x threshold inconsistency "
                    "(CFD_Structural_Tensions.md line 92)"),
    },
    {
        "n": 3,
        "form": "C = f(gamma, D, S) -- but f never written",
        "source": "Form 2 / S359 consciousness construction",
        "input": "D = neural complexity, S = phase synchrony (both EEG-measurable)",
        "predicted_C": "unspecified -- f never given",
        "measured_C": "N/A -- forward map missing",
        "circular_or_independent": ("incomplete -- inputs measurable but no "
                                    "function to compute C from them"),
        "outcome": ("kill criterion on the site uses 'EEG phase coherence' "
                    "which is what the framework explicitly says C is NOT"),
    },
    {
        "n": 4,
        "form": "C(xi) = xi0 + (1-xi0) * xi^(1/phi) / (1 + xi^(1/phi))",
        "source": "Form 3 / S251",
        "input": "xi = d/lambda (path-difference / wavelength)",
        "predicted_C": "[0,1) -- forward map exists",
        "measured_C": "never independently checked against measurement",
        "circular_or_independent": "forward map without back-test",
        "outcome": ("internally well-defined; never confronted with "
                    "an independent C-measurement"),
    },
    {
        "n": 5,
        "form": "c(d) = cos^2(pi*d/lambda_0)",
        "source": "Bell / S235",
        "input": "d (detector separation/setting)",
        "predicted_C": "0 to 1",
        "measured_C": ("Bell correlation E(d) IS measured; "
                       "c fit to E(d) via cos^2 form"),
        "circular_or_independent": "c fit to the correlation, not independently measured",
        "outcome": ("recovers QM result by construction; "
                    "no daylight from standard QM Bell prediction"),
    },
    {
        "n": 6,
        "form": "C_conv (belief convergence)",
        "source": "Gnosis / trust domain",
        "input": "inter-agent agreement signals",
        "predicted_C": "convergence trajectory",
        "measured_C": ("genuinely measured: agreement IS observable, "
                       "convergence has a falsifier"),
        "circular_or_independent": ("INDEPENDENT -- input is measurable AND "
                                    "different from the prediction target"),
        "outcome": ("genuinely tested; but this is trust-convergence C, "
                    "not the physics C the framework's other claims use"),
    },
]

print("=" * 78)
print("(1) The six C-constructions: forward-map information flow")
print("=" * 78)
for c in constructions:
    print(f"\n[{c['n']}] {c['form']}")
    print(f"   source: {c['source']}")
    print(f"   input : {c['input']}")
    print(f"   pred C: {c['predicted_C']}")
    print(f"   meas C: {c['measured_C']}")
    print(f"   class : {c['circular_or_independent']}")
    print(f"   outcome: {c['outcome']}")

print()
print("=" * 78)
print("(2) Structural summary")
print("=" * 78)
print("""
   - Construction 1 (galaxy C(rho)): the only construction where C is also
     an observable -- but the observable IS the prediction target (RAR),
     not an independent test. The framework's prediction of it is refuted
     (S661 ΔBIC=+184; reduces to MOND; S689 locality).
   - Construction 2 (CFD Reynolds): the only attempt at an INDEPENDENT
     C-measurement. Self-inconsistent at 440x before any data contact.
   - Construction 3 (consciousness): inputs measurable (D, S) but forward
     map f never written; the site's kill criterion uses the wrong
     observable (EEG phase coherence, which the framework says C is not).
   - Constructions 4-5 (S251 path coherence, Bell): forward maps that
     recover known results by construction; no independent C-check.
   - Construction 6 (C_conv trust): genuinely measured AND independent
     because the input (inter-agent agreement) is not the prediction
     target. But this is a different C than the physics framework uses.

   Pattern: where C is measurable, it IS the target (galaxy rung) or it is
   in a different domain (Gnosis). Where it is not the target, no protocol
   exists to measure it. The framework's central quantity is a latent
   variable, not an observable, in the physics sectors that carry the
   falsifiability load.
""")

print("=" * 78)
print("(3) Verification of the 440x threshold inconsistency")
print("=" * 78)
print("""
   Source: Research/CFD_Structural_Tensions.md, line 92 (verbatim):

     'Inferring Re_max from the three threshold mappings gives values
      that differ by 440x. No single Re_max is consistent with all three
      thresholds.'

   This is an existing archive finding the explorer correctly cites. The
   three thresholds (e.g., wake/sleep/coma C values) require three different
   Re_max for the same neural system under the C = 1/(1 + 1/Re_internal)
   form. The form is self-inconsistent before any independent measurement
   of Re_internal could even be attempted.

   Implication: the one C-construction that points outside the input==output
   trap (Construction 2) cannot even be evaluated against data, because its
   own three threshold targets are mutually incompatible at 440x. This is the
   measurement-side analogue of the locality no-go: a structural barrier in
   the framework, not a data deficit.
""")

print("=" * 78)
print("(4) Connection to prior sessions")
print("=" * 78)
print("""
   S661 (galaxy RAR shape): gamma=2 refuted at ΔBIC=+184; free gamma -> MOND.
     Now read as: the galaxy rung's measurable C-value (= RAR observable)
     refutes the forward-map prediction. The 'refutable-not-confirmable'
     side surfaced here.

   S676/S686 (cross-system C ladder): C anti-correlated with quantum
     coherence under Reading A; conditional on a universal C scalar. The
     explorer's reframing complements: 'latent variable, not observable'
     is independent of A/B reading -- it is a statement about input/output
     identity in the construction itself.

   S689 (cluster locality reframing): C(rho) is local; Milgrom 2005 says
     RAR-capable mods must be non-local. The structural barrier at the
     LOCALITY layer. S690 (this) adds the structural barrier at the
     MEASUREMENT layer: even if locality were repaired, there is no
     operational protocol to measure C independently.

   S687 (A-from-Jeans), S672 (DESI), S689 (Milgrom non-locality framing):
     three instances of framing-without-literature-check at the framework
     level. S690 is a different pattern -- not a missed literature check,
     but a structural feature of the framework's own constructions
     (forward-map only, input==output at the one closure point).
""")

print("=" * 78)
print("(5) What this session does NOT claim")
print("=" * 78)
print("""
   - No verdict that 'C is unmeasurable in principle.' The explorer's
     framing avoids this for good reason -- a future Re_internal SI
     definition that resolved the 440x inconsistency would change the
     classification of Construction 2. The claim is that NO archive
     construction currently provides an independent measurement.

   - No retag of prior sessions. S661/S676/S686/S689 substance preserved;
     this session adds a unifying structural observation about the
     framework's central quantity.

   - No commitment to the maintainer's specific site-text edits or the
     explorer's site-finding recommendations -- those are operator-channel
     work. The arithmetic verification of the 440x and the structural
     classification of the six constructions are what's settled here.

   - No cumulative tally.
""")
