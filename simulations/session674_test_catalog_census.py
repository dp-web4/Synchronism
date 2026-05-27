#!/usr/bin/env python3
"""
Session 674: Complete census of the Experimental Test Catalog (24 tests) against the
two discriminators developed this arc:
  - DERIVED AMPLITUDE? (S673): is the predicted magnitude derived from first
    principles, or read off the data / given as an order-of-magnitude / calibrated?
  - TIER-1? (S670): central object === a pre-existing quantity (relabel) vs a genuine
    new derived prediction.

Honors "unconfirmed != wrong": distinguishes EXECUTED-AND-COLLAPSED (refuted/
disfavored/degenerate by actual analysis, with session ref) from UNTESTED (nobody ran
it). For untested tests with a claimed specific number, records whether the number's
derived-status is VERIFIED-not-derived, or merely UNVERIFIED by me.

NOTE: the site/proposal numbering (TEST-17 'cluster gamma-gradient', TEST-21 'BAO
sub-peaks') does NOT match this archive catalog (TEST-17 = scale-dependent c,
TEST-21 = entanglement across scales). Numbering is inconsistent across documents.
"""

# (id, short, category, derived_amplitude, note)
# category: EXECUTED_COLLAPSED | SELF_DEGENERATE | NO_DERIVED_AMP | UNTESTED_FRONTIER
# derived_amplitude: "no(executed)" | "no(self-admit)" | "no(order-of-mag)"
#                    | "UNVERIFIED" | "self-flagged-coincidence"
tests = [
 ("01","TDG age-DM",            "UNTESTED_FRONTIER","UNVERIFIED",
    "tau~1.6 Gyr, 0.3/Gyr likely calibrated; standard says ~0 DM in TDGs; untested"),
 ("02","UDG max DM",            "SELF_DEGENERATE","no(self-admit)",
    "catalog: MEDIUM; LCDM also predicts high f_DM for low-mass"),
 ("03","compact-E min DM",      "SELF_DEGENERATE","no(self-admit)",
    "catalog: MEDIUM; tests C(rho)->1 saturation; LCDM also low f_DM dense"),
 ("04","BAO/fsigma8 (04a)",     "EXECUTED_COLLAPSED","no(executed)",
    "S668/S672: disfavored ~2sigma, kill triggered, post-hoc"),
 ("05","CMB coldspot-density",  "SELF_DEGENERATE","no(self-admit)",
    "catalog: MEDIUM; ISW degenerate"),
 ("06","variable alpha_em",     "UNTESTED_FRONTIER","UNVERIFIED",
    "beta~1e-5 order-of-mag; spatial-not-temporal; untested"),
 ("07","cosmic interference 500Mpc","UNTESTED_FRONTIER","UNVERIFIED",
    "catalog touts VERY HIGH/unique; 500 Mpc asserted not derived (quick check); untested"),
 ("08","SPARC environment",     "EXECUTED_COLLAPSED","no(executed)",
    "S637: signal 120x below SPARC floor; S654: MOND+EFE degenerate"),
 ("09","photosynth coherence",  "SELF_DEGENERATE","no(self-admit)",
    "catalog: MEDIUM, 'not uniquely predicted by Synchronism'"),
 ("10","enzyme KIE-gamma",      "SELF_DEGENERATE","no(self-admit)",
    "catalog: LOW-MEDIUM, 'KIE correlates with many things'"),
 ("11","EEG anesthesia Phi=3.5","UNTESTED_FRONTIER","UNVERIFIED",
    "IIT also predicts a threshold; 3.5 derived-status unverified; untested"),
 ("12","qubit optimal C*=0.79", "UNTESTED_FRONTIER","self-flagged-coincidence",
    "framework's OWN OQ doc: 0.79 vs 0.5 'without clear relationship'; 0.79 also "
    "appears as an r-coefficient & 'mean C' elsewhere; untested"),
 ("13","circadian gamma",       "SELF_DEGENERATE","no(self-admit)",
    "catalog: MEDIUM; many things show circadian variation"),
 ("14","wide binary (Gaia)",    "EXECUTED_COLLAPSED","no(executed)",
    "MOND-degenerate (threshold vs density); galaxy sector closed S654/S661"),
 ("15","GW speed-DM column",    "EXECUTED_COLLAPSED","no(executed)",
    "S673: alpha read off GW170817, no derived amplitude; Case-3 -> remove"),
 ("16","BH ringdown shift",     "NO_DERIVED_AMP","no(order-of-mag)",
    "delta~1e-4 to 1e-5 range given, not derived"),
 ("17","scale-dependent c",     "UNTESTED_FRONTIER","UNVERIFIED",
    "km/s deviations provenance unverified; CONTRADICTS S667 (parabolic) & S641 "
    "(Lorentz gap); untested"),
 ("18","hot superconductor",    "EXECUTED_COLLAPSED","no(self-admit)",
    "catalog SELF-ADMITS: Tc formula wrong, eta = standard Abrikosov-Gor'kov, "
    "'not a unique Synchronism prediction' (S616 audit)"),
 ("19","microtubule coherence", "SELF_DEGENERATE","no(self-admit)",
    "catalog: MEDIUM; Orch-OR also predicts MT coherence"),
 ("20","consciousness Phi-scaling","UNTESTED_FRONTIER","UNVERIFIED",
    "same Phi_crit=3.5 as TEST-11; cross-species; untested"),
 ("21","entanglement across scales","UNTESTED_FRONTIER","UNVERIFIED",
    "claims Bell S>2.5 vs QFT S<2.1; S>2.5 asserted not derived; untested"),
 ("22","virus decoherence tau", "UNTESTED_FRONTIER","UNVERIFIED",
    "tau~1e6 s vs Penrose 1e3, std 1e10; derived-status unverified; in-progress Vienna"),
 ("23","SGWB anisotropy",       "NO_DERIVED_AMP","no(order-of-mag)",
    "anisotropy 'follows DM' but no amplitude; needs LISA/ET"),
 ("24","void expansion rate",   "NO_DERIVED_AMP","no(order-of-mag)",
    "epsilon~1e-3 order given, not derived"),
]

from collections import Counter
cats = Counter(t[2] for t in tests)
print("=" * 78)
print("CENSUS of the 24-test Experimental Test Catalog")
print("=" * 78)
for cat in ["EXECUTED_COLLAPSED","SELF_DEGENERATE","NO_DERIVED_AMP","UNTESTED_FRONTIER"]:
    rows = [t for t in tests if t[2]==cat]
    print(f"\n## {cat}  (n={len(rows)})")
    for tid, short, _, da, note in rows:
        print(f"   TEST-{tid} {short:<26} [amp:{da:<22}] {note}")

print()
print("=" * 78)
print("TALLY")
print("=" * 78)
for cat,n in cats.most_common():
    print(f"   {cat:<22} {n}")
n_exec = cats["EXECUTED_COLLAPSED"]
n_closed = n_exec + cats["SELF_DEGENERATE"] + cats["NO_DERIVED_AMP"]
n_front = cats["UNTESTED_FRONTIER"]
print(f"\n   Executed-and-collapsed: {n_exec}  (ALL collapsed; confirmed-discriminating = 0)")
print(f"   Effectively closed (executed + self-degenerate + no-derived-amp): {n_closed}/24")
print(f"   Genuinely UNTESTED frontier (unconfirmed != wrong): {n_front}/24")
print()
print("   Of the untested frontier, derived-amplitude status:")
print("     - 0 verified as first-principles-derived")
print("     - 1 self-flagged by the framework as an unexplained coincidence (C*=0.79)")
print("     - 1 contradicts the framework's own substrate analysis (TEST-17 vs S667/S641)")
print("     - rest: order-of-magnitude or asserted; derived-status UNVERIFIED by me")
print()
print("   HONEST VERDICT: every test that has been EXECUTED collapsed (import/")
print("   degenerate/no-derived-amplitude). The catalog is NOT fully closed --")
print("   ~8/24 remain untested. But the quick provenance check found NO untested")
print("   test with a clean derived amplitude; the live shots share the exact")
print("   structural property (no derived amplitude) that has collapsed every")
print("   executed test. Whether they have derived amplitudes is the open question")
print("   that would settle sterile-vs-generative (S671) -- answerable per-test by")
print("   checking provenance. Novel CONFIRMED-discriminating count: 0, by execution.")
