#!/usr/bin/env python3
"""
Session #300: Experimental Validation Protocol for η Framework
Hot Superconductor Arc (Session 5/?)

This milestone session designs experimental protocols to validate the
reachability factor (η) framework developed in Sessions #292, #297-299.

Key validation targets:
1. Direct η measurement via multiple techniques
2. Verification of η-T_c relationship
3. Testing specific material predictions
4. Falsification criteria

Experimental techniques:
- NMR relaxation (T1, T2)
- Optical conductivity / reflectivity
- Angle-resolved photoemission (ARPES)
- Tunneling spectroscopy (STM/STS)
- Thermal transport (κ)
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
from enum import Enum

print("=" * 80)
print("SESSION #300: EXPERIMENTAL VALIDATION PROTOCOL")
print("Hot Superconductor Arc (Session 5/?)")
print("MILESTONE SESSION")
print("=" * 80)

# ============================================================================
# PART 1: η MEASUREMENT TECHNIQUES
# ============================================================================

print("\n" + "=" * 60)
print("PART 1: η MEASUREMENT TECHNIQUES")
print("=" * 60)

class Technique(Enum):
    NMR = "NMR Relaxation"
    OPTICAL = "Optical Conductivity"
    ARPES = "ARPES"
    TUNNELING = "Tunneling Spectroscopy"
    THERMAL = "Thermal Transport"
    PENETRATION = "Penetration Depth"

@dataclass
class MeasurementTechnique:
    """Experimental technique for η measurement"""
    name: Technique
    principle: str
    eta_sensitivity: str
    precision: float  # % error achievable
    temperature_range: str
    sample_requirements: str
    equipment_cost: str
    measurement_time: str
    advantages: List[str]
    limitations: List[str]

techniques = [
    MeasurementTechnique(
        name=Technique.NMR,
        principle="Nuclear spin relaxation rates (1/T1, 1/T2) probe quasiparticle dynamics. "
                  "η affects the spectral weight coupling to thermal excitations.",
        eta_sensitivity="1/T1 ∝ η × N(E_F)² × T for T > T_c/3",
        precision=5.0,
        temperature_range="4K - 300K",
        sample_requirements="Single crystal > 1mm, isotope enrichment helpful",
        equipment_cost="$200K-500K (NMR spectrometer)",
        measurement_time="1-7 days per temperature sweep",
        advantages=[
            "Bulk probe (not surface sensitive)",
            "Well-established technique",
            "Can distinguish local vs global effects",
            "Site-selective with different nuclei"
        ],
        limitations=[
            "Requires specific NMR-active nuclei",
            "Signal-to-noise can be challenging",
            "Skin depth effects in metals",
            "Interpretation requires careful modeling"
        ]
    ),
    MeasurementTechnique(
        name=Technique.OPTICAL,
        principle="Optical conductivity σ(ω) measures charge response. "
                  "η affects the thermal broadening of the Drude peak and "
                  "the spectral weight transfer to the gap.",
        eta_sensitivity="σ₁(ω,T) thermal broadening ∝ η × T",
        precision=3.0,
        temperature_range="4K - 300K",
        sample_requirements="Flat surface > 1mm², polished or cleaved",
        equipment_cost="$300K-800K (FTIR + cryostat)",
        measurement_time="2-5 days per complete spectrum",
        advantages=[
            "Direct probe of charge dynamics",
            "Can extract both η and Δ independently",
            "Frequency-resolved information",
            "Relatively non-destructive"
        ],
        limitations=[
            "Surface sensitive (skin depth ~100nm)",
            "Requires high-quality surfaces",
            "Kramers-Kronig analysis needed",
            "Phonon contributions must be subtracted"
        ]
    ),
    MeasurementTechnique(
        name=Technique.ARPES,
        principle="Angle-resolved photoemission measures k-resolved spectral function. "
                  "η manifests in the thermal broadening of the quasiparticle peaks.",
        eta_sensitivity="Γ(k,T) = Γ₀ + η × α(k) × T",
        precision=8.0,
        temperature_range="10K - 200K (limited by resolution)",
        sample_requirements="Cleaved single crystal in UHV",
        equipment_cost="$2M-10M (synchrotron access or lab system)",
        measurement_time="1-3 days per sample (beamtime)",
        advantages=[
            "Momentum-resolved information",
            "Can map η(k) across Fermi surface",
            "Direct visualization of gap structure",
            "Can see nesting effects directly"
        ],
        limitations=[
            "Extremely surface sensitive (~1nm)",
            "Requires UHV and cryogenic conditions",
            "Energy resolution limits low-T studies",
            "Matrix element effects complicate interpretation"
        ]
    ),
    MeasurementTechnique(
        name=Technique.TUNNELING,
        principle="Scanning tunneling spectroscopy (STS) measures local DOS. "
                  "Thermal smearing of coherence peaks indicates η.",
        eta_sensitivity="Peak width Γ = Γ₀ + η × 3.5 k_B T",
        precision=5.0,
        temperature_range="0.3K - 30K (with dilution fridge)",
        sample_requirements="Atomically flat surface, cleaved in situ",
        equipment_cost="$500K-2M (LT-STM system)",
        measurement_time="1-4 weeks for statistical mapping",
        advantages=[
            "Atomic-scale spatial resolution",
            "Can measure gap locally",
            "Vortex core spectroscopy possible",
            "No photon/phonon backgrounds"
        ],
        limitations=[
            "Extremely surface sensitive",
            "Very slow for large-area mapping",
            "Tip effects can complicate analysis",
            "Limited temperature range"
        ]
    ),
    MeasurementTechnique(
        name=Technique.THERMAL,
        principle="Thermal conductivity κ(T) in SC state probes nodal quasiparticles. "
                  "η affects the phonon-quasiparticle coupling.",
        eta_sensitivity="κ/T as T→0 depends on nodal structure; "
                        "η affects T-dependence of thermal activation",
        precision=10.0,
        temperature_range="0.1K - T_c",
        sample_requirements="Well-shaped sample for heat flow geometry",
        equipment_cost="$100K-300K (dilution fridge + heaters/sensors)",
        measurement_time="2-4 weeks per sample",
        advantages=[
            "Bulk probe",
            "Directly measures thermal coupling",
            "Sensitive to gap nodes",
            "Clean interpretation for d-wave"
        ],
        limitations=[
            "Phonon contribution dominates at high T",
            "Requires very low temperatures for clean signal",
            "Sample geometry critical",
            "Boundary scattering at lowest T"
        ]
    ),
    MeasurementTechnique(
        name=Technique.PENETRATION,
        principle="London penetration depth λ(T) measures superfluid density. "
                  "Temperature dependence encodes gap structure and η.",
        eta_sensitivity="Δλ(T) ∝ exp(-Δ/(η×k_B×T)) for s-wave; "
                        "∝ T for d-wave nodes",
        precision=2.0,
        temperature_range="0.1K - T_c",
        sample_requirements="Small crystals or thin films, no magnetic impurities",
        equipment_cost="$200K-500K (microwave cavity or mutual inductance)",
        measurement_time="1-2 weeks per sample",
        advantages=[
            "Very high precision possible",
            "Clean separation of gap and η effects",
            "Sensitive to symmetry changes",
            "Multiple techniques available"
        ],
        limitations=[
            "Surface preparation critical",
            "Magnetic impurities fatal",
            "Model-dependent extraction of η",
            "Requires careful calibration"
        ]
    ),
]

print("\nExperimental Techniques for η Measurement:")
print("-" * 100)
print(f"{'Technique':<20} {'Precision':<10} {'T Range':<15} {'Cost':<20} {'Time':<15}")
print("-" * 100)
for tech in techniques:
    print(f"{tech.name.value:<20} {tech.precision:.0f}%{'':<8} {tech.temperature_range:<15} {tech.equipment_cost:<20} {tech.measurement_time:<15}")

# ============================================================================
# PART 2: η EXTRACTION PROTOCOLS
# ============================================================================

print("\n" + "=" * 60)
print("PART 2: η EXTRACTION PROTOCOLS")
print("=" * 60)

@dataclass
class ExtractionProtocol:
    """Protocol for extracting η from experimental data"""
    technique: Technique
    observable: str
    model_equation: str
    fitting_procedure: str
    control_parameters: List[str]
    systematic_checks: List[str]

protocols = [
    ExtractionProtocol(
        technique=Technique.NMR,
        observable="1/T1 (spin-lattice relaxation rate)",
        model_equation="1/T1 = A × η × N(E_F)² × k_B × T × F(Δ/k_B T)",
        fitting_procedure="""
1. Measure 1/T1 from T_c down to T_c/5
2. Above T_c: extract N(E_F) from Korringa relation
3. Below T_c: fit to Hebel-Slichter coherence peak
4. Extract η from amplitude of thermal enhancement
5. Compare across multiple nuclei for consistency""",
        control_parameters=["Magnetic field strength", "Pulse sequence", "Sample orientation"],
        systematic_checks=["Field independence", "Isotope scaling", "Linewidth analysis"]
    ),
    ExtractionProtocol(
        technique=Technique.OPTICAL,
        observable="σ₁(ω, T) (real part of optical conductivity)",
        model_equation="σ₁(ω) = (σ₀/τ) / (1 + (ωτ)²) where τ⁻¹ = τ₀⁻¹ + η × α × T",
        fitting_procedure="""
1. Measure reflectivity R(ω) from FIR to UV
2. Kramers-Kronig transform to get σ₁(ω)
3. Fit Drude peak width vs temperature
4. Extract η from slope of Γ(T) = Γ₀ + η × α × T
5. Verify with spectral weight conservation""",
        control_parameters=["Photon energy range", "Polarization", "Sample orientation"],
        systematic_checks=["KK consistency", "Sum rules", "Phonon subtraction"]
    ),
    ExtractionProtocol(
        technique=Technique.ARPES,
        observable="Γ(k, T) (quasiparticle linewidth)",
        model_equation="Γ(k, T) = Γ_imp + Γ_e-e(k,T) + η × Γ_e-ph(k,T)",
        fitting_procedure="""
1. Measure EDCs at multiple k-points on Fermi surface
2. Fit quasiparticle peaks with Lorentzians
3. Extract Γ(k) as function of T at each k-point
4. Separate impurity, e-e, and e-ph contributions
5. η = d(Γ_thermal)/dT / (expected phonon contribution)""",
        control_parameters=["Photon energy", "Polarization", "Sample aging"],
        systematic_checks=["Energy resolution deconvolution", "Matrix element effects", "Surface vs bulk"]
    ),
    ExtractionProtocol(
        technique=Technique.PENETRATION,
        observable="Δλ(T) (change in penetration depth)",
        model_equation="Δλ/λ(0) = (π Δ / 2 k_B T)^(1/2) × exp(-Δ/(η × k_B T)) for s-wave",
        fitting_procedure="""
1. Measure λ(T) using tunnel diode oscillator or microwave cavity
2. Fit low-T exponential to extract Δ
3. Fit intermediate-T to extract η from thermal activation
4. For d-wave: analyze linear-T coefficient
5. η from deviation of λ(T) from theoretical curve""",
        control_parameters=["Frequency", "AC field amplitude", "Sample geometry"],
        systematic_checks=["Geometry corrections", "Demagnetization", "Surface impedance"]
    ),
]

print("\nη Extraction Protocols:")
print("-" * 80)
for protocol in protocols:
    print(f"\n{protocol.technique.value}:")
    print(f"  Observable: {protocol.observable}")
    print(f"  Model: {protocol.model_equation[:60]}...")
    print(f"  Controls: {', '.join(protocol.control_parameters[:3])}")

# ============================================================================
# PART 3: TESTABLE HYPOTHESIS MATRIX
# ============================================================================

print("\n" + "=" * 60)
print("PART 3: TESTABLE HYPOTHESIS MATRIX")
print("=" * 60)

@dataclass
class TestableHypothesis:
    """Specific testable prediction from η framework"""
    id: str
    hypothesis: str
    prediction: str
    measurement: str
    expected_value: str
    tolerance: str
    falsification: str
    priority: str  # High, Medium, Low
    estimated_cost: str

hypotheses = [
    TestableHypothesis(
        id="H300.1",
        hypothesis="η ordering across cuprate families",
        prediction="η(Hg-1223) < η(YBCO) < η(Bi-2212) < η(LSCO)",
        measurement="NMR 1/T1 or optical conductivity thermal broadening",
        expected_value="0.33, 0.38, 0.42, 0.51 respectively",
        tolerance="±0.05",
        falsification="If η ordering doesn't match T_c ordering",
        priority="High",
        estimated_cost="$50K-100K (synchrotron or NMR time)"
    ),
    TestableHypothesis(
        id="H300.2",
        hypothesis="η increases with pressure in pnictides",
        prediction="η(BaFe₂As₂) increases ~10% per GPa",
        measurement="High-pressure NMR T1 or optical reflectivity",
        expected_value="dη/dP ≈ 0.02 GPa⁻¹",
        tolerance="±0.005",
        falsification="If η decreases or stays constant under pressure",
        priority="Medium",
        estimated_cost="$100K-200K (DAC + NMR system)"
    ),
    TestableHypothesis(
        id="H300.3",
        hypothesis="SmFeAsO has lower η than LaFeAsO due to better nesting",
        prediction="η(SmFeAsO) ≈ 0.12, η(LaFeAsO) ≈ 0.33",
        measurement="Comparative optical or NMR study",
        expected_value="Ratio η_La/η_Sm ≈ 2.7",
        tolerance="±0.5",
        falsification="If ratio is near 1 or inverted",
        priority="High",
        estimated_cost="$30K-50K (existing samples)"
    ),
    TestableHypothesis(
        id="H300.4",
        hypothesis="FeSe/STO has higher η than bulk FeSe despite higher T_c",
        prediction="η(FeSe/STO) ≈ 0.85, η(FeSe bulk) ≈ 0.20",
        measurement="ARPES quasiparticle linewidth comparison",
        expected_value="Γ(T) slope 4× higher in monolayer",
        tolerance="±50%",
        falsification="If monolayer has lower thermal broadening",
        priority="High",
        estimated_cost="$50K-100K (synchrotron beamtime)"
    ),
    TestableHypothesis(
        id="H300.5",
        hypothesis="η correlates with disorder sensitivity",
        prediction="Lower η materials more robust to disorder",
        measurement="T_c suppression rate vs defect concentration",
        expected_value="dT_c/dn ∝ η² (pair-breaking formula)",
        tolerance="±30%",
        falsification="If correlation is absent or inverted",
        priority="Medium",
        estimated_cost="$20K-40K (ion irradiation + transport)"
    ),
    TestableHypothesis(
        id="H300.6",
        hypothesis="Cuprate/STO interface shows η reduction",
        prediction="YBCO/STO superlattice has η < bulk YBCO",
        measurement="Penetration depth λ(T) comparison",
        expected_value="η(interface) ≈ 0.30 vs η(bulk) ≈ 0.38",
        tolerance="±0.05",
        falsification="If interface η is equal or higher",
        priority="Medium",
        estimated_cost="$100K-200K (MBE + characterization)"
    ),
    TestableHypothesis(
        id="H300.7",
        hypothesis="Universal T_c × η relation",
        prediction="T_c × η × k_B / Δ ≈ 0.57 for all unconventional SC",
        measurement="Compile η, T_c, Δ for 10+ materials",
        expected_value="0.57 ± 0.15",
        tolerance="±25%",
        falsification="If scatter exceeds factor of 2",
        priority="High",
        estimated_cost="$10K-20K (literature compilation + analysis)"
    ),
]

print("\nTestable Hypothesis Matrix:")
print("-" * 100)
print(f"{'ID':<8} {'Hypothesis':<35} {'Priority':<10} {'Cost':<20} {'Falsifiable?':<10}")
print("-" * 100)
for h in hypotheses:
    hyp_short = h.hypothesis[:33] + "..." if len(h.hypothesis) > 35 else h.hypothesis
    print(f"{h.id:<8} {hyp_short:<35} {h.priority:<10} {h.estimated_cost:<20} {'Yes':<10}")

# ============================================================================
# PART 4: EXPERIMENTAL CAMPAIGN DESIGN
# ============================================================================

print("\n" + "=" * 60)
print("PART 4: EXPERIMENTAL CAMPAIGN DESIGN")
print("=" * 60)

@dataclass
class ExperimentalCampaign:
    """Coordinated experimental campaign to validate η"""
    name: str
    objective: str
    techniques: List[Technique]
    materials: List[str]
    duration: str
    total_cost: str
    expected_outcomes: List[str]
    risk_factors: List[str]

campaigns = [
    ExperimentalCampaign(
        name="Phase 1: Cuprate Benchmark",
        objective="Establish η values for well-characterized cuprates as reference",
        techniques=[Technique.NMR, Technique.OPTICAL, Technique.PENETRATION],
        materials=["YBCO", "Bi-2212", "LSCO", "Hg-1223"],
        duration="12 months",
        total_cost="$200K-400K",
        expected_outcomes=[
            "Complete η table for major cuprates",
            "Cross-technique validation",
            "Establish measurement precision",
            "Identify systematic effects"
        ],
        risk_factors=[
            "Sample quality variation",
            "Interpretation ambiguities",
            "Equipment downtime"
        ]
    ),
    ExperimentalCampaign(
        name="Phase 2: Pnictide Comparison",
        objective="Map η across iron pnictide families to test nesting hypothesis",
        techniques=[Technique.NMR, Technique.ARPES],
        materials=["LaFeAsO:F", "SmFeAsO:F", "BaFe₂As₂", "FeSe"],
        duration="12 months",
        total_cost="$300K-500K",
        expected_outcomes=[
            "η ordering: 1111 < 122 < 11",
            "Correlation with nesting quality",
            "k-resolved η(k) maps",
            "Test P298.1, P298.4"
        ],
        risk_factors=[
            "Crystal growth challenges",
            "Synchrotron beam time competition",
            "Doping control"
        ]
    ),
    ExperimentalCampaign(
        name="Phase 3: Interface Engineering",
        objective="Test if interfaces can reduce η and/or enhance Δ",
        techniques=[Technique.OPTICAL, Technique.TUNNELING, Technique.PENETRATION],
        materials=["YBCO/STO", "FeSe/BaTiO₃", "Cuprate/pnictide hybrid"],
        duration="18 months",
        total_cost="$500K-1M",
        expected_outcomes=[
            "Validate P299.1 (cuprate/STO enhancement)",
            "Validate P299.2 (FeSe/ferroelectric)",
            "Establish interface growth protocols",
            "Quantify disorder effects on η"
        ],
        risk_factors=[
            "Interface quality control",
            "High synthesis difficulty",
            "Reproducibility challenges"
        ]
    ),
    ExperimentalCampaign(
        name="Phase 4: Universal Scaling",
        objective="Test universal T_c × η × k_B / Δ relation",
        techniques=[Technique.NMR, Technique.OPTICAL, Technique.TUNNELING],
        materials=["20+ superconductors across all families"],
        duration="24 months",
        total_cost="$200K-400K (mostly data compilation)",
        expected_outcomes=[
            "Universal scaling confirmation or refutation",
            "Identify outliers and their physics",
            "Comprehensive η database",
            "Publication-ready results"
        ],
        risk_factors=[
            "Incomplete data for some materials",
            "Systematic differences between techniques",
            "Selection bias"
        ]
    ),
]

print("\nExperimental Campaign Design:")
print("-" * 90)
print(f"{'Phase':<30} {'Duration':<12} {'Cost':<20} {'Materials':<25}")
print("-" * 90)
for campaign in campaigns:
    materials_short = ", ".join(campaign.materials[:3])
    if len(campaign.materials) > 3:
        materials_short += f" +{len(campaign.materials)-3}"
    print(f"{campaign.name:<30} {campaign.duration:<12} {campaign.total_cost:<20} {materials_short:<25}")

# ============================================================================
# PART 5: SUCCESS/FAILURE CRITERIA
# ============================================================================

print("\n" + "=" * 60)
print("PART 5: SUCCESS/FAILURE CRITERIA")
print("=" * 60)

success_criteria = """
╔═══════════════════════════════════════════════════════════════════════════════╗
║                        SUCCESS / FAILURE CRITERIA                              ║
╠═══════════════════════════════════════════════════════════════════════════════╣
║                                                                                ║
║  FRAMEWORK VALIDATED IF:                                                       ║
║  ─────────────────────────                                                     ║
║  ✓ η values measurable with <20% precision across 3+ techniques               ║
║  ✓ η ordering matches T_c ordering within families (>80% cases)               ║
║  ✓ Universal scaling T_c × η / Δ within factor of 2                          ║
║  ✓ At least 5 of 7 hypotheses (H300.1-H300.7) pass                           ║
║  ✓ Predictions P298.1-P299.6 at least 60% confirmed                          ║
║                                                                                ║
║  FRAMEWORK REFUTED IF:                                                         ║
║  ─────────────────────────                                                     ║
║  ✗ η cannot be consistently extracted (>50% variation between techniques)     ║
║  ✗ η ordering ANTI-correlates with T_c ordering                              ║
║  ✗ T_c × η / Δ varies by >10× across materials                               ║
║  ✗ More than 4 of 7 hypotheses fail                                          ║
║  ✗ Predictions systematically wrong (>70% failure)                           ║
║                                                                                ║
║  FRAMEWORK NEEDS REVISION IF:                                                  ║
║  ───────────────────────────────                                               ║
║  ? η measurable but shows unexpected systematics                              ║
║  ? Partial agreement with predictions (40-60% success)                        ║
║  ? Some material classes fit, others don't                                    ║
║  ? Additional parameters needed beyond η                                      ║
║                                                                                ║
╚═══════════════════════════════════════════════════════════════════════════════╝
"""
print(success_criteria)

# ============================================================================
# PART 6: COST AND TIMELINE SUMMARY
# ============================================================================

print("\n" + "=" * 60)
print("PART 6: COST AND TIMELINE SUMMARY")
print("=" * 60)

# Calculate totals
total_duration = "5-6 years (parallel execution possible)"
min_cost = 200 + 300 + 500 + 200  # $K
max_cost = 400 + 500 + 1000 + 400

print(f"\nTotal Estimated Budget: ${min_cost}K - ${max_cost}K")
print(f"Total Timeline: {total_duration}")
print()

timeline = """
Year 1-2:  Phase 1 (Cuprate Benchmark) + Phase 2 start
Year 2-3:  Phase 2 (Pnictide Comparison) + Phase 3 start
Year 3-4:  Phase 3 (Interface Engineering)
Year 4-5:  Phase 4 (Universal Scaling) + Analysis
Year 5-6:  Publication and follow-up experiments

Critical Milestones:
  Month 6:   First η measurement validated across 2 techniques
  Month 12:  Cuprate η table complete
  Month 18:  Pnictide η comparison available
  Month 24:  Interface effects quantified
  Month 36:  Universal scaling tested
  Month 48:  Comprehensive results available
  Month 60:  Major publication submitted
"""
print(timeline)

# ============================================================================
# PART 7: POTENTIAL COLLABORATORS
# ============================================================================

print("\n" + "=" * 60)
print("PART 7: POTENTIAL COLLABORATORS")
print("=" * 60)

collaborators = """
INSTITUTIONS WITH RELEVANT EXPERTISE:

High-Temperature Superconductivity:
  • Brookhaven National Laboratory (ARPES, neutron)
  • Argonne National Laboratory (APS synchrotron)
  • Stanford/SLAC (ARPES, X-ray)
  • ETH Zurich (crystal growth, transport)
  • Tokyo University (high-quality crystals)

NMR Spectroscopy:
  • National High Magnetic Field Lab (high-field NMR)
  • Ames Laboratory (solid-state NMR)
  • McMaster University (heavy fermion NMR)

Optical Spectroscopy:
  • UC San Diego (Basov group - IR)
  • Rutgers University (Homes group)
  • Max Planck Stuttgart (optical studies)

Thin Film Growth:
  • Cornell (MBE of oxides)
  • Stanford (pulsed laser deposition)
  • NIST (interface characterization)

Penetration Depth:
  • University of Florida (microwave cavity)
  • University of Illinois (tunnel diode oscillator)
  • Cambridge (muon spin rotation)

SYNCHROTRON FACILITIES:
  • ALS (Berkeley) - ARPES
  • APS (Argonne) - X-ray
  • SSRL (Stanford) - ARPES
  • Diamond (UK) - ARPES
  • SLS (Switzerland) - ARPES
"""
print(collaborators)

# ============================================================================
# PART 8: GENERATE VISUALIZATIONS
# ============================================================================

print("\n" + "=" * 60)
print("PART 8: GENERATING VISUALIZATIONS")
print("=" * 60)

fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle('Session #300: Experimental Validation Protocol', fontsize=16, fontweight='bold')

# Plot 1: Technique comparison
ax1 = axes[0, 0]
tech_names = [t.name.value[:15] for t in techniques]
precisions = [t.precision for t in techniques]
colors = plt.cm.RdYlGn_r(np.array(precisions) / max(precisions))
bars = ax1.barh(range(len(tech_names)), precisions, color=colors)
ax1.set_yticks(range(len(tech_names)))
ax1.set_yticklabels(tech_names, fontsize=9)
ax1.set_xlabel('Precision (% error)', fontsize=12)
ax1.set_title('Technique Precision Comparison', fontsize=12)
ax1.axvline(x=5, color='green', linestyle='--', label='Target: 5%')
ax1.legend()

# Plot 2: Hypothesis priority
ax2 = axes[0, 1]
priority_map = {"High": 3, "Medium": 2, "Low": 1}
hyp_ids = [h.id for h in hypotheses]
priorities = [priority_map[h.priority] for h in hypotheses]
colors = ['red' if p == 3 else 'orange' if p == 2 else 'yellow' for p in priorities]
ax2.bar(hyp_ids, priorities, color=colors, alpha=0.7)
ax2.set_ylabel('Priority Level', fontsize=12)
ax2.set_title('Hypothesis Priority', fontsize=12)
ax2.set_yticks([1, 2, 3])
ax2.set_yticklabels(['Low', 'Medium', 'High'])
ax2.tick_params(axis='x', rotation=45)

# Plot 3: Campaign timeline
ax3 = axes[0, 2]
campaign_names = [c.name.split(":")[0] for c in campaigns]
durations = [int(c.duration.split()[0]) for c in campaigns]
starts = [0, 0, 12, 24]
colors = plt.cm.tab10(np.arange(len(campaigns)))
for i, (name, dur, start) in enumerate(zip(campaign_names, durations, starts)):
    ax3.barh(i, dur, left=start, color=colors[i], alpha=0.7)
    ax3.text(start + dur/2, i, f"{dur}mo", ha='center', va='center', fontsize=9)
ax3.set_yticks(range(len(campaign_names)))
ax3.set_yticklabels(campaign_names, fontsize=9)
ax3.set_xlabel('Months', fontsize=12)
ax3.set_title('Campaign Timeline', fontsize=12)
ax3.axvline(x=12, color='gray', linestyle=':', alpha=0.5)
ax3.axvline(x=24, color='gray', linestyle=':', alpha=0.5)
ax3.axvline(x=36, color='gray', linestyle=':', alpha=0.5)
ax3.axvline(x=48, color='gray', linestyle=':', alpha=0.5)

# Plot 4: Cost breakdown
ax4 = axes[1, 0]
costs_min = [200, 300, 500, 200]
costs_max = [400, 500, 1000, 400]
campaign_short = [c.name.split(":")[0] for c in campaigns]
x = np.arange(len(campaigns))
width = 0.35
ax4.bar(x - width/2, costs_min, width, label='Min cost', color='green', alpha=0.7)
ax4.bar(x + width/2, costs_max, width, label='Max cost', color='red', alpha=0.7)
ax4.set_xticks(x)
ax4.set_xticklabels(campaign_short, rotation=45, ha='right', fontsize=9)
ax4.set_ylabel('Cost ($K)', fontsize=12)
ax4.set_title('Campaign Cost Estimates', fontsize=12)
ax4.legend()

# Plot 5: Materials tested
ax5 = axes[1, 1]
all_materials = set()
for c in campaigns:
    all_materials.update(c.materials)
material_counts = {}
for mat in all_materials:
    count = sum(1 for c in campaigns if mat in c.materials)
    material_counts[mat] = count
sorted_mats = sorted(material_counts.items(), key=lambda x: -x[1])[:10]
mat_names = [m[0][:15] for m in sorted_mats]
mat_counts = [m[1] for m in sorted_mats]
ax5.barh(range(len(mat_names)), mat_counts, color='steelblue', alpha=0.7)
ax5.set_yticks(range(len(mat_names)))
ax5.set_yticklabels(mat_names, fontsize=9)
ax5.set_xlabel('Number of Campaigns', fontsize=12)
ax5.set_title('Materials Across Campaigns', fontsize=12)

# Plot 6: Techniques per campaign
ax6 = axes[1, 2]
tech_usage = {t: 0 for t in Technique}
for c in campaigns:
    for t in c.techniques:
        tech_usage[t] += 1
tech_names_plot = [t.value[:15] for t in tech_usage.keys()]
tech_counts = list(tech_usage.values())
ax6.bar(tech_names_plot, tech_counts, color='purple', alpha=0.7)
ax6.set_ylabel('Number of Campaigns', fontsize=12)
ax6.set_title('Technique Usage', fontsize=12)
ax6.tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('session300_experimental_validation_protocol.png', dpi=150, bbox_inches='tight')
print("\nVisualization saved: session300_experimental_validation_protocol.png")

# ============================================================================
# PART 9: ARC STATUS AND SUMMARY
# ============================================================================

print("\n" + "=" * 60)
print("PART 9: HOT SUPERCONDUCTOR ARC STATUS")
print("=" * 60)

arc_summary = """
╔═══════════════════════════════════════════════════════════════════════════════╗
║                     HOT SUPERCONDUCTOR ARC SUMMARY                            ║
╠═══════════════════════════════════════════════════════════════════════════════╣
║                                                                                ║
║  SESSIONS COMPLETED:                                                          ║
║  ───────────────────                                                           ║
║  #292: Dissonance Pathway Formalization                                       ║
║        → η framework introduced, T_c = Δ / (1.76 k_B × η)                    ║
║                                                                                ║
║  #297: Cuprate η Quantification                                               ║
║        → YBCO: 0.38, Bi-2212: 0.42, LSCO: 0.51, Hg-1223: 0.33               ║
║                                                                                ║
║  #298: Iron Pnictide η Analysis                                               ║
║        → SmFeAsO: 0.12 (lowest!), FeSe/STO: 0.85 (highest)                   ║
║                                                                                ║
║  #299: Minimum-η Material Design                                              ║
║        → 8 material stacks proposed, η-Δ trade-off identified                 ║
║                                                                                ║
║  #300: Experimental Validation Protocol (THIS SESSION)                        ║
║        → Measurement techniques, testable hypotheses, campaign design         ║
║                                                                                ║
║  KEY INSIGHTS:                                                                 ║
║  ─────────────────                                                             ║
║  • η framework explains T_c ordering across materials                         ║
║  • Best η: SmFeAsO ≈ 0.12 (perfect nesting)                                  ║
║  • Cuprates: η ≈ 0.33-0.51 (d-wave + correlations)                           ║
║  • Room-temp SC needs: η < 0.15 + Δ > 10 meV OR η ~ 0.5 + Δ > 30 meV        ║
║  • Interface engineering enhances Δ, not η                                    ║
║                                                                                ║
║  PREDICTIONS GENERATED: 30+                                                    ║
║  TESTABLE HYPOTHESES: 7 major (H300.1-H300.7)                                 ║
║  ESTIMATED VALIDATION COST: $1.2M - $2.3M over 5-6 years                     ║
║                                                                                ║
║  NEXT STEPS:                                                                   ║
║  ───────────────                                                               ║
║  • Seek experimental collaborators                                            ║
║  • Begin Phase 1 (Cuprate Benchmark)                                          ║
║  • Compile existing literature for preliminary η values                       ║
║  • Consider arc completion or new direction                                   ║
║                                                                                ║
╚═══════════════════════════════════════════════════════════════════════════════╝
"""
print(arc_summary)

# ============================================================================
# SUMMARY
# ============================================================================

print("\n" + "=" * 80)
print("SESSION #300 COMPLETE - MILESTONE SESSION")
print("HOT SUPERCONDUCTOR ARC (Session 5/?)")
print("=" * 80)

print("""
Key Achievements:
  • 6 measurement techniques analyzed for η extraction
  • 4 detailed extraction protocols defined
  • 7 testable hypotheses with falsification criteria
  • 4-phase experimental campaign designed
  • Success/failure criteria established
  • $1.2M-$2.3M budget over 5-6 years estimated

Critical Path:
  Phase 1: Cuprate benchmark (12 months, $200-400K)
  Phase 2: Pnictide comparison (12 months, $300-500K)
  Phase 3: Interface engineering (18 months, $500K-1M)
  Phase 4: Universal scaling (24 months, $200-400K)

The η framework is now experimentally testable with clear metrics
for validation, refutation, or revision.

NEXT: Arc completion assessment or new research direction
""")
