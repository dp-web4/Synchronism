#!/usr/bin/env python3
"""
Chemistry Session #28: Validation Priority Matrix

Ranks all framework predictions by:
1. Impact (if validated, how much does it advance the framework?)
2. Feasibility (can it be tested with available methods?)
3. Falsification Power (does failure kill the framework?)

Goal: Identify the 10 highest-priority experiments to propose.
"""

import numpy as np
import matplotlib.pyplot as plt

def create_prediction_database():
    """
    Compile all predictions with their scores.

    Impact: 1-5 (1=incremental, 5=paradigm-shifting)
    Feasibility: 1-5 (1=requires new tech, 5=straightforward)
    Falsification: 1-5 (1=doesn't affect framework, 5=kills framework if fails)
    """
    predictions = {
        # Category 1: Superconductivity
        'P1.1': {'name': 'BCS Gap Ratio', 'impact': 5, 'feasibility': 5,
                 'falsification': 5, 'status': 'VALIDATED', 'category': 'SC'},
        'P1.2': {'name': 'Cuprate Gap Ratios', 'impact': 4, 'feasibility': 4,
                 'falsification': 4, 'status': 'testable', 'category': 'SC'},
        'P1.3': {'name': 'Doping Dome', 'impact': 3, 'feasibility': 3,
                 'falsification': 3, 'status': 'testable', 'category': 'SC'},
        'P1.5': {'name': 'Pressure Effects', 'impact': 4, 'feasibility': 3,
                 'falsification': 4, 'status': 'testable', 'category': 'SC'},
        'P1.6': {'name': 'Disorder Effects', 'impact': 3, 'feasibility': 4,
                 'falsification': 3, 'status': 'testable', 'category': 'SC'},
        'P1.7': {'name': 'Hydride Gap Ratios', 'impact': 4, 'feasibility': 2,
                 'falsification': 4, 'status': 'testable', 'category': 'SC'},
        'P1.8': {'name': 'Hydride Tc Predictions', 'impact': 5, 'feasibility': 2,
                 'falsification': 5, 'status': 'testable', 'category': 'SC'},

        # Category 2: Catalysis
        'P2.4': {'name': 'KIE-γ Correlation', 'impact': 5, 'feasibility': 5,
                 'falsification': 5, 'status': 'VALIDATED', 'category': 'Catalysis'},
        'P2.5': {'name': 'High-KIE Enzymes', 'impact': 4, 'feasibility': 3,
                 'falsification': 4, 'status': 'testable', 'category': 'Catalysis'},
        'P2.6': {'name': 'Enzyme Mutations', 'impact': 4, 'feasibility': 4,
                 'falsification': 4, 'status': 'testable', 'category': 'Catalysis'},
        'P2.7': {'name': 'KIE Temperature Dependence', 'impact': 3, 'feasibility': 4,
                 'falsification': 3, 'status': 'testable', 'category': 'Catalysis'},

        # Category 3: Bonding
        'P3.2': {'name': "Hückel's Rule", 'impact': 4, 'feasibility': 5,
                 'falsification': 4, 'status': 'VALIDATED', 'category': 'Bonding'},
        'P3.3': {'name': 'Lone Pair Interference', 'impact': 3, 'feasibility': 3,
                 'falsification': 3, 'status': 'testable', 'category': 'Bonding'},
        'P3.4': {'name': 'Bond Angles', 'impact': 2, 'feasibility': 5,
                 'falsification': 2, 'status': 'partial', 'category': 'Bonding'},

        # Category 4: Phase Transitions
        'P4.1': {'name': 'Glass Fragility', 'impact': 3, 'feasibility': 4,
                 'falsification': 3, 'status': 'testable', 'category': 'Phase'},
        'P4.2': {'name': 'Melting Point Model', 'impact': 2, 'feasibility': 5,
                 'falsification': 2, 'status': 'FAILED', 'category': 'Phase'},

        # Category 5: Photosynthesis
        'P5.1': {'name': 'Efficiency-γ Correlation', 'impact': 4, 'feasibility': 3,
                 'falsification': 4, 'status': 'testable', 'category': 'Photo'},
        'P5.2': {'name': 'Protein Mutation Effects', 'impact': 4, 'feasibility': 3,
                 'falsification': 4, 'status': 'testable', 'category': 'Photo'},
        'P5.4': {'name': 'Artificial Light Harvesting', 'impact': 5, 'feasibility': 2,
                 'falsification': 4, 'status': 'testable', 'category': 'Photo'},

        # Category 6: Cross-Domain
        'P6.1': {'name': 'Universal γ Reduction', 'impact': 5, 'feasibility': 4,
                 'falsification': 5, 'status': 'testable', 'category': 'Universal'},
        'P6.2': {'name': 'N_corr Mechanism', 'impact': 5, 'feasibility': 3,
                 'falsification': 5, 'status': 'testable', 'category': 'Universal'},

        # Category 9: Universal Synthesis
        'P9.1': {'name': 'Universal γ Bound', 'impact': 4, 'feasibility': 3,
                 'falsification': 5, 'status': 'testable', 'category': 'Universal'},
        'P9.3': {'name': 'Universal Temperature Scaling', 'impact': 5, 'feasibility': 4,
                 'falsification': 5, 'status': 'testable', 'category': 'Universal'},

        # Category 10: Quantum Computing
        'P10.1': {'name': 'Error Correction Scaling', 'impact': 4, 'feasibility': 3,
                  'falsification': 4, 'status': 'testable', 'category': 'QC'},
        'P10.4': {'name': 'Topological Size Limit', 'impact': 4, 'feasibility': 2,
                  'falsification': 4, 'status': 'testable', 'category': 'QC'},

        # Category 11: Magnetism
        'P11.1': {'name': 'Critical Exponent Relation', 'impact': 5, 'feasibility': 4,
                  'falsification': 5, 'status': 'testable', 'category': 'Magnetism'},
        'P11.4': {'name': 'Cuprate-AF Connection', 'impact': 4, 'feasibility': 3,
                  'falsification': 4, 'status': 'testable', 'category': 'Magnetism'},

        # Category 12: Thermodynamics
        'P12.2': {'name': 'Entropy Reduction', 'impact': 4, 'feasibility': 4,
                  'falsification': 4, 'status': 'testable', 'category': 'Thermo'},

        # New from Sessions 26-27
        'P26.1': {'name': 'N_corr from Coherence Length', 'impact': 4, 'feasibility': 4,
                  'falsification': 4, 'status': 'testable', 'category': 'Measurement'},
        'P26.2': {'name': 'Market Crash N_corr Signature', 'impact': 3, 'feasibility': 4,
                  'falsification': 3, 'status': 'testable', 'category': 'Economics'},
        'P26.3': {'name': 'Consciousness N_corr', 'impact': 4, 'feasibility': 3,
                  'falsification': 3, 'status': 'testable', 'category': 'Neuro'},
        'P27.1': {'name': 'α from Mechanism', 'impact': 5, 'feasibility': 4,
                  'falsification': 5, 'status': 'testable', 'category': 'Catalysis'},
        'P27.2': {'name': 'Multi-Proton α > 1.5', 'impact': 4, 'feasibility': 4,
                  'falsification': 4, 'status': 'testable', 'category': 'Catalysis'},
        'P27.3': {'name': 'α-KIE Correlation', 'impact': 4, 'feasibility': 4,
                  'falsification': 4, 'status': 'testable', 'category': 'Catalysis'},
    }

    return predictions

def calculate_priority_score(pred):
    """
    Calculate composite priority score.

    Score = Impact × Feasibility × Falsification^0.5
    (Falsification weighted less to avoid avoiding risky tests)
    """
    if pred['status'] == 'VALIDATED' or pred['status'] == 'FAILED':
        return 0  # Already resolved
    return pred['impact'] * pred['feasibility'] * np.sqrt(pred['falsification'])

def main():
    """Create validation priority analysis."""

    predictions = create_prediction_database()

    # Calculate scores
    for pid, pred in predictions.items():
        pred['score'] = calculate_priority_score(pred)

    # Sort by score
    sorted_preds = sorted(predictions.items(), key=lambda x: x[1]['score'], reverse=True)

    fig, axes = plt.subplots(2, 3, figsize=(14, 10))
    fig.suptitle('Chemistry Session #28: Validation Priority Matrix', fontsize=14, fontweight='bold')

    # Part 1: Top 15 predictions by priority score
    ax1 = axes[0, 0]

    top_n = 15
    top_preds = [p for p in sorted_preds if p[1]['status'] == 'testable'][:top_n]
    names = [f"{p[0]}: {p[1]['name'][:20]}" for p in top_preds]
    scores = [p[1]['score'] for p in top_preds]

    colors = plt.cm.viridis(np.linspace(0.3, 0.9, top_n))
    bars = ax1.barh(range(top_n), scores, color=colors, edgecolor='black')
    ax1.set_yticks(range(top_n))
    ax1.set_yticklabels(names, fontsize=8)
    ax1.set_xlabel('Priority Score')
    ax1.set_title('Top 15 Predictions to Test')
    ax1.invert_yaxis()

    # Part 2: Impact vs Feasibility scatter
    ax2 = axes[0, 1]

    testable = {k: v for k, v in predictions.items() if v['status'] == 'testable'}

    impacts = [v['impact'] for v in testable.values()]
    feasibilities = [v['feasibility'] for v in testable.values()]
    falsifications = [v['falsification'] for v in testable.values()]

    scatter = ax2.scatter(feasibilities, impacts, s=[f*30 for f in falsifications],
                         c=falsifications, cmap='RdYlGn', edgecolor='black', alpha=0.7)

    # Add labels for high-priority
    for pid, pred in testable.items():
        if pred['score'] > 40:
            ax2.annotate(pid, (pred['feasibility'], pred['impact']),
                        fontsize=7, alpha=0.8)

    ax2.set_xlabel('Feasibility')
    ax2.set_ylabel('Impact')
    ax2.set_title('Impact vs Feasibility\n(size/color = falsification power)')
    ax2.set_xlim(0.5, 5.5)
    ax2.set_ylim(0.5, 5.5)

    # Add quadrant labels
    ax2.axhline(y=3, color='gray', linestyle='--', alpha=0.3)
    ax2.axvline(x=3, color='gray', linestyle='--', alpha=0.3)
    ax2.text(4.5, 4.5, 'HIGH\nPRIORITY', ha='center', va='center', fontsize=10, color='green')
    ax2.text(1.5, 4.5, 'HIGH IMPACT\nLOW FEASIBLE', ha='center', va='center', fontsize=8, color='orange')
    ax2.text(4.5, 1.5, 'INCREMENTAL', ha='center', va='center', fontsize=8, color='gray')

    # Part 3: By category
    ax3 = axes[0, 2]

    categories = {}
    for pid, pred in predictions.items():
        cat = pred['category']
        if cat not in categories:
            categories[cat] = {'testable': 0, 'validated': 0, 'failed': 0}
        if pred['status'] == 'testable':
            categories[cat]['testable'] += 1
        elif pred['status'] == 'VALIDATED':
            categories[cat]['validated'] += 1
        elif pred['status'] == 'FAILED':
            categories[cat]['failed'] += 1

    cat_names = list(categories.keys())
    x = np.arange(len(cat_names))
    width = 0.25

    testable_counts = [categories[c]['testable'] for c in cat_names]
    validated_counts = [categories[c]['validated'] for c in cat_names]
    failed_counts = [categories[c]['failed'] for c in cat_names]

    ax3.bar(x - width, testable_counts, width, label='Testable', color='steelblue')
    ax3.bar(x, validated_counts, width, label='Validated', color='green')
    ax3.bar(x + width, failed_counts, width, label='Failed', color='red')

    ax3.set_xticks(x)
    ax3.set_xticklabels(cat_names, rotation=45, ha='right', fontsize=8)
    ax3.set_ylabel('Count')
    ax3.set_title('Predictions by Category')
    ax3.legend(fontsize=8)

    # Part 4: Top 10 detailed
    ax4 = axes[1, 0]

    ax4.axis('off')

    top10_text = "TOP 10 HIGHEST-PRIORITY EXPERIMENTS\n"
    top10_text += "=" * 50 + "\n\n"

    for i, (pid, pred) in enumerate(top_preds[:10], 1):
        top10_text += f"{i}. {pid}: {pred['name']}\n"
        top10_text += f"   Impact: {pred['impact']} | Feasibility: {pred['feasibility']} | "
        top10_text += f"Falsification: {pred['falsification']}\n"
        top10_text += f"   Category: {pred['category']} | Score: {pred['score']:.1f}\n\n"

    ax4.text(0.05, 0.95, top10_text, transform=ax4.transAxes, fontsize=8,
             verticalalignment='top', family='monospace')

    # Part 5: Framework-critical predictions
    ax5 = axes[1, 1]

    ax5.axis('off')

    critical_text = "FRAMEWORK-CRITICAL PREDICTIONS\n"
    critical_text += "(Falsification Power = 5)\n"
    critical_text += "=" * 40 + "\n\n"

    critical = [p for p in sorted_preds if p[1]['falsification'] == 5 and p[1]['status'] == 'testable']

    for pid, pred in critical[:8]:
        critical_text += f"• {pid}: {pred['name']}\n"
        critical_text += f"  If FAILS: Major framework revision needed\n"
        critical_text += f"  If PASSES: Strong confirmation\n\n"

    ax5.text(0.05, 0.95, critical_text, transform=ax5.transAxes, fontsize=9,
             verticalalignment='top', family='monospace')

    # Part 6: Proposed validation roadmap
    ax6 = axes[1, 2]

    ax6.axis('off')

    roadmap = """
PROPOSED VALIDATION ROADMAP

PHASE 1: Immediate (Existing Data)
--------------------------------
1. P27.1: α from mechanism
   - Use published enzyme data
   - Calculate α from mechanism

2. P9.3: Universal Tc scaling
   - Compare Tc/(2/γ) ratios
   - Across SC, magnets, etc.

3. P11.1: Critical exponent β = 1/(2γ)
   - Check literature values
   - Multiple magnetic materials

PHASE 2: Near-Term (Feasible Now)
---------------------------------
4. P2.6: Enzyme mutations
   - Express mutants
   - Measure KIE changes

5. P6.1: Universal γ reduction
   - Survey known coherent systems
   - Verify γ < γ_standard

6. P12.2: Entropy reduction
   - Calorimetry of correlated systems
   - Compare S to prediction

PHASE 3: Medium-Term (New Experiments)
--------------------------------------
7. P1.5: Pressure effects on cuprates
   - High-pressure gap measurements
   - Test γ(P) predictions

8. P5.2: Photosynthesis mutations
   - Express LH complex mutants
   - Measure coherence time

PHASE 4: Long-Term (New Materials)
----------------------------------
9. P1.8: Hydride Tc predictions
   - Synthesize MgH₆, BeH₈
   - Measure Tc

10. P5.4: Artificial light harvesting
    - Design γ < 0.5 chromophore arrays
    - Measure efficiency
"""
    ax6.text(0.05, 0.95, roadmap, transform=ax6.transAxes, fontsize=8,
             verticalalignment='top', family='monospace')

    plt.tight_layout()
    plt.savefig('/mnt/c/exe/projects/ai-agents/Synchronism/simulations/chemistry/validation_priority.png',
                dpi=150, bbox_inches='tight')
    plt.close()

    # Print summary
    print("=" * 70)
    print("Chemistry Session #28: Validation Priority Matrix")
    print("=" * 70)
    print()
    print("FRAMEWORK STATUS")
    print("-" * 70)

    total = len(predictions)
    validated = sum(1 for p in predictions.values() if p['status'] == 'VALIDATED')
    failed = sum(1 for p in predictions.values() if p['status'] == 'FAILED')
    testable = sum(1 for p in predictions.values() if p['status'] == 'testable')

    print(f"Total Predictions: {total}")
    print(f"  Validated: {validated} ({100*validated/total:.1f}%)")
    print(f"  Failed: {failed} ({100*failed/total:.1f}%)")
    print(f"  Testable: {testable} ({100*testable/total:.1f}%)")
    print()
    print("-" * 70)
    print("TOP 10 HIGHEST-PRIORITY EXPERIMENTS")
    print("-" * 70)
    print()

    for i, (pid, pred) in enumerate(top_preds[:10], 1):
        print(f"{i:2}. {pid}: {pred['name']}")
        print(f"    Score: {pred['score']:.1f} | Impact: {pred['impact']} | "
              f"Feasibility: {pred['feasibility']} | Falsification: {pred['falsification']}")
        print()

    print("-" * 70)
    print("FRAMEWORK-CRITICAL PREDICTIONS (Falsification = 5)")
    print("-" * 70)
    print()

    for pid, pred in critical[:6]:
        print(f"• {pid}: {pred['name']}")
        print(f"  Category: {pred['category']} | Score: {pred['score']:.1f}")
        print()

    print("-" * 70)
    print("RECOMMENDED NEXT STEPS")
    print("-" * 70)
    print()
    print("1. Validate P27.1 (α from mechanism) using existing enzyme literature")
    print("2. Check P9.3 (universal Tc scaling) across material classes")
    print("3. Test P11.1 (β = 1/2γ) using published critical exponents")
    print()
    print("These three can be done with EXISTING DATA - no new experiments needed.")
    print()
    print("=" * 70)
    print("PRIORITY MATRIX ESTABLISHED")
    print("=" * 70)

if __name__ == "__main__":
    main()
