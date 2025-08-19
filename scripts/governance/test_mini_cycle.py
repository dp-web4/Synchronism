#!/usr/bin/env python3
"""Test a mini governance cycle with arbiter selection"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from governance_cycle import GovernanceOrchestrator

def test_mini_cycle():
    """Test a single governance cycle with limited proposals"""
    
    print("Testing mini governance cycle...")
    print("=" * 50)
    
    # Initialize orchestrator
    orchestrator = GovernanceOrchestrator()
    
    # Get initial state
    print("\nInitial state:")
    print(f"  Current cycle ID: {orchestrator.current_cycle_id}")
    history = orchestrator.load_history()
    if history:
        print(f"  Total cycles run: {len(history.get('cycles', []))}")
    else:
        print(f"  Total cycles run: 0 (no history file)")
    
    # Define test sections (just a couple to keep it small)
    test_sections = [
        "04-fundamental-concepts/01-universe-grid",
        "05-quantum-macro/01-crt-analogy"
    ]
    
    print(f"\nTest sections: {test_sections}")
    print(f"Max proposals per participant: 1 (configured)")
    
    print("\n" + "=" * 50)
    print("Running test cycle...")
    print("=" * 50)
    
    try:
        # Run a single cycle with our test sections
        result = orchestrator.run_cycle(
            sections=test_sections,
            max_proposals=1  # Already defaulted to 1
        )
        
        print("\nCycle completed successfully!")
        print(f"  Cycle ID: {result.get('cycle_id')}")
        
        # Extract nested phase data
        phases = result.get('phases', {})
        
        # Arbiter info
        arb_phase = phases.get('arbitration', {})
        arbiter = arb_phase.get('arbiter', {})
        print(f"  Arbiter: {arbiter.get('name', 'Unknown')} ({arbiter.get('lct_id', 'N/A')[:8]}...)")
        
        # Show proposals generated
        prop_phase = phases.get('proposal', {})
        proposals = prop_phase.get('proposals_created', [])
        print(f"\nProposals generated: {len(proposals)}")
        for prop_id in proposals[:3]:
            if isinstance(prop_id, str):
                print(f"  - {prop_id[:20]}...")
            else:
                print(f"  - {str(prop_id)[:50]}...")
        
        # Show reviews
        review_phase = phases.get('review', {})
        reviews = review_phase.get('reviews_submitted', [])
        print(f"\nReviews submitted: {len(reviews)}")
        
        # Show arbiter decisions
        decisions = arb_phase.get('decisions', [])
        print(f"\nArbiter decisions: {len(decisions)}")
        for dec_id in decisions[:3]:
            if isinstance(dec_id, str):
                print(f"  - {dec_id[:20]}...")
            else:
                print(f"  - {str(dec_id)[:50]}...")
        
        # Summary
        summary = result.get('summary', {})
        if summary:
            print(f"\nSummary:")
            print(f"  Accepted: {summary.get('accepted', 0)}")
            print(f"  Rejected: {summary.get('rejected', 0)}")
            print(f"  Deferred: {summary.get('deferred', 0)}")
        
    except Exception as e:
        print(f"\nError during cycle: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 50)
    print("Test complete!")

if __name__ == "__main__":
    test_mini_cycle()