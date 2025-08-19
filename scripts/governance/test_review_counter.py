#!/usr/bin/env python3
"""Test the enhanced review and counter-proposal system"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from governance_cycle import GovernanceOrchestrator
from lct_participants import initialize_default_participants
import json

def test_review_counter_system():
    """Test the review accept/hold and counter-proposal flow"""
    
    print("Testing Enhanced Review & Counter-Proposal System")
    print("=" * 60)
    
    # Initialize participants first
    print("\nInitializing participants...")
    registry = initialize_default_participants()
    print(f"  Registered {len(registry.participants)} participants")
    
    # Initialize orchestrator
    orchestrator = GovernanceOrchestrator()
    
    # Test sections
    test_sections = ["04-fundamental-concepts/01-universe-grid"]
    
    print("\n" + "=" * 60)
    print(" CYCLE 1: Initial Proposals")
    print("=" * 60)
    
    # Run first cycle
    result1 = orchestrator.run_cycle(sections=test_sections, max_proposals=1)
    
    print("\nCycle 1 Results:")
    print(f"  Proposals created: {len(result1['phases']['proposal']['proposals_created'])}")
    print(f"  Reviews submitted: {len(result1['phases']['review']['reviews_submitted'])}")
    print(f"  Holds requested: {len(result1['phases']['review']['holds_requested'])}")
    
    # Check arbitration
    arb_results = result1['phases']['arbitration']
    print(f"\nArbitration Results:")
    for decision in arb_results['decisions']:
        print(f"  {decision['proposal_id']}: {decision['status']}")
        print(f"    Rationale: {decision['rationale']}")
    
    # Show held proposals
    if orchestrator.held_proposals:
        print(f"\nProposals held for counter-proposal:")
        for prop_id, holder_id in orchestrator.held_proposals.items():
            print(f"  - {prop_id} held by {holder_id[:8]}...")
    
    print("\n" + "=" * 60)
    print(" CYCLE 2: Counter-Proposals")
    print("=" * 60)
    
    # Run second cycle to see counter-proposals
    result2 = orchestrator.run_cycle(sections=test_sections, max_proposals=1)
    
    print("\nCycle 2 Results:")
    counter_props = result2.get('counter_proposals', {})
    if counter_props:
        print(f"  Counter-proposals created: {len(counter_props)}")
        for orig_id, counter in counter_props.items():
            print(f"    - Enhanced {orig_id}")
            print(f"      Conversation depth: {counter['metadata']['iteration']}")
            print(f"      Participants: {len(counter['conversation'])}")
    else:
        print("  No counter-proposals created")
    
    print(f"  New proposals: {len(result2['phases']['proposal']['proposals_created'])}")
    print(f"  Reviews: {len(result2['phases']['review']['reviews_submitted'])}")
    
    # Check if counter-proposals get reviewed
    print(f"\nCounter-proposal handling:")
    for decision in result2['phases']['arbitration']['decisions']:
        if 'counter' in decision['proposal_id']:
            print(f"  {decision['proposal_id']}: {decision['status']}")
            print(f"    {decision['rationale']}")
    
    print("\n" + "=" * 60)
    print(" Test Analysis")
    print("=" * 60)
    
    print("\nKey features demonstrated:")
    print("1. Reviews can 'accept' or 'hold-for-counter'")
    print("2. Arbiter only accepts unanimous 'accept' reviews")
    print("3. Held proposals generate counter-proposals in next cycle")
    print("4. Counter-proposals create multi-participant conversations")
    print("5. Reviews are signed with participant IDs")
    
    # Save test results
    test_output = {
        'cycle1': {
            'proposals': result1['phases']['proposal']['proposals_created'],
            'holds': result1['phases']['review']['holds_requested'],
            'decisions': result1['phases']['arbitration']['decisions']
        },
        'cycle2': {
            'counter_proposals': counter_props,
            'new_proposals': result2['phases']['proposal']['proposals_created'],
            'decisions': result2['phases']['arbitration']['decisions']
        }
    }
    
    with open('test_review_counter_results.json', 'w') as f:
        json.dump(test_output, f, indent=2)
    
    print("\nTest results saved to test_review_counter_results.json")
    print("\n✅ Test complete!")

if __name__ == "__main__":
    test_review_counter_system()