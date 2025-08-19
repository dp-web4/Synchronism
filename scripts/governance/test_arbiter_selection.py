#!/usr/bin/env python3
"""Test arbiter selection logic"""

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from lct_participants import LCTRegistry, initialize_default_participants

def test_arbiter_selection():
    """Test the arbiter selection with unavailable human"""
    
    print("Testing arbiter selection logic...")
    print("=" * 50)
    
    # Initialize registry with default participants
    registry = initialize_default_participants()
    
    print(f"\nTotal participants: {len(registry.participants)}")
    for pid, lct in registry.participants.items():
        avail = lct.check_availability()
        print(f"  {lct.name} ({lct.type.value}): Available={avail}")
    
    print("\n" + "=" * 50)
    print("Testing arbiter selection scenarios:")
    print("=" * 50)
    
    # Scenario 1: Try to select human as arbiter (should fail and fallback)
    print("\n1. Trying to select human arbiter (Dennis)...")
    human_id = None
    for pid, lct in registry.participants.items():
        if lct.name == "Dennis":
            human_id = pid
            break
    
    if human_id:
        arbiter = registry.select_arbiter_with_fallback(preferred_id=human_id)
        if arbiter:
            print(f"   Selected: {arbiter.name} ({arbiter.type.value})")
        else:
            print("   ERROR: No arbiter selected!")
    
    # Scenario 2: Default arbiter selection
    print("\n2. Default arbiter selection (no preference)...")
    arbiter = registry.select_arbiter_with_fallback()
    if arbiter:
        print(f"   Selected: {arbiter.name} ({arbiter.type.value})")
    else:
        print("   ERROR: No arbiter selected!")
    
    # Scenario 3: Check available arbiters
    print("\n3. Checking all available arbiters...")
    available = registry.get_available_arbiters()
    print(f"   Found {len(available)} available arbiters:")
    for arb in available:
        print(f"     - {arb.name} (trust: {arb.trust_scores.aggregate():.2f})")
    
    print("\n" + "=" * 50)
    print("Test complete!")

if __name__ == "__main__":
    test_arbiter_selection()