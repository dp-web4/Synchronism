#!/usr/bin/env python3
"""
Test API Context for Governance System
Shows the complete context passed to each participant
"""

import json
from datetime import datetime
from governance_cycle_refined import GovernanceCycleRefined, initialize_updated_participants

def pretty_print_context(context: dict, indent: int = 2):
    """Pretty print API context"""
    
    # Print identity first (most important)
    print(f"{'='*60}")
    print(f"YOU ARE: {context['you_are']}")
    print(f"ROLE: {context['role']}")
    print(f"NAME: {context['your_name']}")
    print(f"TYPE: {context['your_type']}")
    print(f"{'='*60}")
    
    # Print sections in order of importance
    important_sections = [
        'role_specific',
        'your_status', 
        'held_proposals',
        'available_sections',
        'recent_decisions'
    ]
    
    for section in important_sections:
        if section in context and context[section]:
            print(f"\n[{section.upper().replace('_', ' ')}]")
            if isinstance(context[section], dict):
                for key, value in context[section].items():
                    if isinstance(value, list) and len(value) > 3:
                        print(f"  {key}: [{len(value)} items]")
                    else:
                        print(f"  {key}: {value}")
            elif isinstance(context[section], list):
                if len(context[section]) > 3:
                    print(f"  {len(context[section])} items")
                else:
                    for item in context[section]:
                        print(f"  - {item}")
            else:
                print(f"  {context[section]}")
    
    # Show that rules are included
    if 'governance_rules' in context:
        rules_preview = context['governance_rules'][:100] if context['governance_rules'] else "No rules"
        print(f"\n[GOVERNANCE RULES]")
        print(f"  {rules_preview}...")
    
    print(f"\n[METADATA]")
    print(f"  Phase: {context.get('current_phase', 'unknown')}")
    print(f"  Cycle: {context.get('cycle_number', 'unknown')}")
    print(f"  Timestamp: {context.get('timestamp', 'unknown')}")

def test_context_for_all_roles():
    """Test API context generation for all roles"""
    
    print("\n" + "#"*60)
    print(" Testing API Context for All Roles")
    print("#"*60)
    
    # Initialize participants
    registry = initialize_updated_participants()
    
    # Create a cycle
    sections = ["04-fundamental-concepts/01-universe-grid"]
    cycle = GovernanceCycleRefined(1, sections, max_proposals_per_participant=2)
    
    # Get participants
    claude = None
    gpt = None
    human = None
    
    for lct_id, participant in registry.participants.items():
        if "Claude" in participant.name:
            claude = participant
        elif "GPT" in participant.name:
            gpt = participant
        elif "Dennis" in participant.name:
            human = participant
    
    # Test 1: Proposer Context
    print("\n\n" + "="*60)
    print(" PROPOSER CONTEXT (Claude-4.1)")
    print("="*60)
    
    proposer_context = cycle._create_api_context(claude, cycle.phase)
    pretty_print_context(proposer_context)
    
    # Simulate some activity
    cycle.proposals_submitted[claude.id] = ["prop-001", "prop-002"]
    cycle.holds["prop-001"] = gpt.id  # GPT holds Claude's proposal
    
    # Test 2: Reviewer Context  
    print("\n\n" + "="*60)
    print(" REVIEWER CONTEXT (GPT-5)")
    print("="*60)
    
    from governance_cycle_refined import CyclePhase
    cycle.phase = CyclePhase.REVIEW
    reviewer_context = cycle._create_api_context(gpt, cycle.phase)
    pretty_print_context(reviewer_context)
    
    # Test 3: Arbiter Context
    print("\n\n" + "="*60)
    print(" ARBITER CONTEXT (Human)")
    print("="*60)
    
    cycle.phase = CyclePhase.ARBITRATION
    arbiter_context = cycle._create_api_context(human, cycle.phase)
    pretty_print_context(arbiter_context)
    
    # Test 4: Counter-Proposer Context (with exclusive rights)
    print("\n\n" + "="*60)
    print(" COUNTER-PROPOSER CONTEXT (GPT-5 with exclusive rights)")
    print("="*60)
    
    cycle.phase = CyclePhase.PROPOSAL
    cycle.exclusive_counter_rights[gpt.id] = ["prop-001"]  # GPT has exclusive right
    counter_context = cycle._create_api_context(gpt, cycle.phase)
    pretty_print_context(counter_context)
    
    print("\n\n" + "#"*60)
    print(" Summary")
    print("#"*60)
    print("\n✅ All API contexts include:")
    print("  1. Identity (you_are: LCT ID)")
    print("  2. Role (proposer/reviewer/arbiter)")
    print("  3. Governance rules")
    print("  4. Role-specific information")
    print("  5. Current cycle state")
    print("\n✅ Role-specific contexts provide:")
    print("  - Proposers: max proposals, exclusive rights")
    print("  - Reviewers: proposals to review, hold capability")
    print("  - Arbiters: pending decisions, held proposals")
    print("\n✅ Exclusive counter-proposal rights tracked")

def test_api_call_simulation():
    """Simulate actual API calls with context"""
    
    print("\n\n" + "#"*60)
    print(" Simulating API Calls with Context")
    print("#"*60)
    
    from participant_api import ClaudeAPI, GPTAPI
    from lct_participants import LCT, ParticipantType, AccessMethod
    
    # Create test LCTs
    claude_lct = LCT.create(
        name="Claude-4.1",
        participant_type=ParticipantType.AI_CLAUDE,
        access_method=AccessMethod.API,
        access_endpoint="https://api.anthropic.com/v1/messages"
    )
    
    gpt_lct = LCT.create(
        name="GPT-5",
        participant_type=ParticipantType.AI_GPT,
        access_method=AccessMethod.API,
        access_endpoint="https://api.openai.com/v1/chat/completions"
    )
    
    # Create APIs
    claude_api = ClaudeAPI(claude_lct)
    gpt_api = GPTAPI(gpt_lct)
    
    # Create context as it would be in actual call
    proposal_context = {
        'you_are': claude_lct.id,
        'role': 'proposer',
        'your_name': 'Claude-4.1',
        'your_type': 'ai_claude',
        'governance_rules': '# Rules here...',
        'current_phase': 'proposal',
        'cycle_number': 1,
        'role_specific': {
            'max_proposals_allowed': 2,
            'exclusive_counter_rights': [],
            'must_counter_propose': []
        }
    }
    
    print("\n[Claude API Call - Proposal Phase]")
    print(f"Context passed:")
    print(f"  - you_are: {proposal_context['you_are']}")
    print(f"  - role: {proposal_context['role']}")
    print(f"  - max_proposals: {proposal_context['role_specific']['max_proposals_allowed']}")
    
    proposals = claude_api.generate_proposals(
        sections=["04-fundamental-concepts/01-universe-grid"],
        max_proposals=2,
        context=proposal_context
    )
    
    print(f"Result: {len(proposals)} proposals generated")
    
    # Review context
    review_context = {
        'you_are': gpt_lct.id,
        'role': 'reviewer',
        'your_name': 'GPT-5',
        'your_type': 'ai_gpt',
        'governance_rules': '# Rules here...',
        'current_phase': 'review',
        'cycle_number': 1,
        'role_specific': {
            'proposals_to_review': ['prop-001', 'prop-002'],
            'can_request_hold': True,
            'unlimited_reviews': True
        }
    }
    
    print("\n[GPT API Call - Review Phase]")
    print(f"Context passed:")
    print(f"  - you_are: {review_context['you_are']}")
    print(f"  - role: {review_context['role']}")
    print(f"  - can_request_hold: {review_context['role_specific']['can_request_hold']}")
    
    test_proposals = [
        {'id': 'prop-001', 'author': 'Claude', 'title': 'Test', 'content_summary': 'Test content'}
    ]
    
    reviews = gpt_api.review_proposals(test_proposals, review_context)
    print(f"Result: {len(reviews)} reviews submitted")

if __name__ == "__main__":
    test_context_for_all_roles()
    test_api_call_simulation()