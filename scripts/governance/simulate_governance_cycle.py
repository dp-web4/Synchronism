#!/usr/bin/env python3
"""
Simulate a complete governance cycle with AI participants
Creates proposals, reviews them, and processes them through the system
"""

import json
import sys
from datetime import datetime
from pathlib import Path

# Add parent directory to path
sys.path.append(str(Path(__file__).parent))

from whitepaper_governance_enhanced import (
    EnhancedWhitepaperGovernance,
    ProposalType,
    ProposalStatus,
    WithdrawalReason
)
from ai_participant_api import ClaudeParticipant, GPTParticipant, ReviewRecommendation

def simulate_cycle():
    """Simulate a complete governance cycle"""
    
    print("=" * 60)
    print("SIMULATING GOVERNANCE CYCLE")
    print("=" * 60)
    print(f"Timestamp: {datetime.now().isoformat()}")
    print()
    
    # Initialize components
    governance = EnhancedWhitepaperGovernance()
    claude = ClaudeParticipant()
    gpt = GPTParticipant()
    
    # 1. Check current state
    print("1. CURRENT STATE")
    print("-" * 40)
    active = governance.get_active_proposals()
    print(f"Active proposals: {len(active)}")
    for p in active:
        print(f"  [{p['id']}] {p['title'][:40]}... ({p['status']})")
    print()
    
    # 2. Claude creates a proposal
    print("2. CLAUDE CREATES PROPOSAL")
    print("-" * 40)
    
    # First check for duplicates
    section = "04-fundamental-concepts/07-coherence"
    title = "Enhanced Coherence Metrics"
    content = """
    Add quantitative coherence measurement framework to better evaluate
    system-wide alignment. Include cross-scale coherence tensors and
    phase correlation metrics.
    """
    
    duplicates = governance.check_duplicates(section, title, content)
    if duplicates:
        print(f"  ⚠️  Found {len(duplicates)} similar proposals")
        for dup_id, similarity in duplicates[:3]:
            print(f"     - Proposal {dup_id}: {similarity:.0%} similar")
    
    # Create proposal anyway for simulation
    proposal = claude.create_proposal(
        section_path=section,
        title=title,
        proposal_type=ProposalType.EXPANSION,
        content=content,
        rationale="Coherence is central to Synchronism but lacks quantitative framework",
        specific_changes="Add new subsection on coherence metrics with mathematical definitions"
    )
    
    print(f"  Created proposal: {proposal['id']}")
    print(f"  Title: {proposal['title']}")
    print(f"  Section: {proposal['section']}")
    if 'potential_duplicates' in proposal:
        print(f"  ⚠️  Warning: {proposal['warning']}")
    print()
    
    # 3. GPT reviews the proposal
    print("3. GPT REVIEWS PROPOSAL")
    print("-" * 40)
    
    # GPT evaluates from mathematical perspective
    recommendation, analysis = gpt.evaluate_proposal(proposal)
    
    print(f"  Recommendation: {recommendation.value}")
    print(f"  Strengths: {len(analysis['strengths'])}")
    for strength in analysis['strengths']:
        print(f"    + {strength}")
    print(f"  Concerns: {len(analysis['concerns'])}")
    for concern in analysis['concerns']:
        print(f"    - {concern}")
    
    # Submit review
    review = gpt.review_proposal(
        section_path=section,
        proposal_id=proposal['id'],
        recommendation=recommendation,
        analysis=analysis
    )
    print(f"  Review submitted for proposal {proposal['id']}")
    print()
    
    # 4. Check if proposal needs withdrawal
    print("4. LIFECYCLE MANAGEMENT")
    print("-" * 40)
    
    # Simulate author considering withdrawal
    if len(analysis['concerns']) > 2:
        print("  Claude considers withdrawing due to concerns...")
        result = governance.withdraw_proposal(
            proposal_id=proposal['id'],
            reason=WithdrawalReason.AUTHOR_REQUEST,
            withdrawer=claude.name,
            comment="Will revise based on GPT's mathematical concerns"
        )
        print(f"  Withdrawal: {result['status']}")
        print(f"  Message: {result['message']}")
    else:
        print("  Proposal remains active for further review")
    
    # 5. Check for expired proposals
    print()
    print("5. CHECKING FOR EXPIRED PROPOSALS")
    print("-" * 40)
    expired = governance.expire_old_proposals(dry_run=True)
    if expired:
        print(f"  Would expire {len(expired)} proposals")
        for p in expired:
            print(f"    - [{p['id']}] {p['title'][:40]}...")
    else:
        print("  No proposals to expire")
    
    # 6. Final state
    print()
    print("6. FINAL STATE")
    print("-" * 40)
    active = governance.get_active_proposals()
    print(f"Active proposals: {len(active)}")
    for p in active:
        status = p['status']
        if 'withdrawal' in p:
            status = "withdrawn"
        print(f"  [{p['id']}] {p['title'][:40]}... ({status})")
    
    # Load and show archive
    archive_file = Path(__file__).parent / "config" / "whitepaper_proposals_archive.json"
    if archive_file.exists():
        with open(archive_file, 'r') as f:
            archive = json.load(f)
        print(f"Archived proposals: {len(archive['archived_proposals'])}")
    
    print()
    print("=" * 60)
    print("SIMULATION COMPLETE")
    print("=" * 60)

if __name__ == "__main__":
    simulate_cycle()