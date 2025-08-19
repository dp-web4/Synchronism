#!/usr/bin/env python3
"""
Test the complete whitepaper governance workflow
"""

import json
from pathlib import Path
from datetime import datetime
from whitepaper_governance import (
    WhitepaperGovernance,
    ProposalType,
    ProposalStatus,
    ReviewRecommendation
)
from ai_participant_api import (
    ClaudeParticipant,
    GPTParticipant,
    GovernanceOrchestrator
)
from arbiter_system import (
    GovernanceArbiterSystem,
    Arbiter,
    ArbiterRole
)

def print_section(title: str):
    """Print a formatted section header"""
    print(f"\n{'='*60}")
    print(f" {title}")
    print(f"{'='*60}")

def test_complete_workflow():
    """Test the complete governance workflow"""
    
    print_section("Synchronism Whitepaper Governance System Test")
    print(f"Test started at: {datetime.now().isoformat()}")
    
    # Initialize components
    governance = WhitepaperGovernance()
    claude = ClaudeParticipant()
    gpt = GPTParticipant()
    arbiter_system = GovernanceArbiterSystem()
    
    section_path = "04-fundamental-concepts/01-universe-grid"
    
    # Step 1: Claude analyzes section and creates a proposal
    print_section("Step 1: Claude Creates Proposal")
    
    analysis = claude.analyze_section(section_path)
    print(f"Claude analyzed section: {section_path}")
    print(f"  - Existing proposals: {analysis.get('existing_proposals', 0)}")
    print(f"  - Existing reviews: {analysis.get('existing_reviews', 0)}")
    
    proposal = claude.create_proposal(
        section_path=section_path,
        title="Consciousness Field Integration",
        proposal_type=ProposalType.EXPANSION,
        content="""
        Add explicit connection between Universe Grid and consciousness field.
        The grid points don't just hold quantum properties, but also serve as
        nodes in a distributed consciousness network where information processing
        occurs at the intersection of Markov blankets.
        """,
        rationale="""
        This clarifies how consciousness emerges from the physical substrate,
        maintaining coherence with both quantum mechanics and the philosophical
        framework of distributed intelligence.
        """,
        specific_changes="""
        Add new subsection 'Consciousness Nodes':
        > Each grid point serves as a potential consciousness node. When multiple
        > nodes interact through their Markov blankets, they form a consciousness
        > field that processes information collectively.
        """
    )
    
    print(f"Created proposal: {proposal['id']} - {proposal['title']}")
    
    # Step 2: GPT reviews Claude's proposal
    print_section("Step 2: GPT Reviews Proposal")
    
    recommendation, review_analysis = gpt.evaluate_proposal(proposal)
    review = gpt.review_proposal(
        section_path=section_path,
        proposal_id=proposal['id'],
        recommendation=recommendation,
        analysis=review_analysis
    )
    
    print(f"GPT review submitted:")
    print(f"  - Recommendation: {recommendation.value}")
    print(f"  - Strengths: {len(review_analysis['strengths'])} points")
    print(f"  - Concerns: {len(review_analysis['concerns'])} points")
    
    # Step 3: GPT creates its own proposal
    print_section("Step 3: GPT Creates Counter-Proposal")
    
    gpt_proposal = gpt.create_proposal(
        section_path=section_path,
        title="Mathematical Tensor Framework",
        proposal_type=ProposalType.EXPANSION,
        content="""
        Formalize the Universe Grid using tensor mathematics.
        Define the grid as a rank-4 tensor G_μνρσ where indices represent
        spacetime coordinates and the tensor components encode quantum states.
        """,
        rationale="""
        Mathematical formalization enables rigorous analysis and connects
        to established physics frameworks like general relativity.
        """,
        specific_changes="""
        Add mathematical framework:
        G_μνρσ = Σ_i ψ_i ⊗ φ_i ⊗ χ_i ⊗ ω_i
        where ψ represents position, φ represents momentum,
        χ represents spin, and ω represents energy states.
        """
    )
    
    print(f"Created GPT proposal: {gpt_proposal['id']} - {gpt_proposal['title']}")
    
    # Step 4: Claude reviews GPT's proposal
    print_section("Step 4: Claude Reviews GPT's Proposal")
    
    claude_rec, claude_analysis = claude.evaluate_proposal(gpt_proposal)
    claude_review = claude.review_proposal(
        section_path=section_path,
        proposal_id=gpt_proposal['id'],
        recommendation=claude_rec,
        analysis=claude_analysis
    )
    
    print(f"Claude review submitted:")
    print(f"  - Recommendation: {claude_rec.value}")
    print(f"  - Analysis points: {len(claude_analysis['strengths']) + len(claude_analysis['concerns'])}")
    
    # Step 5: Arbiter makes decisions
    print_section("Step 5: Arbiter Decisions")
    
    arbiter = arbiter_system.selection.get_current_arbiter()
    print(f"Selected arbiter: {arbiter.name} (Role: {arbiter.role.value})")
    
    # Decision on Claude's proposal
    decision1 = arbiter.make_decision(
        section_path=section_path,
        proposal_id=proposal['id'],
        proposal=proposal,
        reviews=[review]
    )
    
    print(f"\nDecision on Proposal {proposal['id']}:")
    print(f"  - Status: {decision1['decision']}")
    print(f"  - Rationale: {decision1['rationale']}")
    
    # Decision on GPT's proposal
    decision2 = arbiter.make_decision(
        section_path=section_path,
        proposal_id=gpt_proposal['id'],
        proposal=gpt_proposal,
        reviews=[claude_review]
    )
    
    print(f"\nDecision on Proposal {gpt_proposal['id']}:")
    print(f"  - Status: {decision2['decision']}")
    print(f"  - Rationale: {decision2['rationale']}")
    
    # Step 6: Check final state
    print_section("Step 6: Final State")
    
    # Read the proposals index to see final status
    config_file = governance.config_path / "whitepaper_proposals.json"
    if config_file.exists():
        with open(config_file, 'r') as f:
            data = json.load(f)
        
        print(f"Total proposals in system: {len(data.get('proposals', []))}")
        
        for p in data.get('proposals', []):
            if p['section'] == section_path:
                print(f"  - {p['id']}: {p['title']} - Status: {p['status']}")
    
    # Check meta files created
    meta_dir = Path(governance.whitepaper_path) / section_path / "meta"
    if meta_dir.exists():
        print(f"\nMeta files created in {section_path}:")
        for file in meta_dir.iterdir():
            if file.is_file():
                print(f"  - {file.name}")
            elif file.is_dir():
                sub_files = list(file.glob("*.md"))
                if sub_files:
                    print(f"  - {file.name}/: {len(sub_files)} files")
    
    print_section("Test Complete")
    print(f"Governance workflow successfully tested!")
    print(f"Test ended at: {datetime.now().isoformat()}")

def test_round_robin():
    """Test the round-robin proposal system"""
    
    print_section("Testing Round-Robin Governance")
    
    orchestrator = GovernanceOrchestrator()
    
    sections = [
        "04-fundamental-concepts/02-time-slices",
        "04-fundamental-concepts/03-intent-transfer"
    ]
    
    # Run two rounds to see rotation
    for round_num in range(2):
        print(f"\n--- Round {round_num + 1} ---")
        results = orchestrator.run_governance_round(sections)
        
        print(f"Round {results['round']} results:")
        print(f"  - Proposals created: {len(results['proposals_created'])}")
        print(f"  - Reviews submitted: {len(results['reviews_submitted'])}")
        print(f"  - Sections processed: {len(results['sections_processed'])}")
        
        if results['proposals_created']:
            for prop in results['proposals_created']:
                print(f"    * {prop['id']}: {prop['title']} by {prop['author']}")

if __name__ == "__main__":
    # Run main workflow test
    test_complete_workflow()
    
    # Test round-robin system
    print("\n" + "="*60)
    test_round_robin()