#!/usr/bin/env python3
"""
Governance Cycle Management with LCT-based participants
Implements controlled cycles with proposal limits and hold mechanisms
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from enum import Enum

from lct_participants import (
    LCT, LCTRegistry, ParticipantType, AccessMethod
)
from whitepaper_governance import (
    WhitepaperGovernance,
    ProposalType,
    ProposalStatus,
    ReviewRecommendation
)

class CyclePhase(Enum):
    PROPOSAL = "proposal"
    REVIEW = "review"
    ARBITRATION = "arbitration"
    IMPLEMENTATION = "implementation"
    COMPLETE = "complete"

class ReviewAction(Enum):
    ACCEPT = "accept"
    ACCEPT_WITH_REVISIONS = "accept_with_revisions"
    REVISE_AND_RESUBMIT = "revise_and_resubmit"
    REJECT = "reject"
    HOLD_FOR_COUNTER = "hold_for_counter"  # New: defer decision for counter-proposal

class GovernanceCycle:
    """Manages a single governance cycle"""
    
    def __init__(self, cycle_id: int, sections: List[str], max_proposals_per_participant: int = 1,
                 held_proposals: Dict = None):
        self.cycle_id = cycle_id
        self.sections = sections
        self.max_proposals_per_participant = max_proposals_per_participant
        self.phase = CyclePhase.PROPOSAL
        self.started_at = datetime.now().isoformat()
        self.completed_at = None
        
        # Track cycle activity
        self.proposals_submitted = {}  # lct_id -> [proposal_ids]
        self.reviews_submitted = {}    # lct_id -> [review_data]
        self.holds = {}                # proposal_id -> lct_id (who wants to counter)
        self.decisions = {}            # proposal_id -> decision_data
        
        # Track held proposals for counter-proposal mechanism
        self.held_proposals = held_proposals or {}  # proposal_id -> original_data
        self.counter_proposals = {}    # original_id -> enhanced_proposal
        
        # Initialize components
        self.governance = WhitepaperGovernance()
        self.registry = LCTRegistry()
    
    def run_proposal_phase(self) -> Dict:
        """Run the proposal phase of the cycle"""
        self.phase = CyclePhase.PROPOSAL
        results = {
            'phase': 'proposal',
            'participants': [],
            'proposals_created': [],
            'counter_proposals': []
        }
        
        # First, handle counter-proposals for held items
        for proposal_id, holder_id in self.held_proposals.items():
            holder = self.registry.get_participant(holder_id)
            if holder and holder.check_availability():
                # Create counter-proposal that enhances the original
                counter = self._create_counter_proposal(holder, proposal_id)
                if counter:
                    self.counter_proposals[proposal_id] = counter
                    results['counter_proposals'].append(counter['id'])
        
        # Get eligible participants for new proposals
        participants = self.registry.get_participants_for_cycle()
        
        for participant in participants:
            if participant.current_cycle_proposed:
                continue  # Already proposed this cycle
            
            # Call participant's proposal API
            proposals = self._call_proposal_api(
                participant,
                self.sections,
                self.max_proposals_per_participant
            )
            
            # Record proposals
            if proposals:
                self.proposals_submitted[participant.id] = [p['id'] for p in proposals]
                participant.record_proposal()
                participant.update_heartbeat()
                
                results['proposals_created'].extend(proposals)
                results['participants'].append({
                    'lct_id': participant.id,
                    'name': participant.name,
                    'proposals': len(proposals)
                })
        
        self.registry.save_registry()
        return results
    
    def run_review_phase(self) -> Dict:
        """Run the review phase of the cycle"""
        self.phase = CyclePhase.REVIEW
        results = {
            'phase': 'review',
            'participants': [],
            'reviews_submitted': [],
            'holds_requested': []
        }
        
        # Get all proposals from this cycle
        all_proposals = []
        for proposal_ids in self.proposals_submitted.values():
            all_proposals.extend(proposal_ids)
        
        # Get eligible reviewers
        participants = self.registry.get_participants_for_cycle()
        
        for participant in participants:
            if participant.current_cycle_reviewed:
                continue  # Already reviewed this cycle
            
            # Don't review own proposals
            own_proposals = self.proposals_submitted.get(participant.id, [])
            proposals_to_review = [p for p in all_proposals if p not in own_proposals]
            
            if not proposals_to_review:
                continue
            
            # Call review API (can review as many as they want)
            reviews = self._call_review_api(
                participant,
                proposals_to_review
            )
            
            # Process reviews
            for review in reviews:
                proposal_id = review['proposal_id']
                
                # Check for hold request
                if review.get('action') == ReviewAction.HOLD_FOR_COUNTER.value:
                    self.holds[proposal_id] = participant.id
                    participant.add_hold_for_counter(proposal_id)
                    results['holds_requested'].append({
                        'proposal_id': proposal_id,
                        'holder': participant.name
                    })
                
                # Record review
                self.reviews_submitted.setdefault(participant.id, []).append(review)
                results['reviews_submitted'].append(review)
            
            if reviews:
                participant.record_review()
                participant.update_heartbeat()
        
        results['participants'] = [
            {
                'lct_id': p.id,
                'name': p.name,
                'reviews': len(self.reviews_submitted.get(p.id, []))
            }
            for p in participants if p.id in self.reviews_submitted
        ]
        
        self.registry.save_registry()
        return results
    
    def run_arbitration_phase(self) -> Dict:
        """Run the arbitration phase"""
        self.phase = CyclePhase.ARBITRATION
        results = {
            'phase': 'arbitration',
            'arbiter': None,
            'decisions': [],
            'deferred': []
        }
        
        # Select arbiter with fallback
        arbiter = self.registry.select_arbiter_with_fallback()
        if not arbiter:
            results['error'] = "No available arbiter"
            return results
        
        results['arbiter'] = {
            'lct_id': arbiter.id,
            'name': arbiter.name,
            'trust_score': arbiter.trust_scores.aggregate()
        }
        
        # Process each proposal
        for proposal_ids in self.proposals_submitted.values():
            for proposal_id in proposal_ids:
                # Check if held for counter-proposal
                if proposal_id in self.holds:
                    holder_id = self.holds[proposal_id]
                    holder = self.registry.get_participant(holder_id)
                    results['deferred'].append({
                        'proposal_id': proposal_id,
                        'reason': f"Held for counter-proposal by {holder.name if holder else 'unknown'}"
                    })
                    continue
                
                # Get reviews for this proposal
                proposal_reviews = []
                for reviews in self.reviews_submitted.values():
                    proposal_reviews.extend([r for r in reviews if r.get('proposal_id') == proposal_id])
                
                # Make arbiter decision
                decision = self._make_arbiter_decision(
                    arbiter,
                    proposal_id,
                    proposal_reviews
                )
                
                self.decisions[proposal_id] = decision
                results['decisions'].append(decision)
        
        arbiter.record_arbiter_decision()
        self.registry.save_registry()
        return results
    
    def complete_cycle(self) -> Dict:
        """Complete the cycle and prepare for next"""
        self.phase = CyclePhase.COMPLETE
        self.completed_at = datetime.now().isoformat()
        
        results = {
            'cycle_id': self.cycle_id,
            'started_at': self.started_at,
            'completed_at': self.completed_at,
            'proposals_created': sum(len(p) for p in self.proposals_submitted.values()),
            'reviews_submitted': sum(len(r) for r in self.reviews_submitted.values()),
            'decisions_made': len(self.decisions),
            'proposals_deferred': len(self.holds)
        }
        
        # Reset cycle flags for all participants
        self.registry.reset_cycle()
        
        return results
    
    def _call_proposal_api(self, participant: LCT, sections: List[str], max_proposals: int) -> List[Dict]:
        """Call participant's API to get proposals"""
        
        # Simulate API call based on participant type
        proposals = []
        
        if participant.type == ParticipantType.AI_CLAUDE:
            # Claude focuses on conceptual enhancements
            if len(proposals) < max_proposals:
                proposal = self.governance.create_proposal(
                    section_path=sections[0],
                    author=f"{participant.name} (LCT: {participant.id[:8]})",
                    title=f"Conceptual Enhancement - Cycle {self.cycle_id}",
                    proposal_type=ProposalType.CLARIFICATION,
                    content="Enhance philosophical coherence",
                    rationale="Maintains Synchronism framework consistency",
                    specific_changes="Details from Claude API"
                )
                proposals.append(proposal)
        
        elif participant.type == ParticipantType.AI_GPT:
            # GPT focuses on mathematical formalization
            if len(proposals) < max_proposals:
                proposal = self.governance.create_proposal(
                    section_path=sections[0],
                    author=f"{participant.name} (LCT: {participant.id[:8]})",
                    title=f"Mathematical Framework - Cycle {self.cycle_id}",
                    proposal_type=ProposalType.EXPANSION,
                    content="Add mathematical rigor",
                    rationale="Provides formal foundation",
                    specific_changes="Equations from GPT API"
                )
                proposals.append(proposal)
        
        # Participant can return 0 to max_proposals
        return proposals[:max_proposals]
    
    def _call_review_api(self, participant: LCT, proposal_ids: List[str]) -> List[Dict]:
        """Call participant's API to get reviews (unlimited)"""
        
        reviews = []
        
        for proposal_id in proposal_ids:
            # Simulate review decision based on participant preferences
            # In production, this would call actual AI API
            
            # Check if participant wants to counter-propose (simulated logic)
            wants_counter = (participant.type == ParticipantType.AI_CLAUDE and "Mathematical" in proposal_id) or \
                          (participant.type == ParticipantType.AI_GPT and "Conceptual" in proposal_id)
            
            if wants_counter:
                # Request hold for counter-proposal
                review = {
                    'proposal_id': proposal_id,
                    'reviewer': participant.name,
                    'reviewer_id': participant.id,
                    'action': ReviewAction.HOLD_FOR_COUNTER.value,
                    'comment': "Requesting hold to prepare counter-proposal in next cycle",
                    'strengths': ["Addresses important topic"],
                    'concerns': ["Alternative approach may be better"],
                    'signed': participant.id,  # Digital signature
                    'timestamp': datetime.now().isoformat()
                }
            else:
                # Normal review - simplified to accept or hold
                review = {
                    'proposal_id': proposal_id,
                    'reviewer': participant.name,
                    'reviewer_id': participant.id,
                    'action': ReviewAction.ACCEPT.value,
                    'comment': "Well-structured proposal with clear impact",
                    'strengths': ["Good approach", "Well reasoned"],
                    'concerns': [],
                    'signed': participant.id,  # Digital signature
                    'timestamp': datetime.now().isoformat()
                }
            
            reviews.append(review)
        
        return reviews
    
    def _create_counter_proposal(self, participant: LCT, original_proposal_id: str) -> Dict:
        """Create a counter-proposal that enhances/modifies the original
        
        This creates a collaborative enhancement where the counter-proposer
        adds to the original without removing review tags, enabling
        multi-participant conversation.
        """
        
        # In production, would fetch original proposal from storage
        original = self.held_proposals.get(original_proposal_id, {})
        
        # Create enhanced version with participant's modifications
        counter_proposal = {
            'id': f"{original_proposal_id}_counter_{participant.id[:8]}",
            'original_id': original_proposal_id,
            'type': 'counter_proposal',
            'author': participant.name,
            'author_id': participant.id,
            'conversation': [
                {
                    'participant': original.get('author', 'Unknown'),
                    'content': original.get('content', 'Original proposal content')
                },
                {
                    'participant': participant.name,
                    'content': f"Enhancement: Building on the original proposal with additional considerations",
                    'modifications': [
                        "Added mathematical framework",
                        "Clarified implementation details",
                        "Addressed edge cases"
                    ]
                }
            ],
            'metadata': {
                'title': f"Enhanced: {original.get('title', 'Untitled')}",
                'section': original.get('section', 'Unknown'),
                'iteration': len(original.get('conversation', [])) + 1
            },
            'timestamp': datetime.now().isoformat(),
            'signed': participant.id
        }
        
        # Remove review tags so it appears fresh to new reviewers
        counter_proposal['reviews_cleared'] = True
        counter_proposal['previous_reviews'] = original.get('reviews', [])
        
        return counter_proposal
    
    def _make_arbiter_decision(self, arbiter: LCT, proposal_id: str, reviews: List[Dict]) -> Dict:
        """Make arbiter decision on a proposal
        
        Rules:
        - Accept only if there's at least one review AND all reviews are 'accept'
        - Otherwise defer for next cycle
        """
        
        if not reviews:
            # No reviews = no decision
            status = ProposalStatus.SUBMITTED
            rationale = "No reviews received - deferred to next cycle"
        else:
            # Check if all reviews are accept
            all_accept = all(
                r.get('action') == ReviewAction.ACCEPT.value 
                for r in reviews
            )
            
            if all_accept:
                status = ProposalStatus.ACCEPTED
                rationale = f"Unanimously accepted by {len(reviews)} reviewer(s)"
                # Record signatures of accepting reviewers
                signatures = [r.get('signed', 'unsigned') for r in reviews]
                rationale += f" [Signed: {', '.join(signatures[:3])}...]"
            else:
                # At least one hold or other action
                hold_count = sum(1 for r in reviews if r.get('action') == ReviewAction.HOLD_FOR_COUNTER.value)
                if hold_count > 0:
                    status = ProposalStatus.SUBMITTED
                    rationale = f"Held for counter-proposal by {hold_count} reviewer(s)"
                else:
                    status = ProposalStatus.SUBMITTED
                    rationale = "Mixed reviews - requires further discussion"
        
        return {
            'proposal_id': proposal_id,
            'arbiter': arbiter.name,
            'arbiter_lct': arbiter.id,
            'status': status.value,
            'rationale': rationale,
            'reviews_count': len(reviews),
            'timestamp': datetime.now().isoformat()
        }

class GovernanceOrchestrator:
    """Orchestrates multiple governance cycles"""
    
    def __init__(self, base_path: str = "/mnt/c/exe/projects/ai-agents/Synchronism"):
        self.base_path = Path(base_path)
        self.cycles_file = self.base_path / "scripts" / "governance" / "config" / "cycles.json"
        self.current_cycle_id = 0
        self.cycles_history = []
        self.held_proposals = {}  # Track proposals held for counter across cycles
        self.load_history()
    
    def load_history(self):
        """Load cycle history"""
        if self.cycles_file.exists():
            with open(self.cycles_file, 'r') as f:
                data = json.load(f)
                self.current_cycle_id = data.get('current_cycle_id', 0)
                self.cycles_history = data.get('history', [])
    
    def save_history(self):
        """Save cycle history"""
        self.cycles_file.parent.mkdir(parents=True, exist_ok=True)
        
        data = {
            'current_cycle_id': self.current_cycle_id,
            'history': self.cycles_history,
            'last_updated': datetime.now().isoformat()
        }
        
        with open(self.cycles_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    def run_cycle(self, sections: List[str], max_proposals: int = 1) -> Dict:
        """Run a complete governance cycle"""
        
        self.current_cycle_id += 1
        cycle = GovernanceCycle(
            cycle_id=self.current_cycle_id,
            sections=sections,
            max_proposals_per_participant=max_proposals,
            held_proposals=self.held_proposals  # Pass held proposals from previous cycle
        )
        
        print(f"\n{'='*60}")
        print(f" Governance Cycle {self.current_cycle_id}")
        print(f"{'='*60}")
        
        # Run phases
        results = {'cycle_id': self.current_cycle_id, 'phases': {}}
        
        # Proposal phase
        print("\n[Phase 1: Proposals]")
        proposal_results = cycle.run_proposal_phase()
        results['phases']['proposal'] = proposal_results
        print(f"  - {len(proposal_results['proposals_created'])} proposals created")
        
        # Review phase
        print("\n[Phase 2: Reviews]")
        review_results = cycle.run_review_phase()
        results['phases']['review'] = review_results
        print(f"  - {len(review_results['reviews_submitted'])} reviews submitted")
        print(f"  - {len(review_results['holds_requested'])} holds requested")
        
        # Arbitration phase
        print("\n[Phase 3: Arbitration]")
        arbitration_results = cycle.run_arbitration_phase()
        results['phases']['arbitration'] = arbitration_results
        if arbitration_results.get('arbiter'):
            print(f"  - Arbiter: {arbitration_results['arbiter']['name']}")
            print(f"  - {len(arbitration_results['decisions'])} decisions made")
            print(f"  - {len(arbitration_results['deferred'])} proposals deferred")
        
        # Complete cycle
        completion = cycle.complete_cycle()
        results['summary'] = completion
        
        # Update held proposals for next cycle
        # Collect proposals that were held for counter
        self.held_proposals.clear()
        for proposal_id, holder_id in cycle.holds.items():
            # In production, would fetch full proposal data from storage
            self.held_proposals[proposal_id] = holder_id
        
        # Add counter-proposals info to results
        if cycle.counter_proposals:
            results['counter_proposals'] = cycle.counter_proposals
        
        # Save to history
        self.cycles_history.append(results)
        self.save_history()
        
        return results
    
    def run_multiple_cycles(self, sections: List[str], num_cycles: int = 3, max_proposals: int = 2):
        """Run multiple governance cycles to show back-and-forth"""
        
        for i in range(num_cycles):
            print(f"\n{'#'*60}")
            print(f" Running Cycle {i+1} of {num_cycles}")
            print(f"{'#'*60}")
            
            results = self.run_cycle(sections, max_proposals)
            
            # Show holds carrying forward
            if results['phases']['review'].get('holds_requested'):
                print("\nProposals held for counter in next cycle:")
                for hold in results['phases']['review']['holds_requested']:
                    print(f"  - {hold['proposal_id']} by {hold['holder']}")

def main():
    """Test the governance cycle system"""
    
    # Initialize participants
    from lct_participants import initialize_default_participants
    registry = initialize_default_participants()
    
    # Create orchestrator
    orchestrator = GovernanceOrchestrator()
    
    # Define sections to govern
    sections = [
        "04-fundamental-concepts/01-universe-grid",
        "04-fundamental-concepts/02-time-slices"
    ]
    
    # Run multiple cycles showing back-and-forth
    orchestrator.run_multiple_cycles(
        sections=sections,
        num_cycles=3,
        max_proposals=2  # Limit proposals per participant per cycle
    )
    
    print("\n" + "="*60)
    print(" Governance Cycles Complete")
    print("="*60)
    print(f"Total cycles run: {orchestrator.current_cycle_id}")
    print(f"History saved to: {orchestrator.cycles_file}")

if __name__ == "__main__":
    main()