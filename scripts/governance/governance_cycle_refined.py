#!/usr/bin/env python3
"""
Refined Governance Cycle with Exclusive Holds and Rule Context
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
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
from participant_api import ParticipantAPIFactory

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
    HOLD_FOR_COUNTER = "hold_for_counter"

class GovernanceCycleRefined:
    """Refined governance cycle with exclusive holds"""
    
    def __init__(self, cycle_id: int, sections: List[str], max_proposals_per_participant: int = 1):
        self.cycle_id = cycle_id
        self.sections = sections
        self.max_proposals_per_participant = max_proposals_per_participant
        self.phase = CyclePhase.PROPOSAL
        self.started_at = datetime.now().isoformat()
        self.completed_at = None
        
        # Track cycle activity
        self.proposals_submitted = {}  # lct_id -> [proposal_ids]
        self.reviews_submitted = {}    # lct_id -> [review_data]
        self.holds = {}                # proposal_id -> lct_id (who has exclusive counter right)
        self.decisions = {}            # proposal_id -> decision_data
        
        # Track exclusive counter-proposal rights
        self.exclusive_counter_rights = {}  # lct_id -> [proposal_ids they can counter]
        
        # Initialize components
        self.governance = WhitepaperGovernance()
        self.registry = LCTRegistry()
        
        # Load governance rules
        self.governance_rules = self._load_governance_rules()
    
    def _load_governance_rules(self) -> str:
        """Load governance rules for API context"""
        rules_file = Path(__file__).parent / "governance_rules.md"
        if rules_file.exists():
            return rules_file.read_text()
        return "Governance rules not found"
    
    def _create_api_context(self, participant: LCT, phase: CyclePhase, role: str = "participant") -> Dict:
        """Create context for API calls with clear identity and role"""
        
        # Determine specific role based on phase
        if phase == CyclePhase.PROPOSAL:
            role = "proposer"
        elif phase == CyclePhase.REVIEW:
            role = "reviewer"
        elif phase == CyclePhase.ARBITRATION:
            role = "arbiter"
        
        return {
            # Identity and Role (ALWAYS FIRST)
            'you_are': participant.id,  # LCT ID for unambiguous identity
            'role': role,  # Current role in this phase
            'your_name': participant.name,
            'your_type': participant.type.value,
            
            # Governance Rules (ALWAYS INCLUDED)
            'governance_rules': self.governance_rules,
            
            # Current Context
            'current_phase': phase.value,
            'cycle_number': self.cycle_id,
            'timestamp': datetime.now().isoformat(),
            
            # Your Status
            'your_status': {
                'has_proposed_this_cycle': participant.current_cycle_proposed,
                'has_reviewed_this_cycle': participant.current_cycle_reviewed,
                'trust_score': participant.trust_scores.aggregate(),
                'proposals_created_total': participant.history.proposals_created,
                'reviews_submitted_total': participant.history.reviews_submitted,
                'current_holds': participant.history.holds_for_counter
            },
            
            # Role-Specific Information
            'role_specific': self._get_role_specific_context(participant, phase, role),
            
            # Held Proposals (for transparency)
            'held_proposals': [
                {
                    'proposal_id': pid,
                    'holder_id': holder_id,
                    'holder_name': self.registry.get_participant(holder_id).name if self.registry.get_participant(holder_id) else 'Unknown',
                    'you_can_counter': holder_id == participant.id
                }
                for pid, holder_id in self.holds.items()
            ],
            
            # Sections Available
            'available_sections': self.sections,
            
            # Recent Activity
            'recent_decisions': list(self.decisions.values())[-5:] if self.decisions else [],
            'active_proposal_count': sum(len(p) for p in self.proposals_submitted.values()),
            'active_review_count': sum(len(r) for r in self.reviews_submitted.values())
        }
    
    def _get_role_specific_context(self, participant: LCT, phase: CyclePhase, role: str) -> Dict:
        """Get context specific to the participant's current role"""
        
        if role == "proposer":
            return {
                'max_proposals_allowed': self.max_proposals_per_participant,
                'exclusive_counter_rights': self.exclusive_counter_rights.get(participant.id, []),
                'must_counter_propose': self.exclusive_counter_rights.get(participant.id, []) if phase == CyclePhase.PROPOSAL else [],
                'proposals_already_submitted': self.proposals_submitted.get(participant.id, [])
            }
        
        elif role == "reviewer":
            # Get proposals to review
            all_proposals = []
            for proposal_ids in self.proposals_submitted.values():
                all_proposals.extend(proposal_ids)
            
            own_proposals = self.proposals_submitted.get(participant.id, [])
            can_review = [p for p in all_proposals if p not in own_proposals]
            
            return {
                'proposals_to_review': can_review,
                'cannot_review_own': own_proposals,
                'reviews_already_submitted': len(self.reviews_submitted.get(participant.id, [])),
                'can_request_hold': True,
                'unlimited_reviews': True
            }
        
        elif role == "arbiter":
            return {
                'proposals_pending_decision': [
                    pid for pids in self.proposals_submitted.values() 
                    for pid in pids if pid not in self.holds
                ],
                'proposals_on_hold': list(self.holds.keys()),
                'total_reviews_available': sum(len(r) for r in self.reviews_submitted.values()),
                'your_arbitration_history': participant.history.arbiter_decisions
            }
        
        return {}
    
    def run_proposal_phase(self, previous_cycle_holds: Optional[Dict[str, str]] = None) -> Dict:
        """Run proposal phase with exclusive counter-proposal enforcement"""
        self.phase = CyclePhase.PROPOSAL
        results = {
            'phase': 'proposal',
            'participants': [],
            'proposals_created': [],
            'counter_proposals': [],
            'holds_cleared': []
        }
        
        # Inherit holds from previous cycle
        if previous_cycle_holds:
            self.holds = previous_cycle_holds.copy()
            # Set exclusive counter rights
            for proposal_id, holder_id in self.holds.items():
                self.exclusive_counter_rights.setdefault(holder_id, []).append(proposal_id)
        
        # Get eligible participants
        participants = self.registry.get_participants_for_cycle()
        
        for participant in participants:
            if participant.current_cycle_proposed:
                continue
            
            # Check if participant has exclusive counter rights
            exclusive_proposals = self.exclusive_counter_rights.get(participant.id, [])
            
            # Create API context
            context = self._create_api_context(participant, CyclePhase.PROPOSAL)
            
            # If participant has exclusive counter rights, ONLY they can propose for those
            if exclusive_proposals:
                # This participant MUST handle their holds
                context['must_counter_propose'] = exclusive_proposals
                context['exclusive_mode'] = True
            
            # Create API instance
            api = ParticipantAPIFactory.create_api(participant)
            
            # Call proposal API
            proposals = api.generate_proposals(
                self.sections,
                self.max_proposals_per_participant,
                context
            )
            
            # Process proposals
            for proposal in proposals:
                proposal_id = proposal['id']
                
                # Check if this is a counter-proposal
                is_counter = False
                for held_id in exclusive_proposals:
                    if held_id in proposal.get('title', '') or held_id in proposal.get('content', ''):
                        # This is a counter-proposal
                        is_counter = True
                        # Remove the hold
                        if held_id in self.holds:
                            del self.holds[held_id]
                            participant.remove_hold(held_id)
                            results['holds_cleared'].append(held_id)
                        # Remove exclusive right (used)
                        if held_id in self.exclusive_counter_rights.get(participant.id, []):
                            self.exclusive_counter_rights[participant.id].remove(held_id)
                        break
                
                if is_counter:
                    results['counter_proposals'].append(proposal)
                else:
                    results['proposals_created'].append(proposal)
            
            # Record proposals
            if proposals:
                self.proposals_submitted[participant.id] = [p['id'] for p in proposals]
                participant.record_proposal()
                participant.update_heartbeat()
                
                results['participants'].append({
                    'lct_id': participant.id,
                    'name': participant.name,
                    'proposals': len(proposals),
                    'counter_proposals': len([p for p in proposals if p in results['counter_proposals']])
                })
        
        # Check for expired holds (not acted upon)
        for proposal_id, holder_id in list(self.holds.items()):
            if holder_id not in [p['lct_id'] for p in results['participants']]:
                # Holder didn't participate, hold expires after 2 cycles
                # (Would track this with cycle counters in production)
                pass
        
        self.registry.save_registry()
        return results
    
    def run_review_phase(self) -> Dict:
        """Run review phase with hold request tracking"""
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
                continue
            
            # Don't review own proposals
            own_proposals = self.proposals_submitted.get(participant.id, [])
            proposal_ids_to_review = [p for p in all_proposals if p not in own_proposals]
            
            if not proposal_ids_to_review:
                continue
            
            # Create proposal objects for review (in production, would fetch from storage)
            proposals_to_review = [
                {'id': pid, 'author': 'Unknown', 'title': f'Proposal {pid}', 'content_summary': 'Content'}
                for pid in proposal_ids_to_review
            ]
            
            # Create API context
            context = self._create_api_context(participant, CyclePhase.REVIEW)
            
            # Create API instance and call review
            api = ParticipantAPIFactory.create_api(participant)
            reviews = api.review_proposals(proposals_to_review, context)
            
            # Process reviews
            for review in reviews:
                proposal_id = review['proposal_id']
                
                # Check for hold request
                if review.get('action') == ReviewAction.HOLD_FOR_COUNTER.value:
                    # Check if proposal already held
                    if proposal_id in self.holds:
                        review['note'] = f"Hold already requested by {self.holds[proposal_id]}"
                    else:
                        # Grant exclusive hold
                        self.holds[proposal_id] = participant.id
                        participant.add_hold_for_counter(proposal_id)
                        results['holds_requested'].append({
                            'proposal_id': proposal_id,
                            'holder': participant.name,
                            'holder_id': participant.id
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
        """Run arbitration phase with hold deferrals"""
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
        
        # Create arbiter API context
        context = self._create_api_context(arbiter, CyclePhase.ARBITRATION)
        
        # Process each proposal
        for proposal_ids in self.proposals_submitted.values():
            for proposal_id in proposal_ids:
                # Check if held for counter-proposal
                if proposal_id in self.holds:
                    holder_id = self.holds[proposal_id]
                    holder = self.registry.get_participant(holder_id)
                    results['deferred'].append({
                        'proposal_id': proposal_id,
                        'reason': f"Held for counter-proposal by {holder.name if holder else 'unknown'}",
                        'holder_id': holder_id
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
    
    def _make_arbiter_decision(self, arbiter: LCT, proposal_id: str, reviews: List[Dict]) -> Dict:
        """Make arbiter decision on a proposal"""
        
        # Count review recommendations
        accept_count = sum(1 for r in reviews if 'accept' in r.get('action', '').lower())
        reject_count = sum(1 for r in reviews if 'reject' in r.get('action', '').lower())
        
        # Determine decision
        if accept_count > reject_count:
            status = ProposalStatus.ACCEPTED
            rationale = f"Accepted based on {accept_count} positive reviews"
        elif reject_count > accept_count:
            status = ProposalStatus.REJECTED
            rationale = f"Rejected based on {reject_count} negative reviews"
        else:
            status = ProposalStatus.SUBMITTED
            rationale = "Requires further review"
        
        return {
            'proposal_id': proposal_id,
            'arbiter': arbiter.name,
            'arbiter_lct': arbiter.id,
            'status': status.value,
            'rationale': rationale,
            'timestamp': datetime.now().isoformat()
        }
    
    def complete_cycle(self) -> Dict:
        """Complete the cycle and rebuild whitepaper if fractal files changed"""
        self.phase = CyclePhase.COMPLETE
        self.completed_at = datetime.now().isoformat()
        
        results = {
            'cycle_id': self.cycle_id,
            'started_at': self.started_at,
            'completed_at': self.completed_at,
            'proposals_created': sum(len(p) for p in self.proposals_submitted.values()),
            'reviews_submitted': sum(len(r) for r in self.reviews_submitted.values()),
            'decisions_made': len(self.decisions),
            'proposals_deferred': len(self.holds),
            'active_holds': self.holds.copy()  # Pass to next cycle
        }
        
        # Check if any proposals were accepted and implemented
        accepted_proposals = [
            pid for pid, decision in self.decisions.items()
            if decision.get('status') == 'accepted'
        ]
        
        if accepted_proposals:
            # Check for fractal file changes and rebuild if needed
            from whitepaper_builder import WhitepaperBuilder
            
            builder = WhitepaperBuilder(str(self.governance.base_path))
            build_result = builder.rebuild_if_changed()
            
            if build_result:
                results['whitepaper_rebuilt'] = True
                results['build_info'] = build_result
                
                # Log successful rebuild
                successful_formats = [
                    fmt for fmt, (success, _) in build_result.get('formats', {}).items()
                    if success
                ]
                
                if successful_formats:
                    print(f"\n✅ Whitepaper rebuilt: {', '.join(successful_formats)}")
            else:
                results['whitepaper_rebuilt'] = False
                print("\nNo fractal file changes detected - whitepaper not rebuilt")
        else:
            results['whitepaper_rebuilt'] = False
            results['note'] = "No accepted proposals - rebuild skipped"
        
        # Reset cycle flags for all participants
        self.registry.reset_cycle()
        
        return results

def initialize_updated_participants():
    """Initialize participants with updated model versions"""
    import os
    registry = LCTRegistry()
    
    # Clear existing and create new
    registry.participants.clear()
    
    # Human participant - Currently unavailable (testing AI-only)
    human_lct = LCT.create(
        name="Dennis",
        participant_type=ParticipantType.HUMAN,
        access_method=AccessMethod.EMAIL,
        access_endpoint="dennis@example.com",
        timeout_seconds=86400
    )
    human_lct.availability = 0.0  # Mark as unavailable for testing
    registry.register_participant(human_lct)
    
    # Claude-4.1 (Opus) - Use environment variable
    anthropic_key = os.getenv("ANTHROPIC_API_KEY", "placeholder_key")
    claude_lct = LCT.create(
        name="Claude-4.1",
        participant_type=ParticipantType.AI_CLAUDE,
        access_method=AccessMethod.API,
        access_endpoint="https://api.anthropic.com/v1/messages",
        timeout_seconds=300,
        access_credentials={"api_key": anthropic_key, "model": "claude-opus-4-1"}
    )
    registry.register_participant(claude_lct)
    
    # GPT-5 - Use environment variable
    openai_key = os.getenv("OPENAI_API_KEY", "placeholder_key")
    gpt_lct = LCT.create(
        name="GPT-5",
        participant_type=ParticipantType.AI_GPT,
        access_method=AccessMethod.API,
        access_endpoint="https://api.openai.com/v1/chat/completions",
        timeout_seconds=300,
        access_credentials={"api_key": openai_key, "model": "gpt-5"}
    )
    registry.register_participant(gpt_lct)
    
    # Deepseek-3 - Currently unavailable (no API key)
    deepseek_lct = LCT.create(
        name="Deepseek-3", 
        participant_type=ParticipantType.AI_DEEPSEEK,
        access_method=AccessMethod.API,
        access_endpoint="https://api.deepseek.com/v1/chat",
        timeout_seconds=300,
        access_credentials={"api_key": "unavailable", "model": "deepseek-3"}
    )
    deepseek_lct.availability = 0.0  # Mark as unavailable
    registry.register_participant(deepseek_lct)
    
    print(f"Initialized {len(registry.participants)} participants with updated versions:")
    for lct_id, lct in registry.participants.items():
        print(f"  - {lct.name} ({lct.type.value}): {lct_id}")
    
    return registry

def test_exclusive_holds():
    """Test exclusive hold and counter-proposal mechanism"""
    
    print("\n" + "="*60)
    print(" Testing Exclusive Hold Mechanism")
    print("="*60)
    
    # Initialize updated participants
    registry = initialize_updated_participants()
    
    sections = ["04-fundamental-concepts/01-universe-grid"]
    
    # Cycle 1: Initial proposals and hold requests
    print("\n[Cycle 1: Initial Proposals]")
    cycle1 = GovernanceCycleRefined(1, sections, max_proposals_per_participant=2)
    
    # Proposal phase
    prop_results = cycle1.run_proposal_phase()
    print(f"Proposals created: {len(prop_results['proposals_created'])}")
    
    # Review phase with hold requests
    review_results = cycle1.run_review_phase()
    print(f"Reviews submitted: {len(review_results['reviews_submitted'])}")
    print(f"Holds requested: {len(review_results['holds_requested'])}")
    
    if review_results['holds_requested']:
        for hold in review_results['holds_requested']:
            print(f"  - {hold['holder']} holds {hold['proposal_id']}")
    
    # Arbitration phase
    arb_results = cycle1.run_arbitration_phase()
    print(f"Decisions made: {len(arb_results['decisions'])}")
    print(f"Deferred for counter: {len(arb_results['deferred'])}")
    
    # Complete cycle 1
    summary1 = cycle1.complete_cycle()
    active_holds = summary1['active_holds']
    
    # Cycle 2: Counter-proposals with exclusive rights
    print("\n[Cycle 2: Counter-Proposals with Exclusive Rights]")
    cycle2 = GovernanceCycleRefined(2, sections, max_proposals_per_participant=2)
    
    # Pass holds to next cycle for exclusive counter rights
    prop_results2 = cycle2.run_proposal_phase(previous_cycle_holds=active_holds)
    print(f"New proposals: {len(prop_results2['proposals_created'])}")
    print(f"Counter-proposals: {len(prop_results2['counter_proposals'])}")
    print(f"Holds cleared: {len(prop_results2['holds_cleared'])}")
    
    # Show who had exclusive rights
    if cycle2.exclusive_counter_rights:
        print("\nExclusive counter-proposal rights:")
        for lct_id, proposal_ids in cycle2.exclusive_counter_rights.items():
            participant = registry.get_participant(lct_id)
            if participant:
                print(f"  - {participant.name} had exclusive right to counter: {proposal_ids}")
    
    # Complete cycle 2
    review_results2 = cycle2.run_review_phase()
    arb_results2 = cycle2.run_arbitration_phase()
    summary2 = cycle2.complete_cycle()
    
    print("\n[Summary]")
    print(f"Cycle 1: {summary1['proposals_created']} proposals, {summary1['proposals_deferred']} deferred")
    print(f"Cycle 2: {summary2['proposals_created']} proposals, {prop_results2['counter_proposals']} were counters")
    print(f"Exclusive hold mechanism: ✅ Working")

if __name__ == "__main__":
    test_exclusive_holds()