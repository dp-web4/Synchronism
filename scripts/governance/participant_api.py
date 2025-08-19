#!/usr/bin/env python3
"""
Participant API for LCT-based Governance
Provides clean API interface for AI and human participants
"""

import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from abc import ABC, abstractmethod

from lct_participants import LCT, ParticipantType, AccessMethod
from whitepaper_governance import (
    WhitepaperGovernance,
    ProposalType,
    ReviewRecommendation
)
from governance_cycle import ReviewAction

class ParticipantAPI(ABC):
    """Abstract base class for participant APIs"""
    
    def __init__(self, lct: LCT):
        self.lct = lct
        self.governance = WhitepaperGovernance()
    
    @abstractmethod
    def generate_proposals(self, 
                          sections: List[str], 
                          max_proposals: int,
                          context: Optional[Dict] = None) -> List[Dict]:
        """Generate 0 to max_proposals proposals"""
        pass
    
    @abstractmethod
    def review_proposals(self,
                        proposals: List[Dict],
                        context: Optional[Dict] = None) -> List[Dict]:
        """Review proposals (unlimited reviews allowed)"""
        pass
    
    def check_timeout(self) -> bool:
        """Check if participant has timed out"""
        return self.lct.check_availability()
    
    def update_activity(self):
        """Update last activity timestamp"""
        self.lct.update_heartbeat()

class ClaudeAPI(ParticipantAPI):
    """Claude-specific API implementation"""
    
    def generate_proposals(self, 
                          sections: List[str], 
                          max_proposals: int,
                          context: Optional[Dict] = None) -> List[Dict]:
        """Claude generates conceptual/philosophical proposals"""
        
        proposals = []
        
        # Check for counter-proposal opportunities
        if context and 'held_proposals' in context:
            for held in context['held_proposals']:
                if self.lct.id in held.get('holder_id', ''):
                    # Generate counter-proposal
                    counter = self._generate_counter_proposal(held, sections[0])
                    if counter:
                        proposals.append(counter)
                        self.lct.remove_hold(held['proposal_id'])
                        if len(proposals) >= max_proposals:
                            break
        
        # Generate new proposals if under limit
        remaining = max_proposals - len(proposals)
        if remaining > 0:
            # Analyze sections for conceptual improvements
            for section in sections[:remaining]:
                analysis = self._analyze_section_conceptual(section)
                if analysis.get('needs_improvement'):
                    proposal = self.governance.create_proposal(
                        section_path=section,
                        author=f"{self.lct.name} (LCT: {self.lct.id[:8]})",
                        title=analysis['improvement_title'],
                        proposal_type=ProposalType.CLARIFICATION,
                        content=analysis['improvement_content'],
                        rationale=analysis['rationale'],
                        specific_changes=analysis['changes']
                    )
                    proposals.append(proposal)
                    
                    if len(proposals) >= max_proposals:
                        break
        
        # Update activity
        self.update_activity()
        
        # Can return 0 to max_proposals
        return proposals[:max_proposals]
    
    def review_proposals(self,
                        proposals: List[Dict],
                        context: Optional[Dict] = None) -> List[Dict]:
        """Claude reviews with focus on philosophical coherence"""
        
        reviews = []
        
        for proposal in proposals:
            # Don't review own proposals
            if self.lct.id[:8] in proposal.get('author', ''):
                continue
            
            review = self._evaluate_philosophical_coherence(proposal)
            
            # Check if counter-proposal would be better
            if self._should_counter_propose(proposal):
                review['action'] = ReviewAction.HOLD_FOR_COUNTER.value
                review['comment'] = "Philosophical framework could be better aligned. Requesting hold for counter-proposal."
                self.lct.add_hold_for_counter(proposal['id'])
            
            reviews.append(review)
        
        self.update_activity()
        return reviews
    
    def _analyze_section_conceptual(self, section: str) -> Dict:
        """Analyze section for conceptual improvements"""
        # Simulate analysis
        return {
            'needs_improvement': True,  # For demo, always find something
            'improvement_title': f"Consciousness Integration for {section.split('/')[-1]}",
            'improvement_content': "Enhance connection to consciousness field theory",
            'rationale': "Strengthens philosophical coherence with Synchronism framework",
            'changes': "Add subsection on consciousness emergence patterns"
        }
    
    def _evaluate_philosophical_coherence(self, proposal: Dict) -> Dict:
        """Evaluate proposal's philosophical alignment"""
        
        # Check for key philosophical concepts
        content = proposal.get('content_summary', '')
        has_consciousness = 'consciousness' in content.lower()
        has_markov = 'markov' in content.lower()
        has_fractal = 'fractal' in content.lower()
        
        strengths = []
        concerns = []
        
        if has_consciousness:
            strengths.append("Addresses consciousness emergence")
        if has_markov:
            strengths.append("Maintains Markov blanket boundaries")
        if has_fractal:
            strengths.append("Preserves fractal self-similarity")
        
        if not (has_consciousness or has_markov or has_fractal):
            concerns.append("Lacks connection to core philosophical concepts")
        
        # Determine recommendation
        if len(strengths) > len(concerns):
            recommendation = ReviewRecommendation.ACCEPT_WITH_REVISIONS
        else:
            recommendation = ReviewRecommendation.REVISE_AND_RESUBMIT
        
        return {
            'proposal_id': proposal['id'],
            'reviewer': f"{self.lct.name} (LCT: {self.lct.id[:8]})",
            'action': ReviewAction.ACCEPT_WITH_REVISIONS.value,
            'strengths': strengths,
            'concerns': concerns,
            'recommendation': recommendation.value,
            'timestamp': datetime.now().isoformat()
        }
    
    def _should_counter_propose(self, proposal: Dict) -> bool:
        """Determine if a counter-proposal would be beneficial"""
        # Counter if it's heavily mathematical without philosophy
        return 'mathematical' in proposal.get('title', '').lower() and \
               'consciousness' not in proposal.get('content_summary', '').lower()
    
    def _generate_counter_proposal(self, held_proposal: Dict, section: str) -> Optional[Dict]:
        """Generate a counter-proposal to a held proposal"""
        return self.governance.create_proposal(
            section_path=section,
            author=f"{self.lct.name} (LCT: {self.lct.id[:8]})",
            title=f"Counter: Philosophical Framework for {held_proposal.get('title', 'Unknown')}",
            proposal_type=ProposalType.EXPANSION,
            content="Alternative approach emphasizing consciousness and emergence",
            rationale="Provides philosophical grounding for technical concepts",
            specific_changes="Add consciousness-based interpretation"
        )

class GPTAPI(ParticipantAPI):
    """GPT-specific API implementation"""
    
    def generate_proposals(self, 
                          sections: List[str], 
                          max_proposals: int,
                          context: Optional[Dict] = None) -> List[Dict]:
        """GPT generates mathematical/technical proposals"""
        
        proposals = []
        
        # Check for counter-proposal opportunities
        if context and 'held_proposals' in context:
            for held in context['held_proposals']:
                if self.lct.id in held.get('holder_id', ''):
                    # Generate counter-proposal
                    counter = self._generate_counter_proposal(held, sections[0])
                    if counter:
                        proposals.append(counter)
                        self.lct.remove_hold(held['proposal_id'])
                        if len(proposals) >= max_proposals:
                            break
        
        # Generate new proposals if under limit
        remaining = max_proposals - len(proposals)
        if remaining > 0:
            for section in sections[:remaining]:
                analysis = self._analyze_section_mathematical(section)
                if analysis.get('needs_formalization'):
                    proposal = self.governance.create_proposal(
                        section_path=section,
                        author=f"{self.lct.name} (LCT: {self.lct.id[:8]})",
                        title=analysis['formalization_title'],
                        proposal_type=ProposalType.EXPANSION,
                        content=analysis['formalization_content'],
                        rationale=analysis['rationale'],
                        specific_changes=analysis['equations']
                    )
                    proposals.append(proposal)
                    
                    if len(proposals) >= max_proposals:
                        break
        
        self.update_activity()
        return proposals[:max_proposals]
    
    def review_proposals(self,
                        proposals: List[Dict],
                        context: Optional[Dict] = None) -> List[Dict]:
        """GPT reviews with focus on mathematical rigor"""
        
        reviews = []
        
        for proposal in proposals:
            # Don't review own proposals
            if self.lct.id[:8] in proposal.get('author', ''):
                continue
            
            review = self._evaluate_mathematical_rigor(proposal)
            
            # Check if mathematical formalization needed
            if self._needs_mathematical_counter(proposal):
                review['action'] = ReviewAction.HOLD_FOR_COUNTER.value
                review['comment'] = "Requires mathematical formalization. Requesting hold for counter-proposal."
                self.lct.add_hold_for_counter(proposal['id'])
            
            reviews.append(review)
        
        self.update_activity()
        return reviews
    
    def _analyze_section_mathematical(self, section: str) -> Dict:
        """Analyze section for mathematical improvements"""
        return {
            'needs_formalization': True,
            'formalization_title': f"Tensor Framework for {section.split('/')[-1]}",
            'formalization_content': "Mathematical formalization using tensor algebra",
            'rationale': "Provides rigorous foundation for analysis",
            'equations': "T_μνρσ = Σ ψ_i ⊗ φ_i with constraints..."
        }
    
    def _evaluate_mathematical_rigor(self, proposal: Dict) -> Dict:
        """Evaluate proposal's mathematical accuracy"""
        
        content = proposal.get('content_summary', '')
        has_equations = 'equation' in content.lower() or 'tensor' in content.lower()
        has_proofs = 'proof' in content.lower() or 'theorem' in content.lower()
        
        strengths = []
        concerns = []
        
        if has_equations:
            strengths.append("Includes mathematical formalization")
        else:
            concerns.append("Lacks mathematical precision")
        
        if has_proofs:
            strengths.append("Provides rigorous proofs")
        
        recommendation = ReviewRecommendation.ACCEPT if len(strengths) >= 2 else \
                        ReviewRecommendation.ACCEPT_WITH_REVISIONS
        
        return {
            'proposal_id': proposal['id'],
            'reviewer': f"{self.lct.name} (LCT: {self.lct.id[:8]})",
            'action': ReviewAction.ACCEPT_WITH_REVISIONS.value,
            'strengths': strengths,
            'concerns': concerns,
            'recommendation': recommendation.value,
            'timestamp': datetime.now().isoformat()
        }
    
    def _needs_mathematical_counter(self, proposal: Dict) -> bool:
        """Determine if mathematical counter-proposal needed"""
        return 'conceptual' in proposal.get('title', '').lower() and \
               'equation' not in proposal.get('content_summary', '').lower()
    
    def _generate_counter_proposal(self, held_proposal: Dict, section: str) -> Optional[Dict]:
        """Generate a mathematical counter-proposal"""
        return self.governance.create_proposal(
            section_path=section,
            author=f"{self.lct.name} (LCT: {self.lct.id[:8]})",
            title=f"Counter: Mathematical Framework for {held_proposal.get('title', 'Unknown')}",
            proposal_type=ProposalType.EXPANSION,
            content="Rigorous mathematical formulation with proofs",
            rationale="Provides testable predictions and formal analysis",
            specific_changes="Add mathematical appendix with derivations"
        )

class HumanAPI(ParticipantAPI):
    """Human participant API (email/manual interface)"""
    
    def generate_proposals(self, 
                          sections: List[str], 
                          max_proposals: int,
                          context: Optional[Dict] = None) -> List[Dict]:
        """Human proposals via email or manual entry"""
        
        # For testing, simulate human not proposing this cycle
        print(f"[Human API] Notification sent to {self.lct.access_endpoint}")
        print(f"[Human API] Sections available: {sections}")
        print(f"[Human API] Max proposals: {max_proposals}")
        
        # Human might not respond immediately
        return []
    
    def review_proposals(self,
                        proposals: List[Dict],
                        context: Optional[Dict] = None) -> List[Dict]:
        """Human reviews via email or manual entry"""
        
        print(f"[Human API] Review request sent to {self.lct.access_endpoint}")
        print(f"[Human API] {len(proposals)} proposals to review")
        
        # For testing, simulate one review
        if proposals:
            return [{
                'proposal_id': proposals[0]['id'],
                'reviewer': f"{self.lct.name} (LCT: {self.lct.id[:8]})",
                'action': ReviewAction.ACCEPT.value,
                'strengths': ["Well reasoned", "Addresses key issues"],
                'concerns': [],
                'recommendation': ReviewRecommendation.ACCEPT.value,
                'timestamp': datetime.now().isoformat()
            }]
        
        return []

class DefaultAPI(ParticipantAPI):
    """Default API implementation for unknown participant types"""
    
    def generate_proposals(self, 
                          sections: List[str], 
                          max_proposals: int,
                          context: Optional[Dict] = None) -> List[Dict]:
        """Default: no proposals"""
        return []
    
    def review_proposals(self,
                        proposals: List[Dict],
                        context: Optional[Dict] = None) -> List[Dict]:
        """Default: no reviews"""
        return []

class ParticipantAPIFactory:
    """Factory for creating appropriate API instances"""
    
    @staticmethod
    def create_api(lct: LCT) -> ParticipantAPI:
        """Create API instance based on participant type"""
        
        if lct.type == ParticipantType.AI_CLAUDE:
            return ClaudeAPI(lct)
        elif lct.type == ParticipantType.AI_GPT:
            return GPTAPI(lct)
        elif lct.type == ParticipantType.HUMAN:
            return HumanAPI(lct)
        else:
            # Default implementation for unknown types
            return DefaultAPI(lct)

def test_api_limits():
    """Test that APIs respect max_proposals limit"""
    
    from lct_participants import LCT, ParticipantType, AccessMethod
    
    # Create test LCTs
    claude_lct = LCT.create(
        name="Claude-Test",
        participant_type=ParticipantType.AI_CLAUDE,
        access_method=AccessMethod.API,
        access_endpoint="test_endpoint"
    )
    
    gpt_lct = LCT.create(
        name="GPT-Test",
        participant_type=ParticipantType.AI_GPT,
        access_method=AccessMethod.API,
        access_endpoint="test_endpoint"
    )
    
    # Create APIs
    claude_api = ClaudeAPI(claude_lct)
    gpt_api = GPTAPI(gpt_lct)
    
    sections = [
        "04-fundamental-concepts/01-universe-grid",
        "04-fundamental-concepts/02-time-slices",
        "04-fundamental-concepts/03-intent-transfer",
        "04-fundamental-concepts/04-emergence"
    ]
    
    print("Testing API proposal limits:")
    print("-" * 40)
    
    # Test with different limits
    for max_proposals in [0, 1, 2, 5]:
        print(f"\nMax proposals: {max_proposals}")
        
        claude_proposals = claude_api.generate_proposals(sections, max_proposals)
        print(f"  Claude generated: {len(claude_proposals)} proposals")
        assert len(claude_proposals) <= max_proposals, f"Claude exceeded limit: {len(claude_proposals)} > {max_proposals}"
        
        gpt_proposals = gpt_api.generate_proposals(sections, max_proposals)
        print(f"  GPT generated: {len(gpt_proposals)} proposals")
        assert len(gpt_proposals) <= max_proposals, f"GPT exceeded limit: {len(gpt_proposals)} > {max_proposals}"
    
    print("\n✅ All APIs respect max_proposals limit")
    
    # Test unlimited reviews
    print("\nTesting unlimited reviews:")
    test_proposals = [
        {'id': f'test-{i}', 'author': 'Other', 'title': f'Proposal {i}', 'content_summary': 'Test content'}
        for i in range(10)
    ]
    
    claude_reviews = claude_api.review_proposals(test_proposals)
    print(f"  Claude submitted {len(claude_reviews)} reviews for {len(test_proposals)} proposals")
    
    gpt_reviews = gpt_api.review_proposals(test_proposals)
    print(f"  GPT submitted {len(gpt_reviews)} reviews for {len(test_proposals)} proposals")
    
    print("\n✅ APIs can submit unlimited reviews")

if __name__ == "__main__":
    test_api_limits()