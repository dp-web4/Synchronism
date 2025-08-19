#!/usr/bin/env python3
"""
Arbiter System for Whitepaper Governance
Manages arbiter selection and decision-making process
"""

import json
import os
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from enum import Enum
from whitepaper_governance import (
    WhitepaperGovernance,
    ProposalStatus,
    ReviewRecommendation
)

class ArbiterRole(Enum):
    HUMAN = "human"
    AI_CONSENSUS = "ai_consensus"
    TOKEN_WEIGHTED = "token_weighted"
    ROTATING = "rotating"

class Arbiter:
    """Represents an arbiter in the governance system"""
    
    def __init__(self, name: str, role: ArbiterRole, token_balance: Dict[str, int] = None):
        self.name = name
        self.role = role
        self.token_balance = token_balance or {}
        self.decisions_made = []
        self.governance = WhitepaperGovernance()
        
    def evaluate_proposal_with_reviews(self, proposal: Dict, reviews: List[Dict]) -> Tuple[ProposalStatus, str]:
        """Evaluate a proposal considering all reviews"""
        
        # Count review recommendations
        recommendations = {
            "accept": 0,
            "accept_with_revisions": 0,
            "revise_and_resubmit": 0,
            "reject": 0
        }
        
        for review in reviews:
            rec = review.get("recommendation", "").lower()
            if "accept_with_revisions" in rec:
                recommendations["accept_with_revisions"] += 1
            elif "accept" in rec:
                recommendations["accept"] += 1
            elif "revise" in rec:
                recommendations["revise_and_resubmit"] += 1
            elif "reject" in rec:
                recommendations["reject"] += 1
        
        # Decision logic based on consensus
        total_reviews = sum(recommendations.values())
        
        if total_reviews == 0:
            return ProposalStatus.SUBMITTED, "No reviews to evaluate"
        
        # If majority accept (with or without revisions)
        accept_votes = recommendations["accept"] + recommendations["accept_with_revisions"]
        if accept_votes > total_reviews / 2:
            if recommendations["accept_with_revisions"] > recommendations["accept"]:
                return ProposalStatus.ACCEPTED, "Accepted with suggested revisions from reviewers"
            else:
                return ProposalStatus.ACCEPTED, "Accepted based on reviewer consensus"
        
        # If majority reject
        if recommendations["reject"] > total_reviews / 2:
            return ProposalStatus.REJECTED, "Rejected based on reviewer consensus"
        
        # Otherwise needs revision
        return ProposalStatus.SUBMITTED, "Requires revision and resubmission based on reviews"
    
    def make_decision(self, 
                      section_path: str,
                      proposal_id: str,
                      proposal: Dict,
                      reviews: List[Dict],
                      override: Optional[ProposalStatus] = None,
                      override_rationale: Optional[str] = None) -> Dict:
        """Make final arbiter decision on a proposal"""
        
        # Get automated evaluation
        auto_decision, auto_rationale = self.evaluate_proposal_with_reviews(proposal, reviews)
        
        # Apply override if provided (human arbiter privilege)
        if override and self.role == ArbiterRole.HUMAN:
            final_decision = override
            final_rationale = override_rationale or f"Arbiter override: {override.value}"
        else:
            final_decision = auto_decision
            final_rationale = auto_rationale
        
        # Record decision
        decision_data = self.governance.make_arbiter_decision(
            section_path=section_path,
            proposal_id=proposal_id,
            arbiter=self.name,
            decision=final_decision,
            rationale=final_rationale,
            implementation_notes=self._generate_implementation_notes(proposal, reviews) if final_decision == ProposalStatus.ACCEPTED else None
        )
        
        # Log decision
        self.decisions_made.append({
            "timestamp": datetime.now().isoformat(),
            "proposal_id": proposal_id,
            "decision": final_decision.value,
            "rationale": final_rationale
        })
        
        return decision_data
    
    def _generate_implementation_notes(self, proposal: Dict, reviews: List[Dict]) -> str:
        """Generate implementation notes from accepted proposal and reviews"""
        
        notes = ["Implementation Notes:"]
        notes.append(f"- Original proposal: {proposal.get('title', 'Unknown')}")
        
        # Collect all suggested revisions
        all_revisions = {}
        for review in reviews:
            revisions = review.get("suggested_revisions", {})
            for key, value in revisions.items():
                if key not in all_revisions:
                    all_revisions[key] = []
                all_revisions[key].append(value)
        
        if all_revisions:
            notes.append("- Consolidated revisions from reviews:")
            for key, values in all_revisions.items():
                notes.append(f"  - {key}: {'; '.join(set(values))}")
        
        return "\n".join(notes)
    
    def implement_and_log_changes(self, proposal: Dict, section_path: str, arbiter_name: str) -> bool:
        """Implement accepted proposal and log if noteworthy"""
        from changelog_manager import ChangelogManager
        
        # Initialize changelog manager
        changelog_mgr = ChangelogManager(str(self.governance.base_path))
        
        # Here the arbiter would make actual changes to fractal files
        # For now, we'll simulate by checking if changes would be made
        
        # Process and log only if noteworthy changes occurred
        change_entry = changelog_mgr.process_implemented_proposal(
            proposal_id=proposal['id'],
            proposal_data=proposal,
            arbiter=arbiter_name,
            section_path=section_path
        )
        
        if change_entry:
            print(f"✅ Noteworthy change logged: {change_entry['title']}")
            return True
        else:
            # No actual changes to fractal files - routine activity
            return False

class ArbiterSelection:
    """Manages arbiter selection based on governance tokens"""
    
    def __init__(self, governance_path: Path):
        self.governance_path = governance_path
        self.config_path = governance_path / "config"
        self.token_file = self.config_path / "tokens.json"
        
    def get_current_arbiter(self) -> Arbiter:
        """Get the current arbiter based on selection mechanism"""
        
        # For now, use rotating system
        arbiters = self._get_registered_arbiters()
        if not arbiters:
            # Default human arbiter
            return Arbiter("Dennis", ArbiterRole.HUMAN)
        
        # Rotate based on decision count
        decision_counts = {}
        for arbiter in arbiters:
            decision_counts[arbiter.name] = len(arbiter.decisions_made)
        
        # Select arbiter with fewest decisions
        selected = min(arbiters, key=lambda a: len(a.decisions_made))
        return selected
    
    def _get_registered_arbiters(self) -> List[Arbiter]:
        """Get list of registered arbiters"""
        
        # Load from configuration
        arbiters = []
        
        # Always include human arbiter
        arbiters.append(Arbiter("Dennis", ArbiterRole.HUMAN))
        
        # Add AI consensus arbiter
        arbiters.append(Arbiter("AI-Consensus", ArbiterRole.AI_CONSENSUS))
        
        # Add token-weighted arbiter if tokens exist
        if self.token_file.exists():
            with open(self.token_file, 'r') as f:
                token_data = json.load(f)
            
            contributors = token_data.get("contributors", {})
            for contributor, tokens in contributors.items():
                if sum(tokens.values()) > 10:  # Threshold for arbiter eligibility
                    arbiters.append(Arbiter(
                        contributor,
                        ArbiterRole.TOKEN_WEIGHTED,
                        tokens
                    ))
        
        return arbiters
    
    def calculate_voting_power(self, arbiter: Arbiter, scale: str) -> float:
        """Calculate arbiter's voting power for a specific scale"""
        
        if arbiter.role == ArbiterRole.HUMAN:
            return 1.0  # Human always has full vote
        
        if arbiter.role == ArbiterRole.AI_CONSENSUS:
            return 0.8  # AI consensus has 80% weight
        
        if arbiter.role == ArbiterRole.TOKEN_WEIGHTED:
            # Voting power based on tokens for that scale
            scale_tokens = arbiter.token_balance.get(scale, 0)
            total_tokens = 100  # Default total supply per scale
            return min(scale_tokens / total_tokens, 0.5)  # Cap at 50%
        
        return 0.1  # Default minimal voting power

class GovernanceArbiterSystem:
    """Complete arbiter system for governance"""
    
    def __init__(self, base_path: str = "/mnt/c/exe/projects/ai-agents/Synchronism"):
        self.base_path = Path(base_path)
        self.governance = WhitepaperGovernance()
        self.selection = ArbiterSelection(self.base_path / "scripts" / "governance")
        
    def process_pending_proposals(self) -> List[Dict]:
        """Process all proposals pending arbitration"""
        
        results = []
        
        # Get pending proposals
        config_file = self.governance.config_path / "whitepaper_proposals.json"
        if not config_file.exists():
            return results
        
        with open(config_file, 'r') as f:
            data = json.load(f)
        
        for proposal in data.get("proposals", []):
            if proposal["status"] == ProposalStatus.UNDER_REVIEW.value:
                # Get current arbiter
                arbiter = self.selection.get_current_arbiter()
                
                # Get reviews for this proposal
                reviews = proposal.get("reviews", [])
                
                # Make decision
                decision = arbiter.make_decision(
                    section_path=proposal["section"],
                    proposal_id=proposal["id"],
                    proposal=proposal,
                    reviews=reviews
                )
                
                results.append({
                    "proposal_id": proposal["id"],
                    "arbiter": arbiter.name,
                    "decision": decision
                })
        
        return results
    
    def implement_accepted_proposals(self, section_path: str) -> List[str]:
        """Implement all accepted proposals for a section"""
        
        implemented = []
        
        # Read proposals
        config_file = self.governance.config_path / "whitepaper_proposals.json"
        if not config_file.exists():
            return implemented
        
        with open(config_file, 'r') as f:
            data = json.load(f)
        
        for proposal in data.get("proposals", []):
            if (proposal["section"] == section_path and 
                proposal["status"] == ProposalStatus.ACCEPTED.value):
                
                # Mark as implemented
                self.governance._update_proposal_status(
                    section_path,
                    proposal["id"],
                    ProposalStatus.IMPLEMENTED
                )
                
                implemented.append(proposal["id"])
                
                # Log noteworthy changes if any
                self.implement_and_log_changes(
                    proposal,
                    section_path,
                    "Arbiter"
                )
        
        return implemented
    
    def trigger_whitepaper_rebuild(self) -> Dict:
        """Trigger whitepaper rebuild after implementations"""
        from whitepaper_builder import WhitepaperBuilder
        
        builder = WhitepaperBuilder(str(self.governance.base_path))
        result = builder.rebuild_if_changed()
        
        if result:
            print("\n" + "="*60)
            print(" Whitepaper Rebuild Triggered")
            print("="*60)
            
            for fmt, (success, msg) in result.get('formats', {}).items():
                if success:
                    print(f"  ✅ {fmt}: {msg}")
                else:
                    print(f"  ❌ {fmt}: {msg}")
        
        return result if result else {'rebuilt': False}

def main():
    """Example arbiter workflow"""
    system = GovernanceArbiterSystem()
    
    # Process pending proposals
    decisions = system.process_pending_proposals()
    
    if decisions:
        print(f"Processed {len(decisions)} proposals:")
        for decision in decisions:
            print(f"  - Proposal {decision['proposal_id']}: decided by {decision['arbiter']}")
    else:
        print("No proposals pending arbitration")
    
    # Implement accepted proposals
    sections = [
        "04-fundamental-concepts/01-universe-grid",
        "04-fundamental-concepts/02-time-slices"
    ]
    
    for section in sections:
        implemented = system.implement_accepted_proposals(section)
        if implemented:
            print(f"Implemented {len(implemented)} proposals in {section}")

if __name__ == "__main__":
    main()