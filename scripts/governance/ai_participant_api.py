#!/usr/bin/env python3
"""
AI Participant API for Whitepaper Governance
Provides interfaces for Claude and GPT to participate in governance
"""

import json
import os
import hashlib
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from whitepaper_governance import (
    WhitepaperGovernance, 
    ProposalType, 
    ProposalStatus,
    ReviewRecommendation
)

class AIParticipant:
    """Base class for AI participants in governance"""
    
    def __init__(self, name: str, model_id: str, api_key: Optional[str] = None):
        self.name = name
        self.model_id = model_id
        self.api_key = api_key
        self.governance = WhitepaperGovernance()
        self.participation_log = []
        
    def analyze_section(self, section_path: str) -> Dict:
        """Analyze a section for potential improvements"""
        section_dir = self.governance.whitepaper_path / section_path
        
        # Read the current content
        content_files = list(section_dir.glob("*.md"))
        if not content_files:
            return {"error": "No content files found"}
        
        analysis = {
            "section": section_path,
            "analyzer": self.name,
            "timestamp": datetime.now().isoformat(),
            "findings": []
        }
        
        # Check for existing proposals and reviews
        meta_dir = section_dir / "meta"
        if meta_dir.exists():
            proposals = list((meta_dir / "proposals").glob("*.md")) if (meta_dir / "proposals").exists() else []
            reviews = list((meta_dir / "reviews").glob("*.md")) if (meta_dir / "reviews").exists() else []
            
            analysis["existing_proposals"] = len(proposals)
            analysis["existing_reviews"] = len(reviews)
        
        # Analyze future considerations
        future_file = meta_dir / "future-considerations.md" if meta_dir.exists() else None
        if future_file and future_file.exists():
            content = future_file.read_text()
            # Parse active considerations
            if "## Active Considerations" in content:
                active_section = content.split("## Active Considerations")[1].split("##")[0]
                considerations = [line.strip() for line in active_section.split("\n") if line.startswith("###")]
                analysis["active_considerations"] = considerations
        
        return analysis
    
    def create_proposal(self, 
                       section_path: str,
                       title: str,
                       proposal_type: ProposalType,
                       content: str,
                       rationale: str,
                       specific_changes: str) -> Dict:
        """Create a proposal for section improvement"""
        
        proposal = self.governance.create_proposal(
            section_path=section_path,
            author=f"{self.name} ({self.model_id})",
            title=title,
            proposal_type=proposal_type,
            content=content,
            rationale=rationale,
            specific_changes=specific_changes
        )
        
        # Log participation
        self.participation_log.append({
            "action": "create_proposal",
            "timestamp": datetime.now().isoformat(),
            "details": proposal
        })
        
        return proposal
    
    def review_proposal(self,
                       section_path: str,
                       proposal_id: str,
                       recommendation: ReviewRecommendation,
                       analysis: Dict) -> Dict:
        """Review another participant's proposal"""
        
        review = self.governance.submit_review(
            section_path=section_path,
            proposal_id=proposal_id,
            reviewer=f"{self.name} ({self.model_id})",
            recommendation=recommendation,
            strengths=analysis.get("strengths", []),
            concerns=analysis.get("concerns", []),
            suggested_revisions=analysis.get("revisions", {})
        )
        
        # Log participation
        self.participation_log.append({
            "action": "submit_review",
            "timestamp": datetime.now().isoformat(),
            "details": review
        })
        
        return review

class ClaudeParticipant(AIParticipant):
    """Claude-specific participant implementation"""
    
    def __init__(self):
        super().__init__(
            name="Claude",
            model_id="claude-3.5-sonnet"
        )
        self.focus_areas = ["conceptual_clarity", "philosophical_coherence", "cross_scale_consistency"]
    
    def generate_proposal_from_analysis(self, section_path: str, analysis: Dict) -> Optional[Dict]:
        """Generate a proposal based on section analysis"""
        
        # Claude focuses on conceptual and philosophical improvements
        if "active_considerations" in analysis:
            for consideration in analysis["active_considerations"]:
                if "conceptual" in consideration.lower() or "philosophical" in consideration.lower():
                    # Generate proposal for this consideration
                    return self.create_proposal(
                        section_path=section_path,
                        title=f"Conceptual Enhancement: {consideration.replace('###', '').strip()}",
                        proposal_type=ProposalType.CLARIFICATION,
                        content="Enhance conceptual clarity and philosophical alignment",
                        rationale="Maintains coherence with Synchronism's philosophical framework",
                        specific_changes="To be detailed based on specific consideration"
                    )
        return None
    
    def evaluate_proposal(self, proposal: Dict) -> Tuple[ReviewRecommendation, Dict]:
        """Evaluate a proposal from Claude's perspective"""
        
        analysis = {
            "strengths": [],
            "concerns": [],
            "revisions": {}
        }
        
        # Check for philosophical coherence
        if "fractal" in proposal.get("content_summary", "").lower():
            analysis["strengths"].append("Maintains fractal self-similarity principle")
        
        if "markov" in proposal.get("content_summary", "").lower():
            analysis["strengths"].append("Preserves Markov blanket boundaries")
        
        # Default to constructive feedback
        if not analysis["concerns"]:
            analysis["concerns"].append("Consider adding more concrete examples")
        
        # Determine recommendation
        if len(analysis["strengths"]) > len(analysis["concerns"]):
            recommendation = ReviewRecommendation.ACCEPT_WITH_REVISIONS
        else:
            recommendation = ReviewRecommendation.REVISE_AND_RESUBMIT
        
        return recommendation, analysis

class GPTParticipant(AIParticipant):
    """GPT-specific participant implementation"""
    
    def __init__(self):
        super().__init__(
            name="GPT-4",
            model_id="gpt-4-turbo"
        )
        self.focus_areas = ["mathematical_rigor", "scientific_accuracy", "technical_precision"]
    
    def generate_proposal_from_analysis(self, section_path: str, analysis: Dict) -> Optional[Dict]:
        """Generate a proposal based on section analysis"""
        
        # GPT focuses on mathematical and technical improvements
        if "active_considerations" in analysis:
            for consideration in analysis["active_considerations"]:
                if "mathematical" in consideration.lower() or "technical" in consideration.lower():
                    # Generate proposal for this consideration
                    return self.create_proposal(
                        section_path=section_path,
                        title=f"Mathematical Formalization: {consideration.replace('###', '').strip()}",
                        proposal_type=ProposalType.EXPANSION,
                        content="Add rigorous mathematical framework",
                        rationale="Provides formal foundation for theoretical concepts",
                        specific_changes="Mathematical equations and proofs to be added"
                    )
        return None
    
    def evaluate_proposal(self, proposal: Dict) -> Tuple[ReviewRecommendation, Dict]:
        """Evaluate a proposal from GPT's perspective"""
        
        analysis = {
            "strengths": [],
            "concerns": [],
            "revisions": {}
        }
        
        # Check for mathematical rigor
        if "equation" in proposal.get("content_summary", "").lower() or "theorem" in proposal.get("content_summary", "").lower():
            analysis["strengths"].append("Includes mathematical formalization")
        else:
            analysis["concerns"].append("Lacks mathematical precision")
        
        # Check for testability
        if "prediction" in proposal.get("content_summary", "").lower() or "observable" in proposal.get("content_summary", "").lower():
            analysis["strengths"].append("Includes testable predictions")
        
        # Suggest improvements
        if "tensor" in proposal.get("content_summary", "").lower():
            analysis["revisions"]["terminology"] = "Ensure tensor notation follows standard conventions"
        
        # Determine recommendation
        if len(analysis["strengths"]) >= 2:
            recommendation = ReviewRecommendation.ACCEPT
        elif len(analysis["concerns"]) > 2:
            recommendation = ReviewRecommendation.REJECT
        else:
            recommendation = ReviewRecommendation.ACCEPT_WITH_REVISIONS
        
        return recommendation, analysis

class GovernanceOrchestrator:
    """Orchestrates the governance process with AI participants"""
    
    def __init__(self):
        self.claude = ClaudeParticipant()
        self.gpt = GPTParticipant()
        self.governance = WhitepaperGovernance()
        self.round_robin_order = [self.claude, self.gpt]
        self.current_round = 0
        
    def run_governance_round(self, section_paths: List[str]):
        """Run a complete governance round"""
        
        results = {
            "round": self.current_round,
            "timestamp": datetime.now().isoformat(),
            "proposals_created": [],
            "reviews_submitted": [],
            "sections_processed": []
        }
        
        for section_path in section_paths:
            # Current proposer in round-robin
            proposer = self.round_robin_order[self.current_round % len(self.round_robin_order)]
            reviewer = self.round_robin_order[(self.current_round + 1) % len(self.round_robin_order)]
            
            # Analyze section
            analysis = proposer.analyze_section(section_path)
            
            # Generate proposal if applicable
            proposal = proposer.generate_proposal_from_analysis(section_path, analysis)
            if proposal:
                results["proposals_created"].append(proposal)
                
                # Have the other AI review it
                recommendation, review_analysis = reviewer.evaluate_proposal(proposal)
                review = reviewer.review_proposal(
                    section_path=section_path,
                    proposal_id=proposal["id"],
                    recommendation=recommendation,
                    analysis=review_analysis
                )
                results["reviews_submitted"].append(review)
            
            results["sections_processed"].append(section_path)
        
        self.current_round += 1
        return results
    
    def get_pending_arbitration(self) -> List[Dict]:
        """Get proposals that need arbiter decision"""
        
        config_file = self.governance.config_path / "whitepaper_proposals.json"
        if not config_file.exists():
            return []
        
        with open(config_file, 'r') as f:
            data = json.load(f)
        
        pending = []
        for proposal in data.get("proposals", []):
            if proposal["status"] == ProposalStatus.UNDER_REVIEW.value:
                # Check if it has reviews
                if proposal.get("reviews", []):
                    pending.append(proposal)
        
        return pending

def main():
    """Example orchestration"""
    orchestrator = GovernanceOrchestrator()
    
    # Run a governance round on fundamental concepts
    sections = [
        "04-fundamental-concepts/01-universe-grid",
        "04-fundamental-concepts/02-time-slices",
        "04-fundamental-concepts/03-intent-transfer"
    ]
    
    results = orchestrator.run_governance_round(sections)
    print(f"Governance round {results['round']} completed")
    print(f"Proposals created: {len(results['proposals_created'])}")
    print(f"Reviews submitted: {len(results['reviews_submitted'])}")
    
    # Check for pending arbitration
    pending = orchestrator.get_pending_arbitration()
    if pending:
        print(f"\nProposals awaiting arbitration: {len(pending)}")
        for proposal in pending:
            print(f"  - {proposal['id']}: {proposal['title']}")

if __name__ == "__main__":
    main()