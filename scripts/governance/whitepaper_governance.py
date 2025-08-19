#!/usr/bin/env python3
"""
Whitepaper Governance System
Manages proposals, reviews, and updates to the Synchronism whitepaper sections
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
from enum import Enum

class ProposalStatus(Enum):
    DRAFT = "draft"
    SUBMITTED = "submitted"
    UNDER_REVIEW = "under_review"
    ACCEPTED = "accepted"
    REJECTED = "rejected"
    IMPLEMENTED = "implemented"

class ProposalType(Enum):
    INITIAL = "INITIAL"
    UPDATE = "UPDATE"
    CLARIFICATION = "CLARIFICATION"
    EXPANSION = "EXPANSION"
    CORRECTION = "CORRECTION"

class ReviewRecommendation(Enum):
    ACCEPT = "ACCEPT"
    ACCEPT_WITH_REVISIONS = "ACCEPT_WITH_REVISIONS"
    REVISE_AND_RESUBMIT = "REVISE_AND_RESUBMIT"
    REJECT = "REJECT"

class WhitepaperGovernance:
    def __init__(self, base_path: str = "/mnt/c/exe/projects/ai-agents/Synchronism"):
        self.base_path = Path(base_path)
        self.whitepaper_path = self.base_path / "whitepaper" / "sections"
        self.governance_path = self.base_path / "scripts" / "governance"
        self.config_path = self.governance_path / "config"
        
    def create_proposal(self, 
                       section_path: str,
                       author: str,
                       title: str,
                       proposal_type: ProposalType,
                       content: str,
                       rationale: str,
                       specific_changes: str) -> Dict:
        """Create a new proposal for a whitepaper section"""
        
        section_dir = self.whitepaper_path / section_path
        meta_dir = section_dir / "meta"
        proposals_dir = meta_dir / "proposals"
        
        # Ensure directories exist
        proposals_dir.mkdir(parents=True, exist_ok=True)
        
        # Get next proposal ID
        existing_proposals = list(proposals_dir.glob("*.md"))
        next_id = len(existing_proposals) + 1
        proposal_id = f"{next_id:03d}"
        
        # Create proposal metadata
        proposal_data = {
            "id": proposal_id,
            "author": author,
            "title": title,
            "date": datetime.now().isoformat(),
            "status": ProposalStatus.SUBMITTED.value,
            "type": proposal_type.value,
            "section": section_path,
            "content_summary": content[:200] + "..." if len(content) > 200 else content,
            "reviews": []
        }
        
        # Create proposal markdown
        proposal_md = f"""# Proposal {proposal_id}: {title}

## Metadata
- **ID**: {proposal_id}
- **Author**: {author}
- **Date**: {datetime.now().date()}
- **Status**: {ProposalStatus.SUBMITTED.value}
- **Type**: {proposal_type.value}
- **Section**: {section_path}

## Proposed Change
{content}

## Rationale
{rationale}

## Specific Text Changes
{specific_changes}

## Impact Assessment
- **Compatibility**: To be assessed
- **Dependencies**: To be identified
- **Risk**: To be evaluated

## Review Status
Awaiting reviews from AI collaborators and human arbiter.
"""
        
        # Save proposal
        proposal_file = proposals_dir / f"{proposal_id}-{title.lower().replace(' ', '-')}.md"
        proposal_file.write_text(proposal_md)
        
        # Update proposals index
        self._update_proposals_index(section_path, proposal_data)
        
        return proposal_data
    
    def submit_review(self,
                      section_path: str,
                      proposal_id: str,
                      reviewer: str,
                      recommendation: ReviewRecommendation,
                      strengths: List[str],
                      concerns: List[str],
                      suggested_revisions: Dict[str, str]) -> Dict:
        """Submit a review for a proposal"""
        
        section_dir = self.whitepaper_path / section_path
        reviews_dir = section_dir / "meta" / "reviews"
        reviews_dir.mkdir(parents=True, exist_ok=True)
        
        review_data = {
            "proposal_id": proposal_id,
            "reviewer": reviewer,
            "date": datetime.now().isoformat(),
            "recommendation": recommendation.value,
            "strengths": strengths,
            "concerns": concerns,
            "suggested_revisions": suggested_revisions
        }
        
        # Create review markdown
        review_md = f"""# Review of Proposal {proposal_id}

## Reviewer Information
- **Reviewer**: {reviewer}
- **Date**: {datetime.now().date()}
- **Review ID**: {proposal_id}-{reviewer}

## Summary Recommendation
**{recommendation.value}**

## Detailed Review

### Strengths
{chr(10).join(f"{i+1}. {s}" for i, s in enumerate(strengths))}

### Concerns
{chr(10).join(f"{i+1}. {c}" for i, c in enumerate(concerns))}

### Suggested Revisions
{chr(10).join(f"- {k}: {v}" for k, v in suggested_revisions.items())}

## Final Assessment
Review submitted via governance API.
"""
        
        # Save review
        review_file = reviews_dir / f"{proposal_id}-review-{reviewer.lower().replace(' ', '-')}.md"
        review_file.write_text(review_md)
        
        # Update proposal status
        self._update_proposal_status(section_path, proposal_id, ProposalStatus.UNDER_REVIEW)
        
        return review_data
    
    def make_arbiter_decision(self,
                              section_path: str,
                              proposal_id: str,
                              arbiter: str,
                              decision: ProposalStatus,
                              rationale: str,
                              implementation_notes: Optional[str] = None) -> Dict:
        """Record arbiter's final decision on a proposal"""
        
        decision_data = {
            "proposal_id": proposal_id,
            "arbiter": arbiter,
            "date": datetime.now().isoformat(),
            "decision": decision.value,
            "rationale": rationale,
            "implementation_notes": implementation_notes
        }
        
        # Update proposal status
        self._update_proposal_status(section_path, proposal_id, decision)
        
        # If accepted, update changelog
        if decision == ProposalStatus.ACCEPTED:
            self._update_changelog(section_path, proposal_id, arbiter)
            
        # If implemented, update the actual section file
        if decision == ProposalStatus.IMPLEMENTED and implementation_notes:
            self._implement_changes(section_path, implementation_notes)
        
        return decision_data
    
    def _update_proposals_index(self, section_path: str, proposal_data: Dict):
        """Update the proposals index file"""
        index_file = self.config_path / "whitepaper_proposals.json"
        
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)
        else:
            index = {"proposals": []}
        
        index["proposals"].append(proposal_data)
        index["last_updated"] = datetime.now().isoformat()
        
        with open(index_file, 'w') as f:
            json.dump(index, f, indent=2)
    
    def _update_proposal_status(self, section_path: str, proposal_id: str, status: ProposalStatus):
        """Update the status of a proposal"""
        index_file = self.config_path / "whitepaper_proposals.json"
        
        if index_file.exists():
            with open(index_file, 'r') as f:
                index = json.load(f)
            
            for proposal in index["proposals"]:
                if proposal["id"] == proposal_id and proposal["section"] == section_path:
                    proposal["status"] = status.value
                    proposal["last_updated"] = datetime.now().isoformat()
                    break
            
            with open(index_file, 'w') as f:
                json.dump(index, f, indent=2)
    
    def _update_changelog(self, section_path: str, proposal_id: str, arbiter: str):
        """Add entry to section's changelog"""
        section_dir = self.whitepaper_path / section_path
        changelog_file = section_dir / "meta" / "changelog.md"
        
        if changelog_file.exists():
            content = changelog_file.read_text()
        else:
            content = "# Changelog\n\n"
        
        entry = f"\n### {datetime.now().date()} [{arbiter}] [UPDATE]\n- Accepted proposal {proposal_id}\n"
        
        content += entry
        changelog_file.write_text(content)
    
    def _implement_changes(self, section_path: str, changes: str):
        """Apply approved changes to the section file"""
        # This is where the arbiter's actual edits would be applied
        # For now, we'll just log that changes were made
        print(f"Implementing changes to {section_path}")
        print(f"Changes: {changes}")

def main():
    """Example usage of the governance system"""
    gov = WhitepaperGovernance()
    
    # Example: Create a proposal
    proposal = gov.create_proposal(
        section_path="04-fundamental-concepts/01-universe-grid",
        author="Claude-3.5",
        title="Grid Topology Enhancement",
        proposal_type=ProposalType.EXPANSION,
        content="Add non-Euclidean geometry and fractal structure",
        rationale="Aligns with general relativity and scale invariance",
        specific_changes="Add new subsection on fractal structure"
    )
    print(f"Created proposal: {proposal['id']}")
    
    # Example: Submit a review
    review = gov.submit_review(
        section_path="04-fundamental-concepts/01-universe-grid",
        proposal_id=proposal['id'],
        reviewer="GPT-4",
        recommendation=ReviewRecommendation.ACCEPT_WITH_REVISIONS,
        strengths=["Mathematical rigor", "Consistency with framework"],
        concerns=["Needs clearer quantum connection"],
        suggested_revisions={"Line 3": "Use 'metric tensor' instead of 'warping'"}
    )
    print(f"Submitted review from {review['reviewer']}")

if __name__ == "__main__":
    main()