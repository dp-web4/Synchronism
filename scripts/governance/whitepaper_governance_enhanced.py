#!/usr/bin/env python3
"""
Enhanced Whitepaper Governance System
Adds withdrawal, expiration, and duplicate detection to proposal management
"""

import json
import os
import hashlib
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from enum import Enum
from difflib import SequenceMatcher

# Import original classes
from whitepaper_governance import (
    ProposalStatus as BaseProposalStatus,
    ProposalType,
    ReviewRecommendation,
    WhitepaperGovernance as BaseWhitepaperGovernance
)

class ProposalStatus(Enum):
    """Extended proposal statuses"""
    DRAFT = "draft"
    SUBMITTED = "submitted"
    UNDER_REVIEW = "under_review"
    ACCEPTED = "accepted"
    REJECTED = "rejected"
    IMPLEMENTED = "implemented"
    # New statuses
    WITHDRAWN = "withdrawn"
    EXPIRED = "expired"
    DUPLICATE = "duplicate"
    SUPERSEDED = "superseded"

class WithdrawalReason(Enum):
    """Reasons for withdrawing proposals"""
    AUTHOR_REQUEST = "author_request"
    DUPLICATE_DETECTED = "duplicate_detected"
    SUPERSEDED_BY_NEW = "superseded_by_new"
    EXPIRED_TIMEOUT = "expired_timeout"
    ADMINISTRATIVE = "administrative"

class EnhancedWhitepaperGovernance(BaseWhitepaperGovernance):
    """Enhanced governance with withdrawal and lifecycle management"""
    
    def __init__(self, base_path: str = "/mnt/c/exe/projects/ai-agents/Synchronism"):
        super().__init__(base_path)
        self.expiration_days = 30  # Default expiration period
        self.similarity_threshold = 0.85  # For duplicate detection
        
    def load_proposals(self) -> List[Dict]:
        """Load all proposals from config"""
        proposals_file = self.config_path / "whitepaper_proposals.json"
        if proposals_file.exists():
            with open(proposals_file, 'r') as f:
                data = json.load(f)
                return data.get('proposals', [])
        return []
    
    def save_proposals(self, proposals: List[Dict]):
        """Save proposals to config"""
        proposals_file = self.config_path / "whitepaper_proposals.json"
        data = {
            'proposals': proposals,
            'last_updated': datetime.now().isoformat()
        }
        with open(proposals_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    def withdraw_proposal(self, 
                         proposal_id: str,
                         reason: WithdrawalReason,
                         withdrawer: str,
                         comment: Optional[str] = None) -> Dict:
        """
        Withdraw a proposal
        
        Args:
            proposal_id: The proposal to withdraw
            reason: Reason for withdrawal
            withdrawer: Who is withdrawing (author, reviewer, arbiter)
            comment: Optional explanation
            
        Returns:
            Result dictionary with status and message
        """
        proposals = self.load_proposals()
        
        # Find the proposal
        proposal = None
        for p in proposals:
            if p['id'] == proposal_id:
                proposal = p
                break
        
        if not proposal:
            return {
                'status': 'error',
                'message': f'Proposal {proposal_id} not found'
            }
        
        # Check if already final
        if proposal['status'] in ['implemented', 'withdrawn', 'expired']:
            return {
                'status': 'error',
                'message': f'Cannot withdraw proposal in {proposal["status"]} status'
            }
        
        # Check permissions
        if reason == WithdrawalReason.AUTHOR_REQUEST:
            # Only author can withdraw with this reason
            if withdrawer != proposal['author']:
                return {
                    'status': 'error',
                    'message': 'Only the author can withdraw with AUTHOR_REQUEST'
                }
        
        # Update proposal status
        proposal['status'] = ProposalStatus.WITHDRAWN.value
        proposal['withdrawal'] = {
            'reason': reason.value,
            'withdrawer': withdrawer,
            'date': datetime.now().isoformat(),
            'comment': comment
        }
        proposal['last_updated'] = datetime.now().isoformat()
        
        # Save updated proposals
        self.save_proposals(proposals)
        
        # Archive if needed
        self._consider_archival(proposal)
        
        return {
            'status': 'success',
            'message': f'Proposal {proposal_id} withdrawn',
            'proposal': proposal
        }
    
    def check_duplicates(self, 
                        section: str,
                        title: str,
                        content: str) -> List[Tuple[str, float]]:
        """
        Check for duplicate proposals
        
        Returns:
            List of (proposal_id, similarity_score) for potential duplicates
        """
        proposals = self.load_proposals()
        duplicates = []
        
        # Create comparison string
        new_text = f"{title} {content}".lower()
        
        for proposal in proposals:
            # Only check active proposals in same section
            if (proposal.get('section') == section and 
                proposal.get('status') in ['submitted', 'under_review']):
                
                # Compare text similarity
                existing_text = f"{proposal.get('title', '')} {proposal.get('content_summary', '')}".lower()
                similarity = SequenceMatcher(None, new_text, existing_text).ratio()
                
                if similarity > self.similarity_threshold:
                    duplicates.append((proposal['id'], similarity))
        
        # Sort by similarity score
        duplicates.sort(key=lambda x: x[1], reverse=True)
        return duplicates
    
    def expire_old_proposals(self, dry_run: bool = True) -> List[Dict]:
        """
        Expire proposals that have been inactive too long
        
        Args:
            dry_run: If True, only report what would be expired
            
        Returns:
            List of expired proposals
        """
        proposals = self.load_proposals()
        expired = []
        cutoff_date = datetime.now() - timedelta(days=self.expiration_days)
        
        for proposal in proposals:
            # Only check submitted proposals without reviews
            if (proposal.get('status') == 'submitted' and 
                len(proposal.get('reviews', [])) == 0):
                
                # Check last activity date
                last_activity = proposal.get('last_updated', proposal.get('date'))
                if last_activity:
                    try:
                        activity_date = datetime.fromisoformat(last_activity)
                        if activity_date < cutoff_date:
                            if not dry_run:
                                proposal['status'] = ProposalStatus.EXPIRED.value
                                proposal['expiration'] = {
                                    'date': datetime.now().isoformat(),
                                    'days_inactive': (datetime.now() - activity_date).days
                                }
                            expired.append(proposal)
                    except:
                        pass
        
        if not dry_run and expired:
            self.save_proposals(proposals)
            # Archive expired proposals
            for p in expired:
                self._consider_archival(p)
        
        return expired
    
    def supersede_proposal(self,
                          old_id: str,
                          new_id: str,
                          superseder: str) -> Dict:
        """
        Mark a proposal as superseded by another
        
        Args:
            old_id: Proposal being superseded
            new_id: New proposal that supersedes it
            superseder: Who is marking it superseded
            
        Returns:
            Result dictionary
        """
        proposals = self.load_proposals()
        
        old_proposal = None
        new_proposal = None
        
        for p in proposals:
            if p['id'] == old_id:
                old_proposal = p
            elif p['id'] == new_id:
                new_proposal = p
        
        if not old_proposal:
            return {
                'status': 'error',
                'message': f'Old proposal {old_id} not found'
            }
        
        if not new_proposal:
            return {
                'status': 'error',
                'message': f'New proposal {new_id} not found'
            }
        
        # Update old proposal
        old_proposal['status'] = ProposalStatus.SUPERSEDED.value
        old_proposal['superseded_by'] = {
            'proposal_id': new_id,
            'date': datetime.now().isoformat(),
            'marked_by': superseder
        }
        old_proposal['last_updated'] = datetime.now().isoformat()
        
        # Add reference in new proposal
        if 'supersedes' not in new_proposal:
            new_proposal['supersedes'] = []
        new_proposal['supersedes'].append(old_id)
        new_proposal['last_updated'] = datetime.now().isoformat()
        
        self.save_proposals(proposals)
        
        return {
            'status': 'success',
            'message': f'Proposal {old_id} superseded by {new_id}',
            'old_proposal': old_proposal,
            'new_proposal': new_proposal
        }
    
    def get_active_proposals(self, section: Optional[str] = None) -> List[Dict]:
        """Get only active proposals (not withdrawn/expired/etc)"""
        proposals = self.load_proposals()
        active_statuses = ['submitted', 'under_review', 'accepted']
        
        active = [p for p in proposals 
                 if p.get('status') in active_statuses]
        
        if section:
            active = [p for p in active if p.get('section') == section]
        
        return active
    
    def _consider_archival(self, proposal: Dict):
        """Consider if a proposal should be archived"""
        # Import the cleanup module if available
        try:
            from proposal_cleanup import ProposalCleanup, ArchiveReason
            cleanup = ProposalCleanup()
            
            # Determine archive reason based on status
            reason_map = {
                'withdrawn': ArchiveReason.WITHDRAWN,
                'expired': ArchiveReason.EXPIRED,
                'rejected': ArchiveReason.REJECTED,
                'superseded': ArchiveReason.SUPERSEDED
            }
            
            reason = reason_map.get(proposal.get('status'))
            if reason:
                # Has reviews = historical value
                historical = len(proposal.get('reviews', [])) > 0
                cleanup.archive_proposal(proposal, reason, historical)
                cleanup.save_archive()
        except ImportError:
            # Cleanup module not available, skip archival
            pass
    
    def create_enhanced_proposal(self,
                                section_path: str,
                                author: str,
                                title: str,
                                proposal_type: ProposalType,
                                content: str,
                                rationale: str,
                                specific_changes: str) -> Dict:
        """
        Create proposal with duplicate checking
        
        Returns proposal dict with potential duplicates noted
        """
        # Check for duplicates first
        duplicates = self.check_duplicates(section_path, title, content)
        
        if duplicates:
            # Warn about duplicates but allow creation
            result = super().create_proposal(
                section_path, author, title, proposal_type,
                content, rationale, specific_changes
            )
            result['potential_duplicates'] = duplicates
            result['warning'] = f'Found {len(duplicates)} similar proposals'
            return result
        
        # No duplicates, create normally
        return super().create_proposal(
            section_path, author, title, proposal_type,
            content, rationale, specific_changes
        )

def main():
    """Test the enhanced governance system"""
    governance = EnhancedWhitepaperGovernance()
    
    print("Enhanced Whitepaper Governance System")
    print("-" * 60)
    
    # Show current active proposals
    active = governance.get_active_proposals()
    print(f"Active proposals: {len(active)}")
    for p in active:
        print(f"  [{p['id']}] {p['title'][:40]}... ({p['status']})")
    
    print("\nChecking for expired proposals...")
    expired = governance.expire_old_proposals(dry_run=True)
    if expired:
        print(f"Would expire {len(expired)} proposals:")
        for p in expired:
            print(f"  [{p['id']}] {p['title'][:40]}...")
    else:
        print("No proposals to expire")
    
    print("\nEnhanced governance ready for use!")

if __name__ == "__main__":
    main()