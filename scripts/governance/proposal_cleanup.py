#!/usr/bin/env python3
"""
Proposal Cleanup and Archive Management
Handles cleanup of test proposals and archival of historical proposals
"""

import json
import os
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from enum import Enum

class ArchiveReason(Enum):
    """Reasons for archiving proposals"""
    TEST_CLEANUP = "test_cleanup"
    DUPLICATE = "duplicate"
    REJECTED = "rejected"
    WITHDRAWN = "withdrawn"
    EXPIRED = "expired"
    SUPERSEDED = "superseded"
    INTEGRATED = "integrated"

class ProposalCleanup:
    """Manages proposal cleanup and archival"""
    
    def __init__(self, config_path: str = None):
        if config_path is None:
            # Default to scripts/governance/config
            self.config_path = Path(__file__).parent / "config"
        else:
            self.config_path = Path(config_path)
            
        self.proposals_file = self.config_path / "whitepaper_proposals.json"
        self.archive_file = self.config_path / "whitepaper_proposals_archive.json"
        
        # Load current proposals
        self.load_proposals()
        self.load_archive()
    
    def load_proposals(self):
        """Load current proposals from file"""
        if self.proposals_file.exists():
            with open(self.proposals_file, 'r') as f:
                data = json.load(f)
                self.proposals = data.get('proposals', [])
                self.last_updated = data.get('last_updated', '')
        else:
            self.proposals = []
            self.last_updated = ''
    
    def load_archive(self):
        """Load archived proposals from file"""
        if self.archive_file.exists():
            with open(self.archive_file, 'r') as f:
                self.archive_data = json.load(f)
        else:
            self.archive_data = {
                'archived_proposals': [],
                'archive_policy': {
                    'retention_days': 90,
                    'last_purge': None,
                    'total_archived': 0,
                    'total_purged': 0
                }
            }
    
    def save_proposals(self):
        """Save current proposals to file"""
        data = {
            'proposals': self.proposals,
            'last_updated': datetime.now().isoformat()
        }
        with open(self.proposals_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    def save_archive(self):
        """Save archive to file"""
        with open(self.archive_file, 'w') as f:
            json.dump(self.archive_data, f, indent=2)
    
    def should_archive(self, proposal: Dict) -> Tuple[bool, Optional[ArchiveReason]]:
        """
        Determine if proposal should be archived vs deleted
        Returns (should_archive, reason)
        """
        
        # Test proposals - don't archive
        author = proposal.get('author', '')
        if 'Test' in author or 'test' in author.lower():
            return False, ArchiveReason.TEST_CLEANUP
        
        # Duplicate proposals without reviews - don't archive
        if len(proposal.get('reviews', [])) == 0:
            # Check for duplicates (simple check based on same title/section)
            duplicates = [p for p in self.proposals 
                         if p['id'] != proposal['id'] 
                         and p.get('title') == proposal.get('title')
                         and p.get('section') == proposal.get('section')]
            if duplicates:
                return False, ArchiveReason.DUPLICATE
        
        # Has reviews - archive
        if len(proposal.get('reviews', [])) > 0:
            return True, ArchiveReason.REJECTED
        
        # Accepted or integrated - definitely archive
        if proposal.get('status') in ['accepted', 'integrated']:
            return True, ArchiveReason.INTEGRATED
        
        # Old proposals without activity - archive with expiry
        if proposal.get('date'):
            try:
                proposal_date = datetime.fromisoformat(proposal['date'])
                age = datetime.now() - proposal_date
                if age.days > 30:
                    return True, ArchiveReason.EXPIRED
            except:
                pass
        
        # Default to not archiving if uncertain
        return False, None
    
    def archive_proposal(self, proposal: Dict, reason: ArchiveReason, 
                        historical_value: bool = False):
        """Archive a single proposal"""
        archived = {
            **proposal,
            'archive_reason': reason.value,
            'archive_date': datetime.now().isoformat(),
            'historical_value': historical_value
        }
        
        # Set retention period
        if historical_value:
            archived['preserve_until'] = None  # Keep indefinitely
        else:
            retention_days = self.archive_data['archive_policy']['retention_days']
            preserve_until = datetime.now() + timedelta(days=retention_days)
            archived['preserve_until'] = preserve_until.isoformat()
        
        self.archive_data['archived_proposals'].append(archived)
        self.archive_data['archive_policy']['total_archived'] += 1
    
    def cleanup_test_proposals(self, dry_run: bool = True) -> Dict:
        """
        Clean up test proposals
        Returns summary of actions taken
        """
        results = {
            'deleted': [],
            'archived': [],
            'kept': [],
            'errors': []
        }
        
        proposals_to_keep = []
        
        for proposal in self.proposals:
            try:
                should_archive, reason = self.should_archive(proposal)
                
                if not should_archive and reason == ArchiveReason.TEST_CLEANUP:
                    # Delete test proposals
                    results['deleted'].append({
                        'id': proposal['id'],
                        'title': proposal.get('title', 'Unknown'),
                        'author': proposal.get('author', 'Unknown')
                    })
                elif not should_archive and reason == ArchiveReason.DUPLICATE:
                    # Delete duplicates without value
                    results['deleted'].append({
                        'id': proposal['id'],
                        'title': proposal.get('title', 'Unknown'),
                        'reason': 'duplicate'
                    })
                elif should_archive and reason:
                    # Archive valuable proposals
                    if not dry_run:
                        # Determine historical value
                        historical = proposal.get('status') in ['accepted', 'integrated']
                        self.archive_proposal(proposal, reason, historical)
                    results['archived'].append({
                        'id': proposal['id'],
                        'title': proposal.get('title', 'Unknown'),
                        'reason': reason.value
                    })
                else:
                    # Keep active proposals
                    proposals_to_keep.append(proposal)
                    results['kept'].append({
                        'id': proposal['id'],
                        'title': proposal.get('title', 'Unknown')
                    })
            except Exception as e:
                results['errors'].append({
                    'id': proposal.get('id', 'Unknown'),
                    'error': str(e)
                })
                # Keep proposal if error
                proposals_to_keep.append(proposal)
        
        if not dry_run:
            # Update proposals list
            self.proposals = proposals_to_keep
            self.save_proposals()
            self.save_archive()
        
        return results
    
    def generate_report(self, results: Dict) -> str:
        """Generate a human-readable report of cleanup results"""
        report = []
        report.append("=" * 60)
        report.append("PROPOSAL CLEANUP REPORT")
        report.append("=" * 60)
        report.append(f"Timestamp: {datetime.now().isoformat()}")
        report.append("")
        
        # Summary
        report.append("SUMMARY:")
        report.append(f"  Deleted: {len(results['deleted'])}")
        report.append(f"  Archived: {len(results['archived'])}")
        report.append(f"  Kept Active: {len(results['kept'])}")
        report.append(f"  Errors: {len(results['errors'])}")
        report.append("")
        
        # Details
        if results['deleted']:
            report.append("DELETED PROPOSALS:")
            for item in results['deleted']:
                reason = item.get('reason', 'test')
                report.append(f"  - [{item['id']}] {item['title'][:50]}... ({reason})")
                report.append(f"    Author: {item.get('author', 'Unknown')}")
            report.append("")
        
        if results['archived']:
            report.append("ARCHIVED PROPOSALS:")
            for item in results['archived']:
                report.append(f"  - [{item['id']}] {item['title'][:50]}...")
                report.append(f"    Reason: {item['reason']}")
            report.append("")
        
        if results['kept']:
            report.append("KEPT ACTIVE:")
            for item in results['kept']:
                report.append(f"  - [{item['id']}] {item['title'][:50]}...")
            report.append("")
        
        if results['errors']:
            report.append("ERRORS:")
            for item in results['errors']:
                report.append(f"  - [{item['id']}] Error: {item['error']}")
            report.append("")
        
        report.append("=" * 60)
        return "\n".join(report)

def main():
    """Main cleanup function"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Clean up test proposals")
    parser.add_argument('--dry-run', action='store_true', default=True,
                       help='Show what would be done without making changes')
    parser.add_argument('--execute', action='store_true',
                       help='Actually perform the cleanup')
    parser.add_argument('--config', type=str,
                       help='Path to config directory')
    
    args = parser.parse_args()
    
    # If --execute is specified, turn off dry-run
    if args.execute:
        args.dry_run = False
    
    print(f"Proposal Cleanup - {'DRY RUN' if args.dry_run else 'EXECUTING'}")
    print("-" * 60)
    
    # Initialize cleanup
    cleanup = ProposalCleanup(args.config)
    
    # Show current state
    print(f"Current proposals: {len(cleanup.proposals)}")
    print(f"Archived proposals: {len(cleanup.archive_data['archived_proposals'])}")
    print()
    
    # Run cleanup
    results = cleanup.cleanup_test_proposals(dry_run=args.dry_run)
    
    # Generate and print report
    report = cleanup.generate_report(results)
    print(report)
    
    if args.dry_run:
        print("\nThis was a DRY RUN. No changes were made.")
        print("Run with --execute to actually perform the cleanup.")

if __name__ == "__main__":
    main()