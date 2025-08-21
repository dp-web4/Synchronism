#!/usr/bin/env python3
"""
Governance Maintenance Script
Performs regular maintenance on the proposal system:
- Expires old proposals
- Cleans up duplicates
- Archives completed proposals
- Generates maintenance reports
"""

import json
import argparse
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# Import our modules
from proposal_cleanup import ProposalCleanup, ArchiveReason
from whitepaper_governance_enhanced import (
    EnhancedWhitepaperGovernance,
    WithdrawalReason
)

class GovernanceMaintenance:
    """Handles all governance maintenance tasks"""
    
    def __init__(self):
        self.governance = EnhancedWhitepaperGovernance()
        self.cleanup = ProposalCleanup()
        self.report = []
        
    def run_maintenance(self, dry_run: bool = True) -> Dict:
        """
        Run all maintenance tasks
        
        Returns:
            Summary of all actions taken
        """
        results = {
            'timestamp': datetime.now().isoformat(),
            'dry_run': dry_run,
            'expired': [],
            'duplicates_found': [],
            'archived': [],
            'errors': []
        }
        
        self.report.append("=" * 60)
        self.report.append("GOVERNANCE MAINTENANCE REPORT")
        self.report.append("=" * 60)
        self.report.append(f"Timestamp: {results['timestamp']}")
        self.report.append(f"Mode: {'DRY RUN' if dry_run else 'EXECUTING'}")
        self.report.append("")
        
        # 1. Check for expired proposals
        self.report.append("1. CHECKING FOR EXPIRED PROPOSALS")
        self.report.append("-" * 40)
        expired = self.governance.expire_old_proposals(dry_run=dry_run)
        if expired:
            for p in expired:
                self.report.append(f"  - [{p['id']}] {p.get('title', 'Unknown')[:40]}...")
                results['expired'].append(p['id'])
        else:
            self.report.append("  No proposals to expire")
        self.report.append("")
        
        # 2. Check for duplicates in active proposals
        self.report.append("2. CHECKING FOR DUPLICATE PROPOSALS")
        self.report.append("-" * 40)
        active = self.governance.get_active_proposals()
        checked_pairs = set()
        
        for proposal in active:
            # Check each proposal against others
            for other in active:
                if proposal['id'] != other['id']:
                    pair = tuple(sorted([proposal['id'], other['id']]))
                    if pair not in checked_pairs:
                        checked_pairs.add(pair)
                        
                        # Check similarity
                        dups = self.governance.check_duplicates(
                            proposal.get('section', ''),
                            other.get('title', ''),
                            other.get('content_summary', '')
                        )
                        
                        for dup_id, similarity in dups:
                            if dup_id == proposal['id'] and similarity > 0.85:
                                self.report.append(f"  - [{proposal['id']}] and [{other['id']}] are {similarity:.0%} similar")
                                results['duplicates_found'].append({
                                    'ids': [proposal['id'], other['id']],
                                    'similarity': similarity
                                })
        
        if not results['duplicates_found']:
            self.report.append("  No duplicates found")
        self.report.append("")
        
        # 3. Archive old completed proposals
        self.report.append("3. ARCHIVING COMPLETED PROPOSALS")
        self.report.append("-" * 40)
        proposals = self.governance.load_proposals()
        
        for proposal in proposals:
            if proposal.get('status') in ['rejected', 'implemented']:
                # Check age
                try:
                    date = datetime.fromisoformat(proposal.get('last_updated', proposal.get('date', '')))
                    age = (datetime.now() - date).days
                    
                    if age > 7:  # Archive after 7 days
                        if not dry_run:
                            # Archive it
                            reason = (ArchiveReason.INTEGRATED 
                                    if proposal['status'] == 'implemented' 
                                    else ArchiveReason.REJECTED)
                            self.cleanup.archive_proposal(proposal, reason, 
                                                         historical_value=True)
                            # Remove from active
                            proposals.remove(proposal)
                        
                        self.report.append(f"  - [{proposal['id']}] {proposal.get('title', '')[:40]}... ({proposal['status']})")
                        results['archived'].append(proposal['id'])
                except:
                    pass
        
        if not results['archived']:
            self.report.append("  No proposals to archive")
        self.report.append("")
        
        # Save changes if not dry run
        if not dry_run:
            if results['archived']:
                self.governance.save_proposals(proposals)
                self.cleanup.save_archive()
        
        # 4. Summary
        self.report.append("SUMMARY")
        self.report.append("-" * 40)
        self.report.append(f"  Expired: {len(results['expired'])}")
        self.report.append(f"  Duplicates Found: {len(results['duplicates_found'])}")
        self.report.append(f"  Archived: {len(results['archived'])}")
        self.report.append(f"  Errors: {len(results['errors'])}")
        
        if dry_run:
            self.report.append("")
            self.report.append("This was a DRY RUN. No changes were made.")
            self.report.append("Run with --execute to perform maintenance.")
        
        self.report.append("=" * 60)
        
        return results
    
    def get_report(self) -> str:
        """Get the formatted report"""
        return "\n".join(self.report)
    
    def save_report(self, filepath: Optional[str] = None):
        """Save report to file"""
        if filepath is None:
            # Default to maintenance_reports directory
            reports_dir = Path(__file__).parent / "maintenance_reports"
            reports_dir.mkdir(exist_ok=True)
            
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filepath = reports_dir / f"maintenance_{timestamp}.txt"
        
        with open(filepath, 'w') as f:
            f.write(self.get_report())
        
        return filepath

def main():
    """Main maintenance function"""
    parser = argparse.ArgumentParser(description="Governance System Maintenance")
    parser.add_argument('--dry-run', action='store_true', default=True,
                       help='Show what would be done without making changes')
    parser.add_argument('--execute', action='store_true',
                       help='Actually perform maintenance')
    parser.add_argument('--save-report', action='store_true',
                       help='Save report to file')
    
    args = parser.parse_args()
    
    # If --execute is specified, turn off dry-run
    if args.execute:
        args.dry_run = False
    
    # Run maintenance
    maintenance = GovernanceMaintenance()
    results = maintenance.run_maintenance(dry_run=args.dry_run)
    
    # Print report
    print(maintenance.get_report())
    
    # Save report if requested
    if args.save_report:
        filepath = maintenance.save_report()
        print(f"\nReport saved to: {filepath}")

if __name__ == "__main__":
    main()