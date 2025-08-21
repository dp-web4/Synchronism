#!/usr/bin/env python3
"""
Proposal Sync System
Keeps JSON and markdown files in sync for proposal management
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from enum import Enum

class SyncAction(Enum):
    """Actions to take during sync"""
    CREATE_MARKDOWN = "create_markdown"
    CREATE_JSON = "create_json"
    UPDATE_MARKDOWN = "update_markdown"
    UPDATE_JSON = "update_json"
    DELETE_MARKDOWN = "delete_markdown"
    DELETE_JSON = "delete_json"
    IN_SYNC = "in_sync"

class ProposalSync:
    """Manages synchronization between JSON and markdown proposal files"""
    
    def __init__(self, base_path: str = None):
        if base_path is None:
            self.base_path = Path(__file__).parent.parent.parent
        else:
            self.base_path = Path(base_path)
            
        self.config_path = self.base_path / "scripts" / "governance" / "config"
        self.sections_path = self.base_path / "whitepaper" / "sections"
        self.proposals_json = self.config_path / "whitepaper_proposals.json"
        
    def load_json_proposals(self) -> Dict[str, Dict]:
        """Load proposals from JSON, keyed by section/id"""
        proposals = {}
        if self.proposals_json.exists():
            with open(self.proposals_json, 'r') as f:
                data = json.load(f)
                for prop in data.get('proposals', []):
                    key = f"{prop['section']}/{prop['id']}"
                    proposals[key] = prop
        return proposals
    
    def scan_markdown_proposals(self) -> Dict[str, Path]:
        """Scan for markdown proposal files, keyed by section/id"""
        markdown_files = {}
        for filepath in self.sections_path.glob("**/meta/proposals/*.md"):
            # Extract section and ID from path
            parts = filepath.parts
            try:
                sections_idx = parts.index("sections")
                meta_idx = parts.index("meta")
                section_parts = parts[sections_idx+1:meta_idx]
                section = "/".join(section_parts)
                
                # Extract ID from filename (e.g., "002-title.md" -> "002")
                filename = filepath.name
                prop_id = filename.split('-')[0]  # Keep leading zeros
                
                key = f"{section}/{prop_id}"
                markdown_files[key] = filepath
            except (ValueError, IndexError):
                print(f"Warning: Could not parse {filepath}")
                
        return markdown_files
    
    def create_markdown_from_json(self, proposal: Dict) -> Path:
        """Create a markdown file from JSON proposal data"""
        section = proposal['section']
        prop_id = proposal['id']
        
        # Create the directory structure
        section_path = self.sections_path / section
        proposals_dir = section_path / "meta" / "proposals"
        proposals_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate filename
        title_slug = proposal.get('title', 'untitled').lower()
        title_slug = title_slug.replace(' ', '-')[:50]  # Limit length
        # Remove special characters
        title_slug = ''.join(c if c.isalnum() or c == '-' else '-' for c in title_slug)
        filename = f"{prop_id:0>3}-{title_slug}.md"
        filepath = proposals_dir / filename
        
        # Create markdown content
        content = f"""# Proposal {prop_id}: {proposal.get('title', 'Untitled')}

## Metadata
- **ID**: {prop_id}
- **Author**: {proposal.get('author', 'Unknown')}
- **Date**: {proposal.get('date', 'Unknown')}
- **Status**: {proposal.get('status', 'submitted')}
- **Type**: {proposal.get('type', 'UNKNOWN')}
- **Section**: {section}

## Proposed Change
{proposal.get('content_summary', 'No content provided')}

## Rationale
{proposal.get('rationale', 'No rationale provided')}

## Specific Text Changes
{proposal.get('specific_changes', 'No specific changes provided')}

## Impact Assessment
- **Compatibility**: {proposal.get('compatibility', 'To be assessed')}
- **Dependencies**: {proposal.get('dependencies', 'To be identified')}
- **Risk**: {proposal.get('risk', 'To be evaluated')}

## Reviews
"""
        
        # Add reviews if any
        reviews = proposal.get('reviews', [])
        if reviews:
            for review in reviews:
                content += f"""
### Review by {review.get('reviewer', 'Unknown')}
- **Date**: {review.get('date', 'Unknown')}
- **Recommendation**: {review.get('recommendation', 'None')}
- **Comments**: {review.get('comments', 'No comments')}
"""
        else:
            content += "\n*No reviews yet*\n"
        
        # Add history if any
        if 'last_updated' in proposal:
            content += f"\n## History\n- Last Updated: {proposal['last_updated']}\n"
        
        # Write the file
        with open(filepath, 'w') as f:
            f.write(content)
            
        return filepath
    
    def analyze_sync_status(self) -> Dict[str, Tuple[SyncAction, Optional[Dict], Optional[Path]]]:
        """
        Analyze what needs to be synced
        Returns dict of key -> (action, json_data, markdown_path)
        """
        json_proposals = self.load_json_proposals()
        markdown_files = self.scan_markdown_proposals()
        
        all_keys = set(json_proposals.keys()) | set(markdown_files.keys())
        sync_status = {}
        
        for key in all_keys:
            has_json = key in json_proposals
            has_markdown = key in markdown_files
            
            if has_json and has_markdown:
                # Both exist - check if in sync (simplified - just mark as in sync)
                sync_status[key] = (
                    SyncAction.IN_SYNC,
                    json_proposals[key],
                    markdown_files[key]
                )
            elif has_json and not has_markdown:
                # JSON exists but no markdown - need to create markdown
                sync_status[key] = (
                    SyncAction.CREATE_MARKDOWN,
                    json_proposals[key],
                    None
                )
            elif not has_json and has_markdown:
                # Markdown exists but no JSON - orphaned markdown
                sync_status[key] = (
                    SyncAction.DELETE_MARKDOWN,
                    None,
                    markdown_files[key]
                )
            
        return sync_status
    
    def sync_proposals(self, dry_run: bool = True) -> Dict:
        """
        Synchronize JSON and markdown proposals
        """
        sync_status = self.analyze_sync_status()
        
        results = {
            'created_markdown': [],
            'deleted_markdown': [],
            'in_sync': [],
            'errors': []
        }
        
        print("=" * 60)
        print("PROPOSAL SYNC ANALYSIS")
        print("=" * 60)
        
        for key, (action, json_data, md_path) in sync_status.items():
            if action == SyncAction.CREATE_MARKDOWN:
                print(f"  📝 Need to create markdown for: {key}")
                if not dry_run:
                    try:
                        filepath = self.create_markdown_from_json(json_data)
                        results['created_markdown'].append(str(filepath))
                        print(f"     ✅ Created: {filepath.name}")
                    except Exception as e:
                        results['errors'].append(f"{key}: {str(e)}")
                        print(f"     ❌ Error: {str(e)}")
                else:
                    results['created_markdown'].append(f"Would create: {key}")
                    
            elif action == SyncAction.DELETE_MARKDOWN:
                print(f"  🗑️  Orphaned markdown to delete: {key}")
                print(f"     File: {md_path.name}")
                if not dry_run:
                    try:
                        md_path.unlink()
                        results['deleted_markdown'].append(str(md_path))
                        print(f"     ✅ Deleted")
                    except Exception as e:
                        results['errors'].append(f"{key}: {str(e)}")
                        print(f"     ❌ Error: {str(e)}")
                else:
                    results['deleted_markdown'].append(f"Would delete: {md_path.name}")
                    
            elif action == SyncAction.IN_SYNC:
                results['in_sync'].append(key)
                print(f"  ✅ In sync: {key}")
        
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"  In Sync: {len(results['in_sync'])}")
        print(f"  Need Markdown Creation: {len([k for k, (a, _, _) in sync_status.items() if a == SyncAction.CREATE_MARKDOWN])}")
        print(f"  Need Markdown Deletion: {len([k for k, (a, _, _) in sync_status.items() if a == SyncAction.DELETE_MARKDOWN])}")
        print(f"  Errors: {len(results['errors'])}")
        
        if dry_run:
            print("\nThis was a DRY RUN. Use --execute to perform sync.")
        else:
            print("\nSync completed!")
            
        return results
    
    def restore_from_json(self) -> Dict:
        """
        Restore all markdown files from JSON data
        This is useful after accidental deletion
        """
        json_proposals = self.load_json_proposals()
        results = {
            'restored': [],
            'errors': []
        }
        
        print("Restoring proposal markdown files from JSON...")
        
        for key, proposal in json_proposals.items():
            try:
                filepath = self.create_markdown_from_json(proposal)
                results['restored'].append({
                    'key': key,
                    'file': str(filepath),
                    'title': proposal.get('title', 'Untitled')
                })
                print(f"  ✅ Restored: {key} -> {filepath.name}")
            except Exception as e:
                results['errors'].append({
                    'key': key,
                    'error': str(e)
                })
                print(f"  ❌ Failed: {key} - {str(e)}")
        
        print(f"\nRestored {len(results['restored'])} proposals")
        if results['errors']:
            print(f"Failed to restore {len(results['errors'])} proposals")
            
        return results

def main():
    """Main sync function"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Sync proposal JSON and markdown files")
    parser.add_argument('--dry-run', action='store_true', default=True,
                       help='Show what would be done without making changes')
    parser.add_argument('--execute', action='store_true',
                       help='Actually perform the sync')
    parser.add_argument('--restore', action='store_true',
                       help='Restore all markdown files from JSON')
    
    args = parser.parse_args()
    
    if args.execute:
        args.dry_run = False
    
    sync = ProposalSync()
    
    if args.restore:
        # Special mode to restore from JSON
        results = sync.restore_from_json()
        print(f"\nRestore complete: {len(results['restored'])} files created")
    else:
        # Normal sync mode
        results = sync.sync_proposals(dry_run=args.dry_run)
    
    return 0

if __name__ == "__main__":
    exit(main())