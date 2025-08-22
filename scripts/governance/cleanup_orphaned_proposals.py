#!/usr/bin/env python3
"""
Clean up orphaned proposal markdown files that don't have corresponding JSON entries
"""

import json
import os
from pathlib import Path
import re

def cleanup_orphaned_proposals():
    """Remove proposal markdown files that aren't in the JSON"""
    
    # Load current proposals from JSON
    config_path = Path(__file__).parent / "config"
    proposals_file = config_path / "whitepaper_proposals.json"
    
    active_proposals = set()
    if proposals_file.exists():
        with open(proposals_file, 'r') as f:
            data = json.load(f)
            for proposal in data.get('proposals', []):
                # Store as "section/id" 
                key = f"{proposal['section']}/{proposal['id']}"
                active_proposals.add(key)
    
    print(f"Active proposals in JSON: {len(active_proposals)}")
    for key in active_proposals:
        print(f"  - {key}")
    
    # Find all proposal markdown files
    whitepaper_path = Path(__file__).parent.parent.parent / "whitepaper" / "sections"
    proposal_files = list(whitepaper_path.glob("**/meta/proposals/*.md"))
    
    print(f"\nFound {len(proposal_files)} proposal markdown files")
    
    # Check each file
    orphaned = []
    kept = []
    
    for filepath in proposal_files:
        # Extract section and proposal ID from path
        # Path like: .../04-fundamental-concepts/01-universe-grid/meta/proposals/005-consciousness-integration.md
        parts = filepath.parts
        
        # Find the sections index
        try:
            sections_idx = parts.index("sections")
            # Section path is everything between sections and meta
            meta_idx = parts.index("meta")
            section_parts = parts[sections_idx+1:meta_idx]
            section = "/".join(section_parts)
            
            # Extract ID from filename (e.g., "005-consciousness-integration.md" -> "005")
            filename = filepath.name
            match = re.match(r'^(\d+)', filename)
            if match:
                prop_id = match.group(1).lstrip('0')  # Remove leading zeros
                if prop_id == '':
                    prop_id = '0'
                
                # Check if this proposal is active
                key = f"{section}/{prop_id}"
                
                if key not in active_proposals:
                    orphaned.append(filepath)
                    print(f"  ❌ Orphaned: {section}/{filepath.name}")
                else:
                    kept.append(filepath)
                    print(f"  ✓ Active: {section}/{filepath.name}")
        except (ValueError, IndexError):
            print(f"  ⚠️  Could not parse: {filepath}")
            continue
    
    print(f"\nSummary:")
    print(f"  Orphaned files to remove: {len(orphaned)}")
    print(f"  Active files to keep: {len(kept)}")
    
    if orphaned:
        response = input(f"\nRemove {len(orphaned)} orphaned proposal files? (y/n): ")
        if response.lower() == 'y':
            for filepath in orphaned:
                filepath.unlink()
                print(f"  Removed: {filepath.name}")
            print(f"\n✅ Removed {len(orphaned)} orphaned proposal files")
        else:
            print("Aborted - no files removed")
    else:
        print("\nNo orphaned files to remove!")
    
    return len(orphaned)

if __name__ == "__main__":
    cleanup_orphaned_proposals()