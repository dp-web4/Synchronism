#!/usr/bin/env python3
"""
Changelog Manager for Whitepaper Governance
Only logs noteworthy changes when fractal files are actually modified
"""

import json
import hashlib
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from enum import Enum

class ChangeType(Enum):
    """Types of noteworthy changes"""
    CONTENT_ADDITION = "content_addition"
    CONTENT_REVISION = "content_revision"
    CONTENT_DELETION = "content_deletion"
    STRUCTURAL_CHANGE = "structural_change"
    PHILOSOPHICAL_REFINEMENT = "philosophical_refinement"
    MATHEMATICAL_FORMALIZATION = "mathematical_formalization"
    CLARIFICATION = "clarification"
    CORRECTION = "correction"

class ChangelogManager:
    """Manages changelog updates for actual content changes"""
    
    def __init__(self, base_path: str = "/mnt/c/exe/projects/ai-agents/Synchronism"):
        self.base_path = Path(base_path)
        self.whitepaper_path = self.base_path / "whitepaper"
        self.sections_path = self.whitepaper_path / "sections"
        self.global_changelog = self.whitepaper_path / "CHANGELOG.md"
        self.main_readme = self.base_path / "README.md"
        
        # Track file hashes to detect actual changes
        self.file_hashes = {}
        self._load_file_hashes()
    
    def _load_file_hashes(self):
        """Load stored file hashes to detect changes"""
        hash_file = self.base_path / ".governance" / "file_hashes.json"
        if hash_file.exists():
            with open(hash_file, 'r') as f:
                self.file_hashes = json.load(f)
    
    def _save_file_hashes(self):
        """Save file hashes for change detection"""
        hash_file = self.base_path / ".governance" / "file_hashes.json"
        hash_file.parent.mkdir(parents=True, exist_ok=True)
        with open(hash_file, 'w') as f:
            json.dump(self.file_hashes, f, indent=2)
    
    def _compute_file_hash(self, file_path: Path) -> str:
        """Compute hash of file content"""
        if not file_path.exists():
            return ""
        content = file_path.read_bytes()
        return hashlib.sha256(content).hexdigest()
    
    def detect_fractal_changes(self, section_path: str) -> List[Tuple[Path, bool]]:
        """Detect which fractal files have actually changed"""
        section_dir = self.sections_path / section_path
        changed_files = []
        
        # Check all .md files in section (excluding meta/)
        for md_file in section_dir.rglob("*.md"):
            # Skip meta files - they're not fractal content
            if "meta" in md_file.parts:
                continue
            
            # Compute current hash
            current_hash = self._compute_file_hash(md_file)
            
            # Compare with stored hash
            file_key = str(md_file.relative_to(self.base_path))
            old_hash = self.file_hashes.get(file_key, "")
            
            if current_hash != old_hash:
                changed_files.append((md_file, old_hash == ""))  # (file, is_new)
                self.file_hashes[file_key] = current_hash
        
        return changed_files
    
    def create_change_entry(self,
                           proposal_id: str,
                           proposal_title: str,
                           arbiter: str,
                           change_type: ChangeType,
                           files_changed: List[Path],
                           description: str) -> Dict:
        """Create a noteworthy change entry"""
        return {
            "timestamp": datetime.now().isoformat(),
            "proposal_id": proposal_id,
            "title": proposal_title,
            "arbiter": arbiter,
            "type": change_type.value,
            "files_changed": [str(f.relative_to(self.sections_path)) for f in files_changed],
            "description": description
        }
    
    def update_local_changelog(self,
                              section_path: str,
                              change_entry: Dict) -> bool:
        """Update the local section changelog"""
        section_dir = self.sections_path / section_path
        local_changelog = section_dir / "meta" / "changelog.md"
        
        if not local_changelog.parent.exists():
            local_changelog.parent.mkdir(parents=True, exist_ok=True)
        
        # Read existing changelog
        if local_changelog.exists():
            content = local_changelog.read_text()
        else:
            content = f"# Changelog for {section_path}\n\n"
        
        # Add new entry (noteworthy changes only)
        entry_text = f"""
### {change_entry['timestamp'][:10]} [{change_entry['arbiter']}] [{change_entry['type'].upper()}]
- **Proposal**: {change_entry['proposal_id']} - {change_entry['title']}
- **Changes**: {', '.join(change_entry['files_changed'])}
- **Description**: {change_entry['description']}
"""
        
        content += entry_text
        local_changelog.write_text(content)
        return True
    
    def aggregate_to_global_changelog(self, change_entry: Dict) -> bool:
        """Add noteworthy change to global changelog"""
        
        # Read existing global changelog
        if self.global_changelog.exists():
            content = self.global_changelog.read_text()
        else:
            content = """# Synchronism Whitepaper Changelog

This changelog tracks noteworthy changes to the whitepaper content.
Only actual modifications to fractal files are logged here.

## Format
[Date] [Type] Section - Description (Proposal ID)

---

"""
        
        # Add new entry
        date = change_entry['timestamp'][:10]
        change_type = change_entry['type'].replace('_', ' ').title()
        section = change_entry['files_changed'][0].split('/')[0] if change_entry['files_changed'] else 'unknown'
        
        entry = f"- [{date}] [{change_type}] {section} - {change_entry['description']} (#{change_entry['proposal_id']})\n"
        
        # Insert after the format section
        lines = content.split('\n')
        insert_pos = -1
        for i, line in enumerate(lines):
            if line.startswith('---'):
                insert_pos = i + 2
                break
        
        if insert_pos > 0:
            lines.insert(insert_pos, entry)
            content = '\n'.join(lines)
        else:
            content += entry
        
        self.global_changelog.write_text(content)
        return True
    
    def update_main_readme_timestamp(self, last_change: Dict) -> bool:
        """Update main README with last meaningful change timestamp"""
        
        if not self.main_readme.exists():
            return False
        
        content = self.main_readme.read_text()
        
        # Look for the last change line
        timestamp_line = f"**Last meaningful change**: {last_change['timestamp'][:19]} - {last_change['title']}"
        
        # Replace or add the timestamp line
        lines = content.split('\n')
        found = False
        for i, line in enumerate(lines):
            if line.startswith("**Last meaningful change**:"):
                lines[i] = timestamp_line
                found = True
                break
        
        if not found:
            # Add after the main title
            for i, line in enumerate(lines):
                if line.startswith("# Synchronism"):
                    lines.insert(i + 2, timestamp_line)
                    lines.insert(i + 3, "")
                    break
        
        content = '\n'.join(lines)
        self.main_readme.write_text(content)
        return True
    
    def process_implemented_proposal(self,
                                    proposal_id: str,
                                    proposal_data: Dict,
                                    arbiter: str,
                                    section_path: str) -> Optional[Dict]:
        """Process an implemented proposal and update changelogs if noteworthy"""
        
        # Detect actual file changes
        changed_files = self.detect_fractal_changes(section_path)
        
        if not changed_files:
            # No actual changes to fractal files - not noteworthy
            return None
        
        # Determine change type based on proposal
        change_type = self._determine_change_type(proposal_data)
        
        # Create change entry
        change_entry = self.create_change_entry(
            proposal_id=proposal_id,
            proposal_title=proposal_data.get('title', 'Unknown'),
            arbiter=arbiter,
            change_type=change_type,
            files_changed=[f[0] for f in changed_files],
            description=proposal_data.get('rationale', 'No description')[:200]
        )
        
        # Update local changelog
        self.update_local_changelog(section_path, change_entry)
        
        # Aggregate to global changelog
        self.aggregate_to_global_changelog(change_entry)
        
        # Update main README timestamp
        self.update_main_readme_timestamp(change_entry)
        
        # Save updated file hashes
        self._save_file_hashes()
        
        return change_entry
    
    def _determine_change_type(self, proposal_data: Dict) -> ChangeType:
        """Determine the type of change from proposal data"""
        
        title = proposal_data.get('title', '').lower()
        content = proposal_data.get('content', '').lower()
        proposal_type = proposal_data.get('type', '').lower()
        
        # Check for specific change types
        if 'mathematical' in title or 'equation' in content or 'tensor' in content:
            return ChangeType.MATHEMATICAL_FORMALIZATION
        elif 'philosophical' in title or 'consciousness' in content:
            return ChangeType.PHILOSOPHICAL_REFINEMENT
        elif 'clarif' in proposal_type:
            return ChangeType.CLARIFICATION
        elif 'correct' in proposal_type:
            return ChangeType.CORRECTION
        elif 'structural' in title:
            return ChangeType.STRUCTURAL_CHANGE
        elif 'delet' in title or 'remov' in title:
            return ChangeType.CONTENT_DELETION
        elif 'revis' in proposal_type:
            return ChangeType.CONTENT_REVISION
        else:
            return ChangeType.CONTENT_ADDITION
    
    def get_recent_noteworthy_changes(self, limit: int = 10) -> List[Dict]:
        """Get recent noteworthy changes from global changelog"""
        
        if not self.global_changelog.exists():
            return []
        
        content = self.global_changelog.read_text()
        changes = []
        
        # Parse changelog entries
        for line in content.split('\n'):
            if line.startswith('- ['):
                # Parse the entry
                try:
                    parts = line.split('] ')
                    if len(parts) >= 3:
                        date = parts[0].replace('- [', '')
                        change_type = parts[1].replace('[', '')
                        rest = '] '.join(parts[2:])
                        
                        changes.append({
                            'date': date,
                            'type': change_type,
                            'description': rest
                        })
                except:
                    continue
        
        return changes[:limit]
    
    def cleanup_routine_logs(self):
        """Remove routine log entries that aren't noteworthy"""
        
        # Clean up governance logs that don't represent actual changes
        governance_logs = self.base_path / "scripts" / "governance" / "config"
        
        # Keep only recent reports (last 7 days)
        cutoff = datetime.now().timestamp() - (7 * 24 * 60 * 60)
        
        for log_file in governance_logs.glob("*_report_*.json"):
            if log_file.stat().st_mtime < cutoff:
                log_file.unlink()
                print(f"Removed old log: {log_file.name}")
        
        # Clean up empty cycle logs
        cycles_file = governance_logs / "cycles.json"
        if cycles_file.exists():
            with open(cycles_file, 'r') as f:
                data = json.load(f)
            
            # Keep only cycles with actual decisions
            meaningful_cycles = []
            for cycle in data.get('history', []):
                if cycle.get('summary', {}).get('decisions_made', 0) > 0:
                    meaningful_cycles.append(cycle)
            
            data['history'] = meaningful_cycles
            
            with open(cycles_file, 'w') as f:
                json.dump(data, f, indent=2)
            
            print(f"Cleaned up cycles: kept {len(meaningful_cycles)} meaningful cycles")

def test_changelog_manager():
    """Test the changelog management system"""
    
    print("\n" + "="*60)
    print(" Testing Changelog Manager")
    print("="*60)
    
    manager = ChangelogManager()
    
    # Simulate a proposal implementation
    test_proposal = {
        'id': 'test-001',
        'title': 'Add Consciousness Field Integration',
        'type': 'expansion',
        'content': 'Enhance consciousness emergence patterns',
        'rationale': 'Strengthens philosophical coherence with distributed intelligence framework'
    }
    
    # Simulate file changes
    test_section = "04-fundamental-concepts/01-universe-grid"
    test_file = manager.sections_path / test_section / "universe_grid.md"
    
    # Make a small change to detect
    if test_file.exists():
        original_content = test_file.read_text()
        test_file.write_text(original_content + "\n<!-- Test change -->")
        
        # Process the change
        change_entry = manager.process_implemented_proposal(
            proposal_id=test_proposal['id'],
            proposal_data=test_proposal,
            arbiter="Dennis",
            section_path=test_section
        )
        
        if change_entry:
            print(f"✅ Noteworthy change logged:")
            print(f"   Type: {change_entry['type']}")
            print(f"   Files: {change_entry['files_changed']}")
            print(f"   Description: {change_entry['description'][:50]}...")
        else:
            print("❌ No noteworthy changes detected")
        
        # Restore original
        test_file.write_text(original_content)
    
    # Test cleanup
    print("\n[Testing Log Cleanup]")
    manager.cleanup_routine_logs()
    
    # Show recent noteworthy changes
    print("\n[Recent Noteworthy Changes]")
    recent = manager.get_recent_noteworthy_changes(5)
    if recent:
        for change in recent:
            print(f"  - {change['date']}: {change['type']} - {change['description'][:50]}...")
    else:
        print("  No recent noteworthy changes")

if __name__ == "__main__":
    test_changelog_manager()