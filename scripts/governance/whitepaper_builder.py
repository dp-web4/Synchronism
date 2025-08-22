#!/usr/bin/env python3
"""
Whitepaper Builder for Governance System
Automatically rebuilds MD/PDF/Web versions when fractal files change
"""

import subprocess
import hashlib
import shutil
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional

class WhitepaperBuilder:
    """Manages automatic rebuilding of whitepaper formats"""
    
    def __init__(self, base_path: str = "/mnt/c/exe/projects/ai-agents/Synchronism"):
        self.base_path = Path(base_path)
        self.whitepaper_path = self.base_path / "whitepaper"
        self.build_path = self.whitepaper_path / "build"
        self.sections_path = self.whitepaper_path / "sections"
        
        # Final destination for published versions
        self.docs_path = self.base_path / "docs" / "whitepaper"
        
        # Build scripts (using actual script names)
        self.build_scripts = {
            'markdown': self.whitepaper_path / "make-md.sh",
            'pdf': self.whitepaper_path / "make-pdf.sh",
            'web': self.whitepaper_path / "make-web-clean.sh"
        }
        
        # Build outputs (actual filenames in build directory)
        self.build_outputs = {
            'markdown': self.build_path / "Synchronism_Whitepaper_Complete.md",
            'pdf': self.build_path / "Synchronism_Whitepaper.pdf",
            'web': self.build_path / "web-clean"  # Web is a directory
        }
        
        # Final destinations in docs/whitepaper
        self.final_destinations = {
            'markdown': self.docs_path / "Synchronism_Whitepaper_Complete.md",
            'pdf': self.docs_path / "Synchronism_Whitepaper.pdf",
            'web': self.docs_path  # Web files go directly into docs/whitepaper
        }
        
        # Track build history
        self.build_log = self.base_path / ".governance" / "build_log.json"
    
    def detect_changes(self) -> bool:
        """Detect if any fractal files have changed since last build"""
        
        # Get all fractal files (non-meta .md files)
        fractal_files = []
        for md_file in self.sections_path.rglob("*.md"):
            if "meta" not in md_file.parts:
                fractal_files.append(md_file)
        
        if not fractal_files:
            return False
        
        # Compute combined hash of all fractal files
        hasher = hashlib.sha256()
        for file_path in sorted(fractal_files):
            if file_path.exists():
                hasher.update(file_path.read_bytes())
        
        current_hash = hasher.hexdigest()
        
        # Compare with last build hash
        last_hash = self._get_last_build_hash()
        
        if current_hash != last_hash:
            self._save_build_hash(current_hash)
            return True
        
        return False
    
    def _get_last_build_hash(self) -> str:
        """Get hash from last build"""
        import json
        
        if self.build_log.exists():
            with open(self.build_log, 'r') as f:
                data = json.load(f)
                return data.get('last_content_hash', '')
        return ''
    
    def _save_build_hash(self, content_hash: str):
        """Save current build hash"""
        import json
        
        self.build_log.parent.mkdir(parents=True, exist_ok=True)
        
        data = {
            'last_content_hash': content_hash,
            'last_build_check': datetime.now().isoformat()
        }
        
        with open(self.build_log, 'w') as f:
            json.dump(data, f, indent=2)
    
    def build_format(self, format_name: str) -> Tuple[bool, str]:
        """Build a specific format"""
        
        if format_name not in self.build_scripts:
            return False, f"Unknown format: {format_name}"
        
        script = self.build_scripts[format_name]
        if not script.exists():
            return False, f"Build script not found: {script}"
        
        print(f"Building {format_name} version...")
        
        try:
            # Make script executable
            script.chmod(0o755)
            
            # Run build script
            result = subprocess.run(
                [str(script)],
                cwd=str(self.whitepaper_path),
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode == 0:
                output_file = self.build_outputs.get(format_name)
                
                # Copy to docs/whitepaper after successful build
                if output_file and output_file.exists():
                    success_msg = self._copy_to_docs(format_name, output_file)
                    
                    if format_name == 'web':
                        # Web is a directory, check index.html
                        index_file = output_file / "index.html"
                        if index_file.exists():
                            size_kb = index_file.stat().st_size / 1024
                            return True, f"Built {format_name}: web directory ({size_kb:.1f} KB index.html) - {success_msg}"
                    else:
                        size_kb = output_file.stat().st_size / 1024
                        return True, f"Built {format_name}: {output_file.name} ({size_kb:.1f} KB) - {success_msg}"
                else:
                    return True, f"Built {format_name} (output location varies)"
            else:
                return False, f"Build failed: {result.stderr[:200]}"
                
        except subprocess.TimeoutExpired:
            return False, f"Build timeout for {format_name}"
        except Exception as e:
            return False, f"Build error: {str(e)}"
    
    def _copy_to_docs(self, format_name: str, source: Path) -> str:
        """Copy built files to docs/whitepaper directory"""
        
        # Ensure docs/whitepaper exists
        self.docs_path.mkdir(parents=True, exist_ok=True)
        
        destination = self.final_destinations.get(format_name)
        if not destination:
            return "No destination configured"
        
        try:
            if format_name == 'web':
                # Web files are already copied by make-web-clean.sh
                # Just verify they exist
                index_file = destination / "index.html"
                if index_file.exists():
                    return f"Web files in {destination.relative_to(self.base_path)}"
                else:
                    return "Web files may not have been copied correctly"
            else:
                # Single file copy
                shutil.copy2(source, destination)
                return f"Copied to {destination.relative_to(self.base_path)}"
        except Exception as e:
            return f"Copy failed: {str(e)}"
    
    def build_all_formats(self) -> Dict[str, Tuple[bool, str]]:
        """Build all whitepaper formats"""
        
        results = {}
        
        # Build in order: markdown first (base), then PDF, then web
        build_order = ['markdown', 'pdf', 'web']
        
        for format_name in build_order:
            success, message = self.build_format(format_name)
            results[format_name] = (success, message)
            
            if success:
                print(f"  ✅ {message}")
            else:
                print(f"  ❌ {format_name}: {message}")
        
        return results
    
    def rebuild_if_changed(self) -> Optional[Dict]:
        """Rebuild all formats if fractal files have changed"""
        
        if not self.detect_changes():
            print("No changes to fractal files detected - skipping rebuild")
            return None
        
        print("\n" + "="*60)
        print(" Fractal files changed - rebuilding whitepaper")
        print("="*60)
        
        # Record build start
        build_info = {
            'timestamp': datetime.now().isoformat(),
            'trigger': 'fractal_file_changes',
            'formats': {}
        }
        
        # Build all formats
        results = self.build_all_formats()
        
        # Record results
        for format_name, (success, message) in results.items():
            build_info['formats'][format_name] = {
                'success': success,
                'message': message
            }
        
        # Update build log
        self._update_build_log(build_info)
        
        # Summary
        successful = sum(1 for s, _ in results.values() if s)
        print(f"\nBuild complete: {successful}/{len(results)} formats built successfully")
        
        return build_info
    
    def _update_build_log(self, build_info: Dict):
        """Update build log with latest build info"""
        import json
        
        log_file = self.build_log
        
        # Load existing log
        if log_file.exists():
            with open(log_file, 'r') as f:
                data = json.load(f)
                # Ensure builds key exists
                if 'builds' not in data:
                    data['builds'] = []
        else:
            data = {'builds': []}
        
        # Add new build info
        data['builds'].append(build_info)
        data['last_build'] = build_info['timestamp']
        
        # Keep only last 10 builds
        if len(data['builds']) > 10:
            data['builds'] = data['builds'][-10:]
        
        # Save updated log
        with open(log_file, 'w') as f:
            json.dump(data, f, indent=2)
    
    def force_rebuild(self) -> Dict:
        """Force rebuild regardless of changes"""
        
        print("\n" + "="*60)
        print(" Force rebuilding whitepaper (manual trigger)")
        print("="*60)
        
        build_info = {
            'timestamp': datetime.now().isoformat(),
            'trigger': 'force_rebuild',
            'formats': {}
        }
        
        results = self.build_all_formats()
        
        for format_name, (success, message) in results.items():
            build_info['formats'][format_name] = {
                'success': success,
                'message': message
            }
        
        self._update_build_log(build_info)
        
        return build_info

def integrate_with_governance_cycle():
    """Integration point for governance cycle completion"""
    
    from changelog_manager import ChangelogManager
    
    # This would be called at the end of each governance cycle
    def complete_cycle_with_rebuild(cycle_results: Dict) -> Dict:
        """Complete governance cycle and rebuild if needed"""
        
        # Check if any proposals were implemented
        implemented = cycle_results.get('implemented_proposals', [])
        
        if not implemented:
            print("No proposals implemented - skipping rebuild")
            return cycle_results
        
        # Check for actual fractal file changes
        changelog_mgr = ChangelogManager()
        builder = WhitepaperBuilder()
        
        # Detect changes and rebuild if needed
        build_result = builder.rebuild_if_changed()
        
        if build_result:
            cycle_results['whitepaper_rebuilt'] = True
            cycle_results['build_info'] = build_result
            
            # Log in main changelog
            if build_result['formats'].get('pdf', {}).get('success'):
                print("\n✅ Whitepaper PDF regenerated with latest changes")
            if build_result['formats'].get('web', {}).get('success'):
                print("✅ Whitepaper web version updated")
        else:
            cycle_results['whitepaper_rebuilt'] = False
        
        return cycle_results
    
    return complete_cycle_with_rebuild

def main():
    """Test the whitepaper builder"""
    
    builder = WhitepaperBuilder()
    
    print("Whitepaper Builder Test")
    print("-" * 40)
    
    # Check for changes
    has_changes = builder.detect_changes()
    print(f"Fractal files changed: {has_changes}")
    
    if has_changes:
        # Rebuild all formats
        result = builder.rebuild_if_changed()
        if result:
            print(f"\nBuild completed at: {result['timestamp']}")
    else:
        print("\nNo changes detected. Use --force to rebuild anyway")
        
        # For testing, show force rebuild
        import sys
        if '--force' in sys.argv:
            result = builder.force_rebuild()
            print(f"\nForced rebuild completed at: {result['timestamp']}")

if __name__ == "__main__":
    main()