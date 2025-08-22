#!/usr/bin/env python3
"""
Script to properly add a new subsection to the whitepaper
Handles numbering, navigation, file organization, and changelog
"""

import os
import sys
import json
import shutil
from pathlib import Path
from datetime import datetime
import re

class SubsectionManager:
    def __init__(self, base_path: str = None):
        """Initialize the subsection manager"""
        self.base_path = Path(base_path or os.getcwd())
        self.sections_path = self.base_path / "whitepaper" / "sections"
        
    def get_next_subsection_number(self, section_dir: str) -> tuple[int, str]:
        """Get the next available subsection number"""
        section_path = self.sections_path / section_dir
        
        # Find existing numbered subdirectories
        existing_nums = []
        for item in section_path.iterdir():
            if item.is_dir():
                match = re.match(r'^(\d+)-', item.name)
                if match:
                    existing_nums.append(int(match.group(1)))
        
        # Get next number
        next_num = max(existing_nums) + 1 if existing_nums else 1
        
        # Format with section number (e.g., 4.13)
        section_num = section_dir.split('-')[0]
        # Remove leading zero from section number if present
        section_num = str(int(section_num))
        formatted = f"{section_num}.{next_num}"
        
        return next_num, formatted
    
    def add_subsection(self, 
                      section_dir: str,
                      subsection_name: str,
                      title: str,
                      content: str,
                      author: str = "System",
                      rationale: str = ""):
        """Add a new subsection properly"""
        
        # Get numbering
        next_num, formatted_num = self.get_next_subsection_number(section_dir)
        
        # Create directory name (e.g., "13-compression-trust")
        dir_name = f"{next_num:02d}-{subsection_name}"
        subsection_path = self.sections_path / section_dir / dir_name
        
        # Create the directory structure
        subsection_path.mkdir(parents=True, exist_ok=True)
        (subsection_path / "meta").mkdir(exist_ok=True)
        (subsection_path / "archive").mkdir(exist_ok=True)
        
        # Write the main content file with proper heading
        content_file = subsection_path / f"{subsection_name.replace('-', '_')}.md"
        
        # Ensure content starts with proper section number
        if not content.startswith(f"## {formatted_num}"):
            # Extract the title from content if it has one
            lines = content.split('\n')
            if lines[0].startswith('#'):
                # Replace the heading with numbered version
                lines[0] = f"## {formatted_num} {title}"
                content = '\n'.join(lines)
            else:
                # Add numbered heading
                content = f"## {formatted_num} {title}\n\n{content}"
        
        with open(content_file, 'w') as f:
            f.write(content)
        
        # Update the parent section's changelog
        self.update_changelog(section_dir, formatted_num, title, author, rationale)
        
        # Update navigation if needed
        self.update_navigation(section_dir, next_num, subsection_name, title)
        
        print(f"✅ Added subsection {formatted_num}: {title}")
        print(f"   Location: {subsection_path}")
        print(f"   Content file: {content_file.name}")
        
        return subsection_path
    
    def update_changelog(self, section_dir: str, num: str, title: str, 
                        author: str, rationale: str):
        """Update the section's changelog"""
        changelog_path = self.sections_path / section_dir / "meta" / "CHANGELOG.md"
        
        # Read existing or create new
        if changelog_path.exists():
            with open(changelog_path) as f:
                content = f.read()
        else:
            content = "# Section Changelog\n\n## Format\n"
            content += "- **Date** | **Author** | **Type**\n"
            content += "  - Description: Brief description of change\n"
            content += "  - Rationale: Explanation for why the change was made\n\n"
            content += "## Entries\n\n"
        
        # Add new entry at the top of entries
        date = datetime.now().strftime("%Y-%m-%d")
        new_entry = f"### {date} | {author} | ADD\n"
        new_entry += f"- **Description**: Added subsection {num}: {title}\n"
        new_entry += f"- **Rationale**: {rationale}\n\n"
        
        # Insert after "## Entries" line
        lines = content.split('\n')
        for i, line in enumerate(lines):
            if line.strip() == "## Entries":
                lines.insert(i + 2, new_entry)
                break
        
        content = '\n'.join(lines)
        
        # Write back
        changelog_path.parent.mkdir(exist_ok=True)
        with open(changelog_path, 'w') as f:
            f.write(content)
    
    def update_navigation(self, section_dir: str, num: int, 
                         subsection_name: str, title: str):
        """Update navigation structure if needed"""
        # This would update any navigation files or indexes
        # For now, we'll update a section index if it exists
        
        index_path = self.sections_path / section_dir / "index.md"
        if index_path.exists():
            with open(index_path) as f:
                content = f.read()
            
            # Add link to new subsection
            new_link = f"{num:02d}. [{title}](./{num:02d}-{subsection_name}/)\n"
            
            # Find where to insert (at the end of the list)
            lines = content.split('\n')
            inserted = False
            for i in range(len(lines) - 1, -1, -1):
                if re.match(r'^\d+\.\s+\[', lines[i]):
                    lines.insert(i + 1, new_link)
                    inserted = True
                    break
            
            if not inserted:
                # No existing links, add at end
                lines.append(new_link)
            
            with open(index_path, 'w') as f:
                f.write('\n'.join(lines))

def main():
    """Main function for CLI usage"""
    if len(sys.argv) < 5:
        print("Usage: add_subsection.py <section_dir> <subsection_name> <title> <content_file> [author] [rationale]")
        print("Example: add_subsection.py 04-fundamental-concepts compression-trust 'Compression, Trust, and Communication' content.md 'Claude' 'Reorganization per proposal'")
        sys.exit(1)
    
    section_dir = sys.argv[1]
    subsection_name = sys.argv[2]
    title = sys.argv[3]
    content_file = sys.argv[4]
    author = sys.argv[5] if len(sys.argv) > 5 else "System"
    rationale = sys.argv[6] if len(sys.argv) > 6 else "Added new subsection"
    
    # Read content
    with open(content_file) as f:
        content = f.read()
    
    # Get base path (assuming we're in scripts/governance)
    base_path = Path(__file__).parent.parent.parent
    
    manager = SubsectionManager(base_path)
    manager.add_subsection(
        section_dir=section_dir,
        subsection_name=subsection_name,
        title=title,
        content=content,
        author=author,
        rationale=rationale
    )

if __name__ == "__main__":
    main()