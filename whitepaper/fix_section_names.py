#!/usr/bin/env python3
"""
Fix truncated section names by reading actual titles from files
"""

import os
import re

def get_title_from_file(filepath):
    """Extract the first heading from a markdown file"""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('#'):
                    # Remove the # symbols and clean up
                    title = re.sub(r'^#+\s*', '', line).strip()
                    # Remove leading numbers
                    title = re.sub(r'^\d+\s+', '', title)
                    title = re.sub(r'^\d+\.\d+\s+', '', title)
                    return title
    except:
        pass
    return None

def rename_truncated_files(sections_dir):
    """Rename files that were truncated"""
    renamed = []
    
    for root, dirs, files in os.walk(sections_dir):
        for filename in files:
            if filename.endswith('.md') and filename != 'index.md':
                filepath = os.path.join(root, filename)
                title = get_title_from_file(filepath)
                
                if title:
                    # Create new filename from title
                    # Keep the number prefix if it exists
                    prefix_match = re.match(r'^(\d+-)', filename)
                    prefix = prefix_match.group(1) if prefix_match else ''
                    
                    # Clean title for filename
                    clean_title = re.sub(r'[^\w\s-]', '', title.lower())
                    clean_title = re.sub(r'[-\s]+', '-', clean_title)
                    
                    # Limit length but try to keep complete words
                    if len(clean_title) > 50:
                        clean_title = clean_title[:50]
                        # Try to cut at word boundary
                        last_dash = clean_title.rfind('-')
                        if last_dash > 30:
                            clean_title = clean_title[:last_dash]
                    
                    new_filename = f"{prefix}{clean_title}.md"
                    
                    if new_filename != filename:
                        new_filepath = os.path.join(root, new_filename)
                        print(f"Renaming: {filename}")
                        print(f"      to: {new_filename}")
                        print(f"   Title: {title}")
                        print()
                        os.rename(filepath, new_filepath)
                        renamed.append((filename, new_filename, title))
    
    return renamed

if __name__ == "__main__":
    sections_dir = "sections"
    
    print("Fixing truncated section names...")
    print()
    
    renamed = rename_truncated_files(sections_dir)
    
    print(f"\nRenamed {len(renamed)} files")
    
    # Update the reorganize_sections.py file references if needed
    if renamed:
        print("\nYou may need to update reorganize_sections.py with the new filenames")