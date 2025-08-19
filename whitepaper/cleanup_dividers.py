#!/usr/bin/env python3

import os
import glob

def clean_file(filepath):
    """Remove trailing --- and empty lines from a markdown file."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    if not lines:
        return False
    
    # Remove trailing empty lines and ---
    while lines and (lines[-1].strip() == '' or lines[-1].strip() == '---'):
        lines.pop()
    
    # Write back
    with open(filepath, 'w') as f:
        f.writelines(lines)
    
    return True

# Process all markdown files
count = 0
for filepath in glob.glob('sections/**/*.md', recursive=True):
    if clean_file(filepath):
        print(f"Cleaned: {os.path.basename(filepath)}")
        count += 1

print(f"\nTotal files cleaned: {count}")