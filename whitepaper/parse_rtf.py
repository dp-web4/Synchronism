#!/usr/bin/env python3
"""
Parse Synchronism RTF-converted markdown into proper sections
"""

import re
import os

def clean_markdown(text):
    """Clean up RTF conversion artifacts"""
    # Remove RTF artifacts
    text = re.sub(r'\[.*?\]\{#.*?\}', '', text)  # Remove anchor tags
    text = re.sub(r'\{\.underline\}', '', text)  # Remove underline markup
    text = re.sub(r'þ', ':', text)  # Fix special characters
    text = re.sub(r'#{6,}', '###', text)  # Fix excessive heading levels
    text = re.sub(r'\n{3,}', '\n\n', text)  # Remove excessive blank lines
    
    return text.strip()

def parse_sections(filepath):
    """Parse markdown into logical sections"""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Clean the content
    content = clean_markdown(content)
    
    sections = []
    
    # Find major section headings (## or ###)
    pattern = r'^(#{1,3})\s+(.+?)$'
    
    lines = content.split('\n')
    current_section = None
    current_content = []
    
    for line in lines:
        match = re.match(pattern, line)
        if match:
            # Save previous section
            if current_section and current_content:
                sections.append({
                    'title': current_section,
                    'content': '\n'.join(current_content).strip()
                })
            
            # Start new section
            level = len(match.group(1))
            title = match.group(2).strip()
            
            # Skip table of contents entries
            if '...' in title or title.startswith('['):
                continue
                
            current_section = title
            current_content = []
        elif current_section:
            current_content.append(line)
    
    # Add final section
    if current_section and current_content:
        sections.append({
            'title': current_section,
            'content': '\n'.join(current_content).strip()
        })
    
    return sections

def create_section_files(sections):
    """Create clean markdown files for each section"""
    
    # Section name mapping
    section_map = {
        'Introduction': '02-introduction.md',
        'Importance of Perspective': '03-importance-of-perspective.md',
        'Relation to Hermetic Principles': '04-hermetic-principles.md',
        'Fundamental Concepts': '05-fundamental-concepts.md',
        'Alternative Perspective on Quantum': '06-quantum-perspective.md',
        'Implications and Applications': '07-implications.md',
        'Conclusion': '08-conclusion.md',
        'Appendix A': '09-appendix-mathematical.md',
    }
    
    # Keep executive summary and intro as is
    created_files = ['00-executive-summary.md', '01-introduction.md']
    
    for section in sections:
        title = section['title']
        
        # Find matching filename
        filename = None
        for key, fname in section_map.items():
            if key.lower() in title.lower():
                filename = fname
                break
        
        if not filename:
            # Generate filename from title
            clean_title = re.sub(r'[^\w\s-]', '', title.lower())
            clean_title = re.sub(r'[-\s]+', '-', clean_title)
            filename = f"{len(created_files):02d}-{clean_title[:30]}.md"
        
        # Skip if already created
        if filename in created_files:
            continue
            
        filepath = f"sections/{filename}"
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"# {section['title']}\n\n")
            f.write(section['content'])
        
        created_files.append(filename)
        print(f"Created: {filename}")
    
    return created_files

def main():
    # Parse the RTF-converted markdown
    sections = parse_sections('synchronism_from_rtf.md')
    
    print(f"Found {len(sections)} sections")
    
    # Create section files
    files = create_section_files(sections)
    
    print(f"\nCreated {len(files)} section files")
    
    # Update index
    with open('sections/index.md', 'w') as f:
        f.write("# Synchronism Whitepaper Sections\n\n")
        for filename in sorted(files):
            if os.path.exists(f"sections/{filename}"):
                f.write(f"- [{filename}]({filename})\n")

if __name__ == "__main__":
    main()