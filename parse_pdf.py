#!/usr/bin/env python3
"""
Parse Synchronism PDF into markdown sections
"""

import re
import os

def clean_text(text):
    """Clean up extracted text"""
    # Remove page numbers
    text = re.sub(r'^\d+\s*$', '', text, flags=re.MULTILINE)
    # Remove excessive whitespace
    text = re.sub(r'\n{3,}', '\n\n', text)
    # Fix spacing issues
    text = re.sub(r'([a-z])([A-Z])', r'\1 \2', text)
    return text.strip()

def identify_sections(text):
    """Identify major sections in the document"""
    sections = []
    
    # Common section patterns
    patterns = [
        r'^([A-Z][A-Z\s]+):',  # ALL CAPS followed by colon
        r'^(\d+\.\s+[A-Z][^.]+)',  # Numbered sections
        r'^(Chapter\s+\d+[:\s]+[^.]+)',  # Chapter headings
        r'^(Section\s+\d+[:\s]+[^.]+)',  # Section headings
        r'^(Appendix\s+[A-Z][:\s]+[^.]+)',  # Appendix headings
    ]
    
    lines = text.split('\n')
    current_section = {'title': 'Introduction', 'content': [], 'start': 0}
    
    for i, line in enumerate(lines):
        is_section = False
        for pattern in patterns:
            match = re.match(pattern, line.strip())
            if match:
                # Save previous section
                if current_section['content']:
                    current_section['content'] = '\n'.join(current_section['content'])
                    sections.append(current_section)
                
                # Start new section
                current_section = {
                    'title': match.group(1).strip(),
                    'content': [],
                    'start': i
                }
                is_section = True
                break
        
        if not is_section and line.strip():
            current_section['content'].append(line)
    
    # Add final section
    if current_section['content']:
        current_section['content'] = '\n'.join(current_section['content'])
        sections.append(current_section)
    
    return sections

def create_section_files(sections, output_dir='whitepaper/sections'):
    """Create markdown files for each section"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Create section mapping
    section_map = {
        'introduction': '00-introduction.md',
        'executive': '00-executive-summary.md',
        'glossary': '02-glossary.md',
        'unified': '03-unified-model.md',
        'intent': '04-intent-dynamics.md',
        'fractal': '05-fractal-ontology.md',
        'embryo': '06-embryogenic-cosmology.md',
        'quantum': '07-quantum-cosmic-bridge.md',
        'math': '08-mathematical-framework.md',
        'implement': '09-implementation.md',
        'govern': '10-governance-model.md',
        'conclusion': '11-conclusion.md',
        'reference': '12-references.md',
        'appendix': '13-appendices.md'
    }
    
    file_index = []
    
    for i, section in enumerate(sections):
        # Determine filename
        title_lower = section['title'].lower()
        filename = None
        
        for key, fname in section_map.items():
            if key in title_lower:
                filename = fname
                break
        
        if not filename:
            # Generate sequential filename
            filename = f"{i+1:02d}-section-{i+1}.md"
        
        filepath = os.path.join(output_dir, filename)
        
        # Write section content
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(f"# {section['title']}\n\n")
            f.write(clean_text(section['content']))
        
        file_index.append((filename, section['title']))
        print(f"Created: {filename} - {section['title'][:50]}...")
    
    return file_index

def main():
    # Read extracted text
    with open('synchronism_extracted.txt', 'r', encoding='utf-8') as f:
        text = f.read()
    
    print("Parsing Synchronism document...")
    print(f"Total characters: {len(text)}")
    
    # Identify sections
    sections = identify_sections(text)
    print(f"Found {len(sections)} sections")
    
    # Create section files
    file_index = create_section_files(sections)
    
    # Create index file
    with open('whitepaper/sections/index.md', 'w') as f:
        f.write("# Synchronism Whitepaper Sections\n\n")
        for filename, title in file_index:
            f.write(f"- [{title}]({filename})\n")
    
    print("\nSection files created in whitepaper/sections/")
    print("Index created at whitepaper/sections/index.md")

if __name__ == "__main__":
    main()