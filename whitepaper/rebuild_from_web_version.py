#!/usr/bin/env python3
"""
Rebuild clean whitepaper structure from web-version HTML files
"""

import os
import re
from pathlib import Path
from html.parser import HTMLParser

class HTMLToMarkdown(HTMLParser):
    """Simple HTML to Markdown converter"""
    def __init__(self):
        super().__init__()
        self.markdown = []
        self.in_heading = 0
        self.in_paragraph = False
        self.in_list = False
        self.in_code = False
        self.skip_script = False
        self.skip_style = False
        
    def handle_starttag(self, tag, attrs):
        if tag == 'script':
            self.skip_script = True
        elif tag == 'style':
            self.skip_style = True
        elif tag == 'h1':
            self.in_heading = 1
            self.markdown.append('\n# ')
        elif tag == 'h2':
            self.in_heading = 2
            self.markdown.append('\n## ')
        elif tag == 'h3':
            self.in_heading = 3
            self.markdown.append('\n### ')
        elif tag == 'h4':
            self.in_heading = 4
            self.markdown.append('\n#### ')
        elif tag == 'p':
            self.in_paragraph = True
            self.markdown.append('\n')
        elif tag == 'br':
            self.markdown.append('\n')
        elif tag == 'ul' or tag == 'ol':
            self.in_list = True
            self.markdown.append('\n')
        elif tag == 'li':
            self.markdown.append('- ')
        elif tag == 'code':
            self.in_code = True
            self.markdown.append('`')
        elif tag == 'pre':
            self.markdown.append('\n```\n')
        elif tag == 'blockquote':
            self.markdown.append('\n> ')
        elif tag == 'em' or tag == 'i':
            self.markdown.append('*')
        elif tag == 'strong' or tag == 'b':
            self.markdown.append('**')
            
    def handle_endtag(self, tag):
        if tag == 'script':
            self.skip_script = False
        elif tag == 'style':
            self.skip_style = False
        elif tag in ['h1', 'h2', 'h3', 'h4']:
            self.in_heading = 0
            self.markdown.append('\n')
        elif tag == 'p':
            self.in_paragraph = False
            self.markdown.append('\n')
        elif tag == 'ul' or tag == 'ol':
            self.in_list = False
        elif tag == 'li':
            self.markdown.append('\n')
        elif tag == 'code':
            self.in_code = False
            self.markdown.append('`')
        elif tag == 'pre':
            self.markdown.append('\n```\n')
        elif tag == 'em' or tag == 'i':
            self.markdown.append('*')
        elif tag == 'strong' or tag == 'b':
            self.markdown.append('**')
            
    def handle_data(self, data):
        if not self.skip_script and not self.skip_style:
            # Clean up whitespace
            if self.in_heading or self.in_paragraph or self.in_list or self.in_code:
                data = data.strip()
                if data:
                    self.markdown.append(data)
                    
    def get_markdown(self):
        return ''.join(self.markdown).strip()

def extract_html_content(html_file):
    """Extract and convert HTML content to markdown"""
    try:
        with open(html_file, 'r', encoding='utf-8') as f:
            html_content = f.read()
            
        # Remove navigation and other non-content elements
        # Look for main content area
        content_match = re.search(r'<main[^>]*>(.*?)</main>', html_content, re.DOTALL)
        if not content_match:
            content_match = re.search(r'<article[^>]*>(.*?)</article>', html_content, re.DOTALL)
        if not content_match:
            content_match = re.search(r'<div[^>]*class="content"[^>]*>(.*?)</div>', html_content, re.DOTALL)
        if not content_match:
            # Fall back to body content
            content_match = re.search(r'<body[^>]*>(.*?)</body>', html_content, re.DOTALL)
            
        if content_match:
            content = content_match.group(1)
            # Parse HTML to Markdown
            parser = HTMLToMarkdown()
            parser.feed(content)
            return parser.get_markdown()
    except Exception as e:
        print(f"Error processing {html_file}: {e}")
    return None

def create_clean_structure():
    """Create clean directory structure based on web-version"""
    
    web_sections = "/mnt/c/projects/ai-agents/Synchronism/web-version/sections"
    new_sections = "/mnt/c/projects/ai-agents/Synchronism/whitepaper/sections"
    
    # Define the clean structure based on web-version
    structure = {
        "01-introduction": {
            "files": ["introduction.md", "about.md"],
            "source": "01-introduction/index.html"
        },
        "02-perspective": {
            "files": ["perspective.md"],
            "source": "02-perspective/index.html"
        },
        "03-hermetic-principles": {
            "files": ["hermetic-principles.md"],
            "source": "03-hermetic-principles/index.html"
        },
        "04-fundamental-concepts": {
            "subdirs": {
                "01-universe-grid": {"source": "04-fundamental-concepts/01-universe-grid/index.html"},
                "02-time-slices": {"source": "04-fundamental-concepts/02-time-slices/index.html"},
                "03-intent-transfer": {"source": "04-fundamental-concepts/03-intent-transfer/index.html"},
                "04-emergence": {"source": "04-fundamental-concepts/04-emergence/index.html"},
                "05-field-effects": {"source": "04-fundamental-concepts/05-field-effects/index.html"},
                "06-interaction-modes": {"source": "04-fundamental-concepts/06-interaction-modes/index.html"},
                "07-coherence": {"source": "04-fundamental-concepts/07-coherence/index.html"},
                "08-markov-blankets": {"source": "04-fundamental-concepts/08-markov-blankets/index.html"},
                "09-markov-relevancy": {"source": "04-fundamental-concepts/09-mrh/index.html"},
                "10-spectral-existence": {"source": "04-fundamental-concepts/10-spectral-existence/index.html"},
                "11-abstraction": {"source": "04-fundamental-concepts/11-abstraction/index.html"},
                "12-entity-interactions": {"source": "04-fundamental-concepts/12-entity-interactions/index.html"}
            }
        },
        "05-quantum-macro": {
            "subdirs": {
                "01-crt-analogy": {"source": "05-quantum-macro/01-crt-analogy/index.html"},
                "02-superposition": {"source": "05-quantum-macro/02-superposition/index.html"},
                "03-wave-particle": {"source": "05-quantum-macro/03-wave-particle/index.html"},
                "04-entanglement": {"source": "05-quantum-macro/04-entanglement/index.html"},
                "05-witness-effect": {"source": "05-quantum-macro/05-witness-effect/index.html"},
                "06-relativity": {"source": "05-quantum-macro/06-relativity/index.html"},
                "07-speed-limits": {"source": "05-quantum-macro/07-speed-limits/index.html"},
                "08-decoherence": {"source": "05-quantum-macro/08-decoherence/index.html"},
                "09-temperature": {"source": "05-quantum-macro/09-temperature/index.html"},
                "10-energy": {"source": "05-quantum-macro/10-energy/index.html"},
                "11-universal-field": {"source": "05-quantum-macro/11-universal-field/index.html"},
                "12-chemistry": {"source": "05-quantum-macro/12-chemistry/index.html"},
                "13-life-cognition": {"source": "05-quantum-macro/13-life-cognition/index.html"},
                "14-gravity": {"source": "05-quantum-macro/14-gravity/index.html"},
                "15-dark-matter": {"source": "05-quantum-macro/15-dark-matter/index.html"},
                "16-superconductivity": {"source": "05-quantum-macro/16-superconductivity/index.html"},
                "17-permeability": {"source": "05-quantum-macro/17-permeability/index.html"},
                "18-electromagnetic": {"source": "05-quantum-macro/18-electromagnetic/index.html"},
                "19-energy-refinement": {"source": "05-quantum-macro/19-energy-refinement/index.html"},
                "20-temperature-refinement": {"source": "05-quantum-macro/20-temperature-refinement/index.html"},
                "21-cognition-refinement": {"source": "05-quantum-macro/21-cognition-refinement/index.html"},
                "22-string-theory": {"source": "05-quantum-macro/22-string-theory/index.html"}
            }
        },
        "06-implications": {
            "subdirs": {
                "01-unified-understanding": {"source": "06-implications/01-unified-understanding/index.html"},
                "02-scientific-inquiry": {"source": "06-implications/02-scientific-inquiry/index.html"},
                "03-ethical-philosophical": {"source": "06-implications/03-ethical-philosophical/index.html"},
                "04-open-questions": {"source": "06-implications/04-open-questions/index.html"}
            }
        },
        "07-conclusion": {
            "files": ["conclusion.md"],
            "source": "07-conclusion/index.html"
        },
        "08-glossary": {
            "files": ["glossary.md"],
            "source": "08-glossary/index.html"
        },
        "09-appendix-mathematical": {
            "files": ["mathematical-framework.md"],
            "source": "09-appendix-mathematical/index.html"
        }
    }
    
    # Process each section
    for section_name, section_info in structure.items():
        section_path = os.path.join(new_sections, section_name)
        os.makedirs(section_path, exist_ok=True)
        
        # Handle direct files
        if "files" in section_info and "source" in section_info:
            source_file = os.path.join(web_sections, section_info["source"])
            if os.path.exists(source_file):
                content = extract_html_content(source_file)
                if content:
                    # Save as first file in list
                    output_file = os.path.join(section_path, section_info["files"][0])
                    with open(output_file, 'w', encoding='utf-8') as f:
                        f.write(content)
                    print(f"Created: {section_name}/{section_info['files'][0]}")
        
        # Handle subdirectories
        if "subdirs" in section_info:
            for subdir_name, subdir_info in section_info["subdirs"].items():
                subdir_path = os.path.join(section_path, subdir_name)
                os.makedirs(subdir_path, exist_ok=True)
                
                if "source" in subdir_info:
                    source_file = os.path.join(web_sections, subdir_info["source"])
                    if os.path.exists(source_file):
                        content = extract_html_content(source_file)
                        if content:
                            # Create markdown file with clean name
                            filename = subdir_name.replace('-', '_') + '.md'
                            output_file = os.path.join(subdir_path, filename)
                            with open(output_file, 'w', encoding='utf-8') as f:
                                f.write(content)
                            print(f"Created: {section_name}/{subdir_name}/{filename}")
                else:
                    # Create placeholder
                    output_file = os.path.join(subdir_path, "content.md")
                    with open(output_file, 'w', encoding='utf-8') as f:
                        f.write(f"# {subdir_name.replace('-', ' ').title()}\n\n*Content to be extracted and reviewed*\n")
        
        # Create index.md for navigation
        index_file = os.path.join(section_path, "index.md")
        with open(index_file, 'w', encoding='utf-8') as f:
            f.write(f"# {section_name.replace('-', ' ').title()}\n\n")
            if "subdirs" in section_info:
                f.write("## Sections\n\n")
                for subdir in section_info.get("subdirs", {}).keys():
                    clean_name = subdir.replace('-', ' ').title()
                    f.write(f"- [{clean_name}]({subdir}/)\n")

def main():
    print("Rebuilding clean whitepaper structure from web-version...")
    print("=" * 60)
    
    # Create the clean structure
    create_clean_structure()
    
    print("\n" + "=" * 60)
    print("Clean structure created!")
    print("\nNext steps:")
    print("1. Review extracted content for formatting issues")
    print("2. Merge any useful content from archived chaotic structure")
    print("3. Update build scripts to use new structure")
    print("4. Test PDF, MD, and Web generation")

if __name__ == "__main__":
    main()