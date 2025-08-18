#!/usr/bin/env python3
"""
Extract content from web-version HTML files into markdown
"""

import os
import re
from pathlib import Path

def extract_html_section(html_file):
    """Extract content from HTML file and convert to markdown"""
    try:
        with open(html_file, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Extract content between section tags
        match = re.search(r'<section[^>]*>(.*?)</section>', content, re.DOTALL)
        if not match:
            print(f"No section found in {html_file}")
            return None
            
        html_content = match.group(1)
        
        # Convert HTML to Markdown
        md = html_content
        
        # Headers
        md = re.sub(r'<h1[^>]*>(.*?)</h1>', r'# \1\n', md)
        md = re.sub(r'<h2[^>]*>(.*?)</h2>', r'## \1\n', md)
        md = re.sub(r'<h3[^>]*>(.*?)</h3>', r'### \1\n', md)
        md = re.sub(r'<h4[^>]*>(.*?)</h4>', r'#### \1\n', md)
        
        # Paragraphs and divs
        md = re.sub(r'<p[^>]*>(.*?)</p>', r'\1\n\n', md, flags=re.DOTALL)
        md = re.sub(r'<div class="section-content"[^>]*>(.*?)</div>', r'\1', md, flags=re.DOTALL)
        md = re.sub(r'<div class="key-concept"[^>]*>(.*?)</div>', r'\1', md, flags=re.DOTALL)
        md = re.sub(r'<div class="navigation-hints"[^>]*>(.*?)</div>', r'\n---\n\1\n---\n', md, flags=re.DOTALL)
        md = re.sub(r'<div[^>]*>', '', md)
        md = re.sub(r'</div>', '', md)
        
        # Lists
        md = re.sub(r'<ul[^>]*>(.*?)</ul>', lambda m: process_list(m.group(1)), md, flags=re.DOTALL)
        md = re.sub(r'<ol[^>]*>(.*?)</ol>', lambda m: process_list(m.group(1), ordered=True), md, flags=re.DOTALL)
        
        # Links
        md = re.sub(r'<a href="([^"]*)"[^>]*>(.*?)</a>', r'[\2](\1)', md)
        
        # Bold/Italic
        md = re.sub(r'<strong[^>]*>(.*?)</strong>', r'**\1**', md)
        md = re.sub(r'<b[^>]*>(.*?)</b>', r'**\1**', md)
        md = re.sub(r'<em[^>]*>(.*?)</em>', r'*\1*', md)
        md = re.sub(r'<i[^>]*>(.*?)</i>', r'*\1*', md)
        
        # Code
        md = re.sub(r'<code[^>]*>(.*?)</code>', r'`\1`', md)
        md = re.sub(r'<pre[^>]*>(.*?)</pre>', r'```\n\1\n```', md, flags=re.DOTALL)
        
        # Blockquotes
        md = re.sub(r'<blockquote[^>]*>(.*?)</blockquote>', r'> \1', md, flags=re.DOTALL)
        
        # Remove remaining HTML tags
        md = re.sub(r'<[^>]+>', '', md)
        
        # Clean up whitespace
        md = re.sub(r'\n\s*\n\s*\n', '\n\n', md)
        md = re.sub(r' +', ' ', md)
        
        # HTML entities
        md = md.replace('&lt;', '<')
        md = md.replace('&gt;', '>')
        md = md.replace('&amp;', '&')
        md = md.replace('&quot;', '"')
        md = md.replace('&#39;', "'")
        
        return md.strip()
    except Exception as e:
        print(f"Error extracting {html_file}: {e}")
    return None

def process_list(html_list, ordered=False):
    """Process HTML list items"""
    items = re.findall(r'<li[^>]*>(.*?)</li>', html_list, re.DOTALL)
    result = []
    for i, item in enumerate(items):
        item = item.strip()
        if ordered:
            result.append(f"{i+1}. {item}")
        else:
            result.append(f"- {item}")
    return '\n'.join(result) + '\n'

def main():
    """Main extraction process"""
    web_sections = Path("/mnt/c/projects/ai-agents/Synchronism/web-version/sections")
    whitepaper_sections = Path("/mnt/c/projects/ai-agents/Synchronism/whitepaper/sections")
    
    if not web_sections.exists():
        print(f"Error: Web sections directory not found: {web_sections}")
        return
    
    # Define extraction mapping
    extraction_map = [
        # Simple sections
        ("01-introduction", "01-introduction/index.html", "introduction.md"),
        ("02-perspective", "02-perspective/index.html", "perspective.md"),
        ("03-hermetic-principles", "03-hermetic-principles/index.html", "hermetic-principles.md"),
        ("07-conclusion", "07-conclusion/index.html", "conclusion.md"),
        ("08-glossary", "08-glossary/index.html", "glossary.md"),
        
        # Fundamental concepts subsections
        ("04-fundamental-concepts/01-universe-grid", "04-fundamental-concepts/01-universe-grid/index.html", "universe_grid.md"),
        ("04-fundamental-concepts/02-time-slices", "04-fundamental-concepts/02-time-slices/index.html", "time_slices.md"),
        ("04-fundamental-concepts/03-intent-transfer", "04-fundamental-concepts/03-intent-transfer/index.html", "intent_transfer.md"),
        ("04-fundamental-concepts/04-emergence", "04-fundamental-concepts/04-emergence/index.html", "emergence.md"),
        ("04-fundamental-concepts/05-field-effects", "04-fundamental-concepts/05-field-effects/index.html", "field_effects.md"),
        ("04-fundamental-concepts/06-interaction-modes", "04-fundamental-concepts/06-interaction-modes/index.html", "interaction_modes.md"),
        ("04-fundamental-concepts/07-coherence", "04-fundamental-concepts/07-coherence/index.html", "coherence.md"),
        ("04-fundamental-concepts/08-markov-blankets", "04-fundamental-concepts/08-markov-blankets/index.html", "markov_blankets.md"),
        ("04-fundamental-concepts/09-markov-relevancy", "04-fundamental-concepts/09-mrh/index.html", "markov_relevancy.md"),
        ("04-fundamental-concepts/10-spectral-existence", "04-fundamental-concepts/10-spectral-existence/index.html", "spectral_existence.md"),
        ("04-fundamental-concepts/11-abstraction", "04-fundamental-concepts/11-abstraction/index.html", "abstraction.md"),
        ("04-fundamental-concepts/12-entity-interactions", "04-fundamental-concepts/12-entity-interactions/index.html", "entity_interactions.md"),
        
        # Quantum-macro subsections (first batch)
        ("05-quantum-macro/01-crt-analogy", "05-quantum-macro/01-crt-analogy/index.html", "crt_analogy.md"),
        ("05-quantum-macro/02-superposition", "05-quantum-macro/02-superposition/index.html", "superposition.md"),
        ("05-quantum-macro/03-wave-particle", "05-quantum-macro/03-wave-particle/index.html", "wave_particle.md"),
        ("05-quantum-macro/04-entanglement", "05-quantum-macro/04-entanglement/index.html", "entanglement.md"),
        ("05-quantum-macro/05-witness-effect", "05-quantum-macro/05-witness-effect/index.html", "witness_effect.md"),
        ("05-quantum-macro/06-relativity", "05-quantum-macro/06-relativity/index.html", "relativity.md"),
        ("05-quantum-macro/07-speed-limits", "05-quantum-macro/07-speed-limits/index.html", "speed_limits.md"),
        ("05-quantum-macro/08-decoherence", "05-quantum-macro/08-decoherence/index.html", "decoherence.md"),
        ("05-quantum-macro/09-temperature", "05-quantum-macro/09-temperature/index.html", "temperature.md"),
        ("05-quantum-macro/10-energy", "05-quantum-macro/10-energy/index.html", "energy.md"),
        
        # Quantum-macro subsections (second batch)
        ("05-quantum-macro/11-universal-field", "05-quantum-macro/11-universal-field/index.html", "universal_field.md"),
        ("05-quantum-macro/12-chemistry", "05-quantum-macro/12-chemistry/index.html", "chemistry.md"),
        ("05-quantum-macro/13-life-cognition", "05-quantum-macro/13-life-cognition/index.html", "life_cognition.md"),
        ("05-quantum-macro/14-gravity", "05-quantum-macro/14-gravity/index.html", "gravity.md"),
        ("05-quantum-macro/15-dark-matter", "05-quantum-macro/15-dark-matter/index.html", "dark_matter.md"),
        ("05-quantum-macro/16-superconductivity", "05-quantum-macro/16-superconductivity/index.html", "superconductivity.md"),
        ("05-quantum-macro/17-permeability", "05-quantum-macro/17-permeability/index.html", "permeability.md"),
        ("05-quantum-macro/18-electromagnetic", "05-quantum-macro/18-electromagnetic/index.html", "electromagnetic.md"),
        ("05-quantum-macro/19-energy-refinement", "05-quantum-macro/19-energy-refinement/index.html", "energy_refinement.md"),
        ("05-quantum-macro/20-temperature-refinement", "05-quantum-macro/20-temperature-refinement/index.html", "temperature_refinement.md"),
        ("05-quantum-macro/21-cognition-refinement", "05-quantum-macro/21-cognition-refinement/index.html", "cognition_refinement.md"),
        ("05-quantum-macro/22-string-theory", "05-quantum-macro/22-string-theory/index.html", "string_theory.md"),
        
        # Implications subsections
        ("06-implications/01-unified-understanding", "06-implications/01-unified-understanding/index.html", "unified_understanding.md"),
        ("06-implications/02-scientific-inquiry", "06-implications/02-scientific-inquiry/index.html", "scientific_inquiry.md"),
        ("06-implications/03-ethical-philosophical", "06-implications/03-ethical-philosophical/index.html", "ethical_philosophical.md"),
        ("06-implications/04-open-questions", "06-implications/04-open-questions/index.html", "open_questions.md"),
        
        # Appendix sections
        ("09-appendix-mathematical", "09-appendix-mathematical/index.html", "mathematical_framework.md"),
    ]
    
    success_count = 0
    fail_count = 0
    
    print("Extracting content from web-version HTML files...")
    print("=" * 60)
    
    for output_dir, source_file, output_name in extraction_map:
        source_path = web_sections / source_file
        output_path = whitepaper_sections / output_dir / output_name
        
        # Create directory if needed
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        if source_path.exists():
            content = extract_html_section(source_path)
            if content:
                with open(output_path, 'w', encoding='utf-8') as f:
                    f.write(content)
                print(f"✓ Extracted: {output_dir}/{output_name}")
                success_count += 1
            else:
                print(f"✗ Failed to extract: {source_file}")
                fail_count += 1
        else:
            print(f"✗ Source not found: {source_file}")
            fail_count += 1
    
    # Create index files for main sections
    print("\nCreating section index files...")
    create_section_indexes(whitepaper_sections)
    
    print("\n" + "=" * 60)
    print(f"Extraction complete: {success_count} successful, {fail_count} failed")
    print("\nNext steps:")
    print("1. Review extracted content for formatting issues")
    print("2. Update build scripts to use new structure")
    print("3. Test PDF, MD, and Web generation")

def create_section_indexes(sections_dir):
    """Create index.md files for navigation"""
    indexes = {
        "04-fundamental-concepts": [
            "01-universe-grid", "02-time-slices", "03-intent-transfer",
            "04-emergence", "05-field-effects", "06-interaction-modes",
            "07-coherence", "08-markov-blankets", "09-markov-relevancy",
            "10-spectral-existence", "11-abstraction", "12-entity-interactions"
        ],
        "05-quantum-macro": [
            "01-crt-analogy", "02-superposition", "03-wave-particle",
            "04-entanglement", "05-witness-effect", "06-relativity",
            "07-speed-limits", "08-decoherence", "09-temperature",
            "10-energy", "11-universal-field", "12-chemistry",
            "13-life-cognition", "14-gravity", "15-dark-matter",
            "16-superconductivity", "17-permeability", "18-electromagnetic",
            "19-energy-refinement", "20-temperature-refinement",
            "21-cognition-refinement", "22-string-theory"
        ],
        "06-implications": [
            "01-unified-understanding", "02-scientific-inquiry",
            "03-ethical-philosophical", "04-open-questions"
        ]
    }
    
    for section, subdirs in indexes.items():
        index_path = sections_dir / section / "index.md"
        with open(index_path, 'w', encoding='utf-8') as f:
            title = section.replace('-', ' ').title()
            f.write(f"# {title}\n\n")
            f.write("## Sections\n\n")
            for subdir in subdirs:
                clean_name = subdir[3:].replace('-', ' ').title() if subdir[:2].isdigit() else subdir.replace('-', ' ').title()
                f.write(f"- [{clean_name}]({subdir}/)\n")
        print(f"✓ Created index: {section}/index.md")

if __name__ == "__main__":
    main()