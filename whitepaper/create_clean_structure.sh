#!/bin/bash

# Create clean whitepaper structure based on web-version organization
# This creates the directory structure with placeholder files

echo "Creating clean whitepaper structure..."

# Base sections directory
SECTIONS="sections"

# Create main sections
mkdir -p "$SECTIONS/00-executive-summary"
mkdir -p "$SECTIONS/01-introduction"
mkdir -p "$SECTIONS/02-perspective"
mkdir -p "$SECTIONS/03-hermetic-principles"

# 04-fundamental-concepts with subdirectories
mkdir -p "$SECTIONS/04-fundamental-concepts/01-universe-grid"
mkdir -p "$SECTIONS/04-fundamental-concepts/02-time-slices"
mkdir -p "$SECTIONS/04-fundamental-concepts/03-intent-transfer"
mkdir -p "$SECTIONS/04-fundamental-concepts/04-emergence"
mkdir -p "$SECTIONS/04-fundamental-concepts/05-field-effects"
mkdir -p "$SECTIONS/04-fundamental-concepts/06-interaction-modes"
mkdir -p "$SECTIONS/04-fundamental-concepts/07-coherence"
mkdir -p "$SECTIONS/04-fundamental-concepts/08-markov-blankets"
mkdir -p "$SECTIONS/04-fundamental-concepts/09-markov-relevancy"
mkdir -p "$SECTIONS/04-fundamental-concepts/10-spectral-existence"
mkdir -p "$SECTIONS/04-fundamental-concepts/11-abstraction"
mkdir -p "$SECTIONS/04-fundamental-concepts/12-entity-interactions"

# 05-quantum-macro with subdirectories
mkdir -p "$SECTIONS/05-quantum-macro/01-crt-analogy"
mkdir -p "$SECTIONS/05-quantum-macro/02-superposition"
mkdir -p "$SECTIONS/05-quantum-macro/03-wave-particle"
mkdir -p "$SECTIONS/05-quantum-macro/04-entanglement"
mkdir -p "$SECTIONS/05-quantum-macro/05-witness-effect"
mkdir -p "$SECTIONS/05-quantum-macro/06-relativity"
mkdir -p "$SECTIONS/05-quantum-macro/07-speed-limits"
mkdir -p "$SECTIONS/05-quantum-macro/08-decoherence"
mkdir -p "$SECTIONS/05-quantum-macro/09-temperature"
mkdir -p "$SECTIONS/05-quantum-macro/10-energy"
mkdir -p "$SECTIONS/05-quantum-macro/11-universal-field"
mkdir -p "$SECTIONS/05-quantum-macro/12-chemistry"
mkdir -p "$SECTIONS/05-quantum-macro/13-life-cognition"
mkdir -p "$SECTIONS/05-quantum-macro/14-gravity"
mkdir -p "$SECTIONS/05-quantum-macro/15-dark-matter"
mkdir -p "$SECTIONS/05-quantum-macro/16-superconductivity"
mkdir -p "$SECTIONS/05-quantum-macro/17-permeability"
mkdir -p "$SECTIONS/05-quantum-macro/18-electromagnetic"
mkdir -p "$SECTIONS/05-quantum-macro/19-energy-refinement"
mkdir -p "$SECTIONS/05-quantum-macro/20-temperature-refinement"
mkdir -p "$SECTIONS/05-quantum-macro/21-cognition-refinement"
mkdir -p "$SECTIONS/05-quantum-macro/22-string-theory"

# 06-implications with subdirectories
mkdir -p "$SECTIONS/06-implications/01-unified-understanding"
mkdir -p "$SECTIONS/06-implications/02-scientific-inquiry"
mkdir -p "$SECTIONS/06-implications/03-ethical-philosophical"
mkdir -p "$SECTIONS/06-implications/04-open-questions"

# Simple sections
mkdir -p "$SECTIONS/07-conclusion"
mkdir -p "$SECTIONS/08-glossary"
mkdir -p "$SECTIONS/09-appendix-mathematical"

echo "  ✓ Created directory structure"

# Create placeholder files
echo "Creating placeholder markdown files..."

# Executive Summary (from archived whitepaper)
cat > "$SECTIONS/00-executive-summary/executive-summary.md" << 'EOF'
# Executive Summary

Synchronism presents a comprehensive model of reality through the lens of intent dynamics. This framework bridges scientific understanding with philosophical inquiry by proposing that all observable phenomena emerge from fundamental patterns of "intent" transfer through a discrete spacetime grid.

## Core Principles

1. **Intent as Fundamental Force**: All interactions derive from quantized intent transfer
2. **Single Observer Model**: Reality experienced through unified consciousness
3. **Pattern-Based Existence**: Everything exists as cycling patterns, never static states
4. **Synchronization vs Interaction**: Observation synchronizes with patterns rather than changing them
5. **Scale-Invariant Coherence**: Principles apply uniformly from quantum to cosmic scales

## Key Innovations

- **Markov Blankets**: Define boundaries between scales of existence
- **Spectral Existence**: Entities exist on a spectrum of coherence
- **Witness Effect**: Participatory observation without measurement collapse
- **Coherence-Based Ethics**: Morality as optimization of pattern coherence

This document explores these concepts through mathematical formalism, practical applications, and philosophical implications.
EOF

# Create content extraction script
cat > "$SECTIONS/extract_web_content.py" << 'EOF'
#!/usr/bin/env python3
"""
Extract content from web-version HTML files into markdown
Run this to populate the clean structure with actual content
"""

import os
import re
from pathlib import Path

def extract_html_section(html_file):
    """Extract main content from HTML file"""
    try:
        with open(html_file, 'r', encoding='utf-8') as f:
            content = f.read()
            
        # Extract content between section tags
        match = re.search(r'<section[^>]*>(.*?)</section>', content, re.DOTALL)
        if match:
            html_content = match.group(1)
            
            # Convert HTML to Markdown (basic conversion)
            md = html_content
            
            # Headers
            md = re.sub(r'<h1[^>]*>(.*?)</h1>', r'# \1\n', md)
            md = re.sub(r'<h2[^>]*>(.*?)</h2>', r'## \1\n', md)
            md = re.sub(r'<h3[^>]*>(.*?)</h3>', r'### \1\n', md)
            md = re.sub(r'<h4[^>]*>(.*?)</h4>', r'#### \1\n', md)
            
            # Paragraphs
            md = re.sub(r'<p[^>]*>(.*?)</p>', r'\1\n\n', md)
            
            # Lists
            md = re.sub(r'<li[^>]*>(.*?)</li>', r'- \1\n', md)
            md = re.sub(r'<ul[^>]*>', '', md)
            md = re.sub(r'</ul>', '\n', md)
            
            # Links
            md = re.sub(r'<a href="([^"]*)"[^>]*>(.*?)</a>', r'[\2](\1)', md)
            
            # Bold/Italic
            md = re.sub(r'<strong[^>]*>(.*?)</strong>', r'**\1**', md)
            md = re.sub(r'<b[^>]*>(.*?)</b>', r'**\1**', md)
            md = re.sub(r'<em[^>]*>(.*?)</em>', r'*\1*', md)
            md = re.sub(r'<i[^>]*>(.*?)</i>', r'*\1*', md)
            
            # Code
            md = re.sub(r'<code[^>]*>(.*?)</code>', r'`\1`', md)
            
            # Remove remaining HTML tags
            md = re.sub(r'<[^>]+>', '', md)
            
            # Clean up whitespace
            md = re.sub(r'\n\s*\n\s*\n', '\n\n', md)
            
            return md.strip()
    except Exception as e:
        print(f"Error extracting {html_file}: {e}")
    return None

# Usage: python3 extract_web_content.py
EOF

chmod +x "$SECTIONS/extract_web_content.py"

# Create simple content files as placeholders
for dir in "$SECTIONS"/04-fundamental-concepts/*; do
    dirname=$(basename "$dir")
    clean_name=$(echo "$dirname" | sed 's/-/ /g' | sed 's/\b\(.\)/\u\1/g')
    cat > "$dir/content.md" << EOF
# $clean_name

*Content to be extracted from web-version and merged with archived content*

## Key Concepts

- Intent transfer mechanisms
- Pattern coherence
- Scale boundaries

## Mathematical Framework

See Appendix A for detailed formalism.
EOF
done

for dir in "$SECTIONS"/05-quantum-macro/*; do
    dirname=$(basename "$dir")
    clean_name=$(echo "$dirname" | sed 's/-/ /g' | sed 's/\b\(.\)/\u\1/g')
    cat > "$dir/content.md" << EOF
# $clean_name

*Content to be extracted from web-version and merged with archived content*

## Synchronism Perspective

- Pattern-based interpretation
- No wave function collapse
- Synchronization rather than measurement

## Applications

Practical implications for understanding reality.
EOF
done

echo "  ✓ Created placeholder files"

echo ""
echo "Clean structure created successfully!"
echo ""
echo "Directory structure:"
tree -d "$SECTIONS" -L 2

echo ""
echo "Next steps:"
echo "1. Run extract_web_content.py to populate with web-version content"
echo "2. Selectively merge useful content from sections_archived_chaotic/"
echo "3. Test build scripts with new clean structure"
echo "4. Review and refine content section by section"