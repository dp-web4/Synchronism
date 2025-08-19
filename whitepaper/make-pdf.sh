#!/bin/bash

# make-pdf.sh - Generate PDF from Synchronism whitepaper markdown
# Usage: ./make-pdf.sh

echo "Building Synchronism whitepaper PDF..."

# Check if pandoc is installed
if ! command -v pandoc &> /dev/null; then
    echo "❌ Error: pandoc is not installed"
    echo "Please install pandoc: sudo apt-get install pandoc texlive-xetex"
    exit 1
fi

MD_FILE="build/Synchronism_Whitepaper_Complete.md"
PDF_FILE="build/Synchronism_Whitepaper.pdf"

# Ensure markdown file exists
if [ ! -f "$MD_FILE" ]; then
    echo "⚠️  Markdown file not found. Running make-md.sh first..."
    bash make-md.sh
fi

# Reorder the document to put TOC after Executive Summary
echo "Reordering document structure..."
TEMP_MD="build/Synchronism_Whitepaper_Reordered.md"

python3 << 'PYTHON_SCRIPT'
import re

# Read the complete markdown
with open('build/Synchronism_Whitepaper_Complete.md', 'r') as f:
    content = f.read()

# Split into sections
lines = content.split('\n')
sections = []
current_section = []
current_title = ""

for line in lines:
    if line.startswith('# '):
        if current_section:
            sections.append((current_title, '\n'.join(current_section)))
        current_title = line
        current_section = [line]
    else:
        current_section.append(line)

# Add the last section
if current_section:
    sections.append((current_title, '\n'.join(current_section)))

# Find Executive Summary
exec_summary = None
other_sections = []

for title, content in sections:
    if 'Executive Summary' in title:
        exec_summary = (title, content)
    else:
        other_sections.append((title, content))

# Rebuild document with custom order
with open('build/Synchronism_Whitepaper_Reordered.md', 'w') as f:
    # Check if first section already has the main title
    has_main_title = False
    for title, content in sections:
        if 'Synchronism: A Comprehensive Model' in content:
            has_main_title = True
            break
    
    if not has_main_title:
        # Add title if not present
        f.write('# Synchronism: A Comprehensive Model of Reality\n\n')
        f.write('*A Framework Unifying Scientific, Philosophical, and Spiritual Perspectives*\n\n')
        f.write('---\n\n')
    
    # Executive Summary
    if exec_summary:
        # Write executive summary as-is (it already has proper formatting)
        f.write(exec_summary[1] + '\n\n')
    
    # Add TOC marker for pandoc to generate it here
    f.write('\\newpage\n\n')
    f.write('\\tableofcontents\n\n')
    f.write('\\newpage\n\n')
    
    # All other sections
    for title, content in other_sections:
        # Skip any duplicate title sections
        if 'Synchronism: A Comprehensive Model' not in title:
            f.write(content + '\n\n')

print("✓ Document reordered with TOC after Executive Summary")
PYTHON_SCRIPT

echo "Generating PDF..."

# Generate PDF with pandoc (note: no --toc flag since we're manually placing it)
pandoc "$TEMP_MD" -o "$PDF_FILE" \
    --from markdown+raw_tex \
    --to pdf \
    --pdf-engine=xelatex \
    --toc-depth=3 \
    --highlight-style=tango \
    -V documentclass=article \
    -V geometry:margin=1in \
    -V fontsize=11pt \
    -V linkcolor=blue \
    -V urlcolor=blue \
    -V toccolor=black \
    -V colorlinks=true \
    2>/dev/null

if [ -f "$PDF_FILE" ]; then
    echo "✅ PDF created with TOC after Executive Summary: $PDF_FILE"
    echo ""
    echo "📊 PDF Statistics:"
    echo "   Size: $(du -h $PDF_FILE | cut -f1)"
    echo "   Location: $PDF_FILE"
    
    # Clean up temp file
    rm -f "$TEMP_MD"
    
    # Copy to docs for GitHub Pages access
    DOCS_DIR="../docs/whitepaper"
    if [ ! -d "$DOCS_DIR" ]; then
        mkdir -p "$DOCS_DIR"
        echo "📁 Created docs/whitepaper directory"
    fi
    
    cp "$PDF_FILE" "$DOCS_DIR/"
    echo "📄 Copied PDF to GitHub Pages location: $DOCS_DIR/Synchronism_Whitepaper.pdf"
else
    echo "❌ PDF generation failed"
    # Clean up temp file even on failure
    rm -f "$TEMP_MD"
fi