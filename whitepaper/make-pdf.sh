#!/bin/bash

# make-pdf.sh - Generate PDF from Synchronism whitepaper markdown
# Usage: ./make-pdf.sh

# Pull latest changes before building to avoid conflicts
echo "Checking for updates..."
git fetch

# Check if we're behind the remote
LOCAL=$(git rev-parse @)
REMOTE=$(git rev-parse @{u})
BASE=$(git merge-base @ @{u})

if [ $LOCAL = $REMOTE ]; then
    echo "Already up to date."
elif [ $LOCAL = $BASE ]; then
    echo "Pulling latest changes..."
    git pull
else
    echo "❌ Error: Your branch has diverged from the remote branch."
    echo "Please resolve conflicts manually before building:"
    echo "  1. Review changes with: git status"
    echo "  2. Either stash your changes: git stash"
    echo "  3. Or commit them: git add . && git commit -m 'your message'"
    echo "  4. Then pull: git pull"
    echo "  5. Run this script again"
    exit 1
fi

# Preprocess proposals to demote headers
./preprocess-proposals.sh
./preprocess-sections.sh

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

# Split by Executive Summary to preserve all content
lines = content.split('\n')
before_exec = []
exec_summary_content = []
after_exec = []
current_section = 'before'

for line in lines:
    if line == '# Executive Summary':
        current_section = 'exec'
        exec_summary_content.append(line)
    elif current_section == 'exec' and line.startswith('## ') and 'Synchronism' not in line:
        # Found start of next section after exec summary
        current_section = 'after'
        after_exec.append(line)
    elif current_section == 'before':
        before_exec.append(line)
    elif current_section == 'exec':
        exec_summary_content.append(line)
    else:
        after_exec.append(line)

# Rebuild document with custom order for TOC placement
with open('build/Synchronism_Whitepaper_Reordered.md', 'w') as f:
    # Start with title page (from before_exec content)
    # Skip any empty lines at the start
    title_written = False
    for line in before_exec:
        if line.strip():
            if not title_written and line.startswith('# Synchronism'):
                f.write(line + '\n\n')
                f.write('*A Framework Unifying Scientific, Philosophical, and Spiritual Perspectives*\n\n')
                f.write('---\n\n')
                title_written = True
            elif title_written:
                break  # Don't include other content before exec summary
    
    # Executive Summary (change # to ## so it's not a chapter)
    if exec_summary_content:
        # Replace the # with ## for executive summary
        for line in exec_summary_content:
            if line == '# Executive Summary':
                f.write('## Executive Summary\n')
            else:
                f.write(line + '\n')
        f.write('\n')
    
    # Add TOC after Executive Summary (no newpage before, just after)
    f.write('\\tableofcontents\n')
    f.write('\\newpage\n')
    
    # All content after Executive Summary
    f.write('\n'.join(after_exec) + '\n')

print("✓ Document reordered with TOC after Executive Summary")
PYTHON_SCRIPT

echo "Generating PDF..."

# Generate PDF with pandoc
# Using markdown-compact_definition_lists to ensure proper list spacing
pandoc "$TEMP_MD" -o "$PDF_FILE" \
    --from markdown+raw_tex-compact_definition_lists \
    --to pdf \
    --pdf-engine=xelatex \
    --toc-depth=2 \
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
fi

# Always clean up temp file
rm -f "$TEMP_MD"