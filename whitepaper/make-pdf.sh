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

echo "Generating PDF..."

# Generate PDF with pandoc
pandoc "$MD_FILE" -o "$PDF_FILE" \
    --from markdown \
    --to pdf \
    --pdf-engine=xelatex \
    --toc \
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
    echo "✅ PDF created: $PDF_FILE"
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