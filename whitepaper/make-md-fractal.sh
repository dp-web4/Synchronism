#!/bin/bash

# make-md-fractal.sh - Build monolithic markdown from fractal section structure
# Usage: ./make-md-fractal.sh

echo "Building Synchronism whitepaper from fractal structure..."

OUTPUT_DIR="build"
OUTPUT_FILE="$OUTPUT_DIR/Synchronism_Whitepaper_Complete.md"
SECTIONS_DIR="sections"

# Create build directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Clear output file if it exists
> "$OUTPUT_FILE"

# Function to process a directory recursively
process_directory() {
    local dir=$1
    local level=$2
    
    # Process index.md if it exists (skip it, just for navigation)
    # Process all .md files in this directory (except index.md)
    for file in "$dir"/*.md; do
        if [ -f "$file" ] && [ "$(basename "$file")" != "index.md" ]; then
            echo "$(printf '%*s' $((level*2)) '')✓ Adding $(basename "$file")"
            cat "$file" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
        fi
    done
    
    # Process subdirectories in order
    for subdir in "$dir"/*/; do
        if [ -d "$subdir" ]; then
            dirname=$(basename "$subdir")
            echo "$(printf '%*s' $((level*2)) '')📁 Processing $dirname/"
            process_directory "$subdir" $((level+1))
        fi
    done
}

# Add title page
cat > "$OUTPUT_FILE" << 'TITLE'
# Synchronism: A Comprehensive Model of Reality

**Unified Model of Reality Through Intent Dynamics**

Version: 0.24.09.28.11.00

---

TITLE

echo "Processing fractal section structure..."
echo ""

# Process the sections directory recursively
process_directory "$SECTIONS_DIR" 0

# Add timestamp
echo "" >> "$OUTPUT_FILE"
echo "---" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "*Generated: $(date '+%Y-%m-%d %H:%M:%S')*" >> "$OUTPUT_FILE"
echo "*Source: Synchronism RTF Document*" >> "$OUTPUT_FILE"
echo "*Latest: https://dpcars.net/Synchronism_0.pdf*" >> "$OUTPUT_FILE"

echo ""
echo "✅ Monolithic markdown created: $OUTPUT_FILE"
echo ""

# Show file info
if [ -f "$OUTPUT_FILE" ]; then
    lines=$(wc -l < "$OUTPUT_FILE")
    size=$(du -h "$OUTPUT_FILE" | cut -f1)
    echo "📊 Statistics:"
    echo "   Lines: $lines"
    echo "   Size: $size"
    
    # Count sections
    sections=$(find "$SECTIONS_DIR" -type d | wc -l)
    files=$(find "$SECTIONS_DIR" -name "*.md" -not -name "index.md" | wc -l)
    echo "   Sections: $sections directories"
    echo "   Documents: $files markdown files"
    
    # Copy to docs for GitHub Pages access
    DOCS_DIR="../docs/whitepaper"
    if [ ! -d "$DOCS_DIR" ]; then
        mkdir -p "$DOCS_DIR"
        echo ""
        echo "📁 Created docs/whitepaper directory"
    fi
    
    cp "$OUTPUT_FILE" "$DOCS_DIR/"
    echo ""
    echo "📄 Copied to GitHub Pages: $DOCS_DIR/$(basename "$OUTPUT_FILE")"
fi