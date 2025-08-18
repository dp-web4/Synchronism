#!/bin/bash

# make-md.sh - Combine Synchronism whitepaper sections into monolithic markdown
# Usage: ./make-md.sh

echo "Building monolithic Synchronism whitepaper markdown..."

OUTPUT_DIR="build"
OUTPUT_FILE="$OUTPUT_DIR/Synchronism_Whitepaper_Complete.md"
SECTIONS_DIR="sections"

# Create build directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Clear output file if it exists
> "$OUTPUT_FILE"

# Function to add section with spacing
add_section() {
    local file=$1
    if [ -f "$SECTIONS_DIR/$file" ]; then
        cat "$SECTIONS_DIR/$file" >> "$OUTPUT_FILE"
        echo "" >> "$OUTPUT_FILE"  # Add blank line between sections
        echo "" >> "$OUTPUT_FILE"  # Add another for spacing
        echo "  ✓ Added $file"
    else
        echo "  ⚠ Warning: $file not found"
    fi
}

echo "Combining sections..."

# Add title page
cat > "$OUTPUT_FILE" << 'TITLE'
# Synchronism: A Comprehensive Model of Reality

**Unified Model of Reality Through Intent Dynamics**

Version: 0.24.09.28.11.00

---

TITLE

# Add sections in proper order
add_section "00-executive-summary.md"
add_section "01-introduction.md"

# Add remaining sections in order, skipping the old introduction
for file in $(ls "$SECTIONS_DIR" | grep -E '^[0-9]{2}-' | sort); do
    if [[ "$file" != "00-executive-summary.md" && "$file" != "01-introduction.md" && "$file" != "00-introduction.md" ]]; then
        add_section "$file"
    fi
done

# Add timestamp
echo "" >> "$OUTPUT_FILE"
echo "---" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "*Generated: $(date '+%Y-%m-%d %H:%M:%S')*" >> "$OUTPUT_FILE"
echo "*Source: Synchronism_0.pdf*" >> "$OUTPUT_FILE"
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
    
    # Copy to docs for GitHub Pages access
    DOCS_DIR="../docs/whitepaper"
    if [ ! -d "$DOCS_DIR" ]; then
        mkdir -p "$DOCS_DIR"
        echo "📁 Created docs/whitepaper directory"
    fi
    
    cp "$OUTPUT_FILE" "$DOCS_DIR/"
    echo "📄 Copied markdown to GitHub Pages location: $DOCS_DIR/Synchronism_Whitepaper_Complete.md"
fi