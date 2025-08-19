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

# Function to process markdown files recursively
process_section() {
    local dir=$1
    local depth=$2
    
    # Process all .md files in directory (except index.md)
    for file in "$dir"/*.md; do
        if [ -f "$file" ] && [ "$(basename "$file")" != "index.md" ]; then
            # Add section separator for clarity
            if [ $depth -eq 1 ]; then
                echo "" >> "$OUTPUT_FILE"
                echo "---" >> "$OUTPUT_FILE"
                echo "" >> "$OUTPUT_FILE"
            fi
            
            cat "$file" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            
            # Report progress
            local rel_path=${file#$SECTIONS_DIR/}
            echo "  ✓ Added $rel_path"
        fi
    done
    
    # Process subdirectories in sorted order
    for subdir in $(ls -d "$dir"/*/ 2>/dev/null | sort); do
        if [ -d "$subdir" ]; then
            process_section "$subdir" $((depth + 1))
        fi
    done
}

echo "Combining sections..."

# Add title page
cat > "$OUTPUT_FILE" << 'TITLE'
# Synchronism: A Comprehensive Model of Reality

**Unified Model of Reality Through Intent Dynamics**

TITLE

# Process sections in specific order
sections=(
    "00-executive-summary"
    "01-introduction"
    "02-perspective"
    "03-hermetic-principles"
    "04-fundamental-concepts"
    "05-quantum-macro"
    "06-implications"
    "07-conclusion"
    "08-glossary"
    "09-appendix-mathematical"
)

for section in "${sections[@]}"; do
    if [ -d "$SECTIONS_DIR/$section" ]; then
        echo "Processing $section..."
        process_section "$SECTIONS_DIR/$section" 1
    else
        echo "  ⚠ Warning: $section directory not found"
    fi
done

# Document complete - no footer needed

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