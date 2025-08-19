#!/bin/bash

# make-md.sh - Combine Synchronism whitepaper sections into monolithic markdown
# Usage: ./make-md.sh

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
    
    # Process all .md files in directory (except index.md and meta directories)
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
    
    # Process subdirectories in sorted order (excluding meta and archive)
    for subdir in $(ls -d "$dir"/*/ 2>/dev/null | sort); do
        if [ -d "$subdir" ]; then
            local dirname=$(basename "$subdir")
            # Skip meta and archive directories
            if [ "$dirname" != "meta" ] && [ "$dirname" != "archive" ]; then
                process_section "$subdir" $((depth + 1))
            fi
        fi
    done
}

# Function to collect all proposals for Appendix B
collect_proposals() {
    local proposals_found=false
    
    echo "" >> "$OUTPUT_FILE"
    echo "---" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    echo "# Appendix B: Current Proposals" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    echo "*This appendix contains all active proposals for improvements to the Synchronism whitepaper. These are suggestions under review and not yet integrated into the main text.*" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    
    # Find all proposal files across all sections
    for section_dir in "$SECTIONS_DIR"/*; do
        if [ -d "$section_dir" ]; then
            local section_name=$(basename "$section_dir")
            local proposals_dir="$section_dir/meta/proposals"
            
            if [ -d "$proposals_dir" ] && [ "$(ls -A "$proposals_dir"/*.md 2>/dev/null)" ]; then
                if [ "$proposals_found" = false ]; then
                    proposals_found=true
                fi
                
                echo "" >> "$OUTPUT_FILE"
                echo "## Proposals for Section: $section_name" >> "$OUTPUT_FILE"
                echo "" >> "$OUTPUT_FILE"
                
                # Process each proposal file
                for proposal in "$proposals_dir"/*.md; do
                    if [ -f "$proposal" ]; then
                        local proposal_name=$(basename "$proposal" .md)
                        echo "### $proposal_name" >> "$OUTPUT_FILE"
                        echo "" >> "$OUTPUT_FILE"
                        cat "$proposal" >> "$OUTPUT_FILE"
                        echo "" >> "$OUTPUT_FILE"
                        echo "---" >> "$OUTPUT_FILE"
                        echo "" >> "$OUTPUT_FILE"
                        
                        echo "  ✓ Added proposal: $section_name/$proposal_name"
                    fi
                done
            fi
        fi
    done
    
    if [ "$proposals_found" = false ]; then
        echo "*No active proposals at this time.*" >> "$OUTPUT_FILE"
        echo "" >> "$OUTPUT_FILE"
    fi
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

# Add Appendix B with all proposals
echo "Collecting proposals for Appendix B..."
collect_proposals

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