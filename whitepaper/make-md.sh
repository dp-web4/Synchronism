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

# Preprocess proposals and sections to demote headers
./preprocess-proposals.sh
./preprocess-sections.sh

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

# Function to collect all proposals for Appendix D
collect_proposals() {
    local proposals_found=false
    
    echo "" >> "$OUTPUT_FILE"
    echo "---" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    echo "# Appendix D: Current Proposals" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    echo "*This appendix contains all active proposals for improvements to the Synchronism whitepaper. These are suggestions under review and not yet integrated into the main text.*" >> "$OUTPUT_FILE"
    echo "" >> "$OUTPUT_FILE"
    
    # Find ALL proposal files recursively
    local proposal_files=$(find "$SECTIONS_DIR" -path "*/meta/proposals/*.md" -type f 2>/dev/null | sort)
    
    if [ -n "$proposal_files" ]; then
        proposals_found=true
        local current_section=""
        
        # Skip navigation for markdown/PDF - go straight to proposals
        while IFS= read -r proposal_file; do
            # Extract section path relative to SECTIONS_DIR
            local rel_path="${proposal_file#$SECTIONS_DIR/}"
            local section_path=$(echo "$rel_path" | sed 's|/meta/proposals/.*||')
            
            # Print section header if it's a new section (using #### to stay out of TOC)
            if [ "$section_path" != "$current_section" ]; then
                current_section="$section_path"
                echo "" >> "$OUTPUT_FILE"
                echo "#### Proposals for: $section_path" >> "$OUTPUT_FILE"
                echo "" >> "$OUTPUT_FILE"
                echo "" >> "$OUTPUT_FILE"  # Extra line for better PDF separation
            fi
            
            # Extract key info from the proposal file
            local proposal_id=$(grep "\*\*ID\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
            local title=$(grep -m1 "^#### Proposal" "$proposal_file" | sed 's/^#### Proposal [0-9]*: //')
            # Extract author - handle both formats: "Author" and "Author (LCT: hash)"
            local author=$(grep "\*\*Author\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //' | sed 's/ (LCT:.*)$//')
            local date=$(grep "\*\*Date\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
            local status=$(grep "\*\*Status\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
            local type=$(grep "\*\*Type\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
            
            # Single line with all metadata
            echo "**${proposal_id}. ${title}** — ${author} | ${date} | ${status} | ${type}" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            
            # Extract just the content sections, skip metadata and current state
            awk '
                /^###### Proposed Change/,/^###### Rationale/ { 
                    if (!/^######/) print 
                }
                /^###### Specific Text Changes/,/^###### Impact Assessment/ { 
                    if (!/^###### Impact Assessment/) {
                        if (/^######/) {
                            # Convert heading to bold text
                            gsub(/^###### /, "**", $0)
                            gsub(/:?$/, ":**", $0)
                        }
                        print
                    }
                }
            ' "$proposal_file" >> "$OUTPUT_FILE"
            
            echo "" >> "$OUTPUT_FILE"
            echo "---" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            
            echo "  ✓ Added proposal: $section_path/$proposal_name"
        done <<< "$proposal_files"
    fi
    
    if [ "$proposals_found" = false ]; then
        echo "*No active proposals at this time.*" >> "$OUTPUT_FILE"
        echo "" >> "$OUTPUT_FILE"
    fi
}

echo "Combining sections..."

# Add title page
cat > "$OUTPUT_FILE" << 'TITLE'
# Synchronism: A Computational Framework for Pattern Dynamics

**Non-Anthropocentric Model of Reality (Working Draft)**

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

# Add Appendix D with all proposals
echo "Collecting proposals for Appendix D..."
collect_proposals

# Add Appendix E with changelog
echo "Creating Appendix E for changelog..."
echo "" >> "$OUTPUT_FILE"
echo "---" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "# Appendix E: Change Log" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"
echo "*Version history and evolution of the Synchronism whitepaper.*" >> "$OUTPUT_FILE"
echo "" >> "$OUTPUT_FILE"

# Collect all changelog entries
for section in "$SECTIONS_DIR"/*; do
    if [ -d "$section" ]; then
        changelog_file="$section/meta/CHANGELOG.md"
        if [ -f "$changelog_file" ]; then
            section_name=$(basename "$section")
            echo "#### $section_name" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            # Skip format section and downgrade headers as safety measure
            # Extract only the Entries section onwards (skip Format section)
            awk '/^###### Entries/,EOF' "$changelog_file" | \
            # Downgrade headers (h1->h4, h2->h5, h3->h6)
            sed 's/^###/######/g; s/^##/#####/g; s/^#\([^#]\)/####\1/g' >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            echo "---" >> "$OUTPUT_FILE"
            echo "" >> "$OUTPUT_FILE"
            echo "  ✓ Added changelog for $section_name"
        fi
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