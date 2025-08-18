#!/bin/bash

# fix_all_formatting.sh - Remove all navigation sections and convert ALL ### to bold

echo "Comprehensive formatting cleanup..."
echo "=================================="

cd /mnt/c/projects/ai-agents/Synchronism/whitepaper

# Counter for changes
nav_removed=0
headings_converted=0

# First, remove all Continue Exploring sections
echo "Removing Continue Exploring sections..."
for file in sections/**/*.md sections/**/**/*.md; do
    if [ -f "$file" ]; then
        if grep -q "Continue Exploring" "$file"; then
            # Create temp file
            temp_file="${file}.tmp"
            
            # Process file line by line, skip Continue Exploring sections
            in_continue_section=false
            while IFS= read -r line; do
                # Check if we're entering a Continue Exploring section
                if [[ $line == *"#### Continue Exploring"* ]]; then
                    in_continue_section=true
                    ((nav_removed++))
                    continue
                fi
                
                # Check if we're in the section (skip until we hit non-link content)
                if [ "$in_continue_section" = true ]; then
                    # Skip empty lines and lines with links
                    if [[ -z "$line" || $line == "- ["* || $line == "---"* ]]; then
                        continue
                    else
                        # We've hit content that's not part of the nav section
                        in_continue_section=false
                        # Don't skip this line, it's real content
                    fi
                fi
                
                # Write non-navigation lines
                if [ "$in_continue_section" = false ]; then
                    echo "$line"
                fi
            done < "$file" > "$temp_file"
            
            mv "$temp_file" "$file"
            echo "  ✓ Removed navigation from: $file"
        fi
    fi
done

echo ""
echo "Converting ALL ### headings to bold text..."

# Now convert ALL ### headings to bold
for file in sections/**/*.md sections/**/**/*.md; do
    if [ -f "$file" ]; then
        if grep -q "### " "$file"; then
            # Create temp file
            temp_file="${file}.tmp"
            
            # Convert ### headings to bold (handle leading spaces)
            sed 's/^[[:space:]]*### \(.*\)/**\1**/' "$file" > "$temp_file"
            
            mv "$temp_file" "$file"
            echo "  ✓ Converted ### to bold in: $file"
            ((headings_converted++))
        fi
    fi
done

echo ""
echo "Cleanup Summary:"
echo "================"
echo "Continue Exploring sections removed: $nav_removed"
echo "Files with ### converted to bold: $headings_converted"

echo ""
echo "Rebuilding all outputs..."
./build.sh

echo ""
echo "✅ Comprehensive formatting cleanup complete!"