#!/bin/bash

# cleanup_formatting.sh - Remove navigation hints and convert ### to bold

echo "Cleaning up whitepaper formatting..."
echo "=================================="

# Counter for changes
nav_removed=0
headings_converted=0

# Process all markdown files
for file in sections/**/*.md sections/**/**/*.md; do
    if [ -f "$file" ]; then
        # Check if file has changes to make
        has_changes=false
        
        # Check for Continue Reading sections
        if grep -q "#### Continue Reading" "$file"; then
            # Remove the entire Continue Reading section (from #### to the next empty line or end)
            sed -i '/^---$/,/^---$/d' "$file" 2>/dev/null || sed -i '' '/^---$/,/^---$/d' "$file"
            sed -i '/^#### Continue Reading$/,/^$/d' "$file" 2>/dev/null || sed -i '' '/^#### Continue Reading$/,/^$/d' "$file"
            echo "✓ Removed navigation hints from: $file"
            ((nav_removed++))
            has_changes=true
        fi
        
        # Check for ### headings (these should be bold instead)
        if grep -q "^### " "$file"; then
            # Store the file content
            temp_file="${file}.tmp"
            
            # Convert ### headings to bold
            while IFS= read -r line; do
                if [[ $line =~ ^###\ (.+)$ ]]; then
                    echo "**${BASH_REMATCH[1]}**"
                    echo ""
                else
                    echo "$line"
                fi
            done < "$file" > "$temp_file"
            
            mv "$temp_file" "$file"
            echo "✓ Converted ### to bold in: $file"
            ((headings_converted++))
            has_changes=true
        fi
    fi
done

echo ""
echo "Cleanup Summary:"
echo "================"
echo "Navigation sections removed: $nav_removed"
echo "Headings converted to bold: $headings_converted files"

echo ""
echo "Now running build to update all outputs..."
./build.sh

echo ""
echo "✅ Formatting cleanup complete!"