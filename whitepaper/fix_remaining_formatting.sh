#!/bin/bash

# Fix remaining formatting issues in whitepaper sections

echo "Fixing remaining formatting issues..."
echo ""

# Counter for changes
navigation_removed=0
headings_converted=0
total_files=0

# Process all markdown files
find sections -name "*.md" -type f | while read -r file; do
    total_files=$((total_files + 1))
    temp_file="${file}.tmp"
    changes_made=false
    
    # First pass: Remove navigation arrows and links
    if grep -q "→" "$file"; then
        # Remove lines with navigation arrows (Next:, Application:, Philosophy:, etc.)
        sed '/^\s*-\s*\[.*→\]/d; /^Next:.*→/d; /^Application:.*→/d; /^Philosophy:.*→/d' "$file" > "$temp_file"
        
        # Check if changes were made
        if ! diff -q "$file" "$temp_file" > /dev/null 2>&1; then
            navigation_removed=$((navigation_removed + 1))
            changes_made=true
            cp "$temp_file" "$file"
            echo "  ✓ Removed navigation arrows from: $file"
        fi
    fi
    
    # Second pass: Convert ### headings (with space) to bold
    if grep -q "^### " "$file"; then
        if [ "$changes_made" = true ]; then
            # Work with the already modified file
            sed 's/^### \(.*\)/**\1**/' "$file" > "$temp_file"
        else
            # Work with original file
            sed 's/^### \(.*\)/**\1**/' "$file" > "$temp_file"
        fi
        
        # Check if changes were made
        if ! diff -q "$file" "$temp_file" > /dev/null 2>&1; then
            headings_converted=$((headings_converted + 1))
            cp "$temp_file" "$file"
            echo "  ✓ Converted ### headings to bold in: $file"
        fi
    fi
    
    # Clean up temp file
    rm -f "$temp_file"
done

echo ""
echo "Summary:"
echo "  - Files with navigation arrows removed: $navigation_removed"
echo "  - Files with ### headings converted: $headings_converted"
echo ""

# Also check for #### headings that should be bold
echo "Checking for #### headings..."
heading_count=$(grep -r "^####" sections --include="*.md" | wc -l)
if [ "$heading_count" -gt 0 ]; then
    echo "Found $heading_count #### headings. Converting to bold..."
    find sections -name "*.md" -type f | while read -r file; do
        if grep -q "^#### " "$file"; then
            temp_file="${file}.tmp"
            sed 's/^#### \(.*\)/**\1**/' "$file" > "$temp_file"
            if ! diff -q "$file" "$temp_file" > /dev/null 2>&1; then
                cp "$temp_file" "$file"
                echo "  ✓ Converted #### headings in: $file"
            fi
            rm -f "$temp_file"
        fi
    done
fi

echo ""
echo "✅ Formatting cleanup complete!"