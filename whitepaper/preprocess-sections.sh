#!/bin/bash

# Preprocessor for section files - demotes headers to prevent TOC expansion
# Keeps ## for section headers, demotes ### to ####, #### to #####, etc.

echo "Preprocessing section files - demoting subsection headers..."

# Process all section markdown files
find sections -type f -name "*.md" \
    ! -path "*/meta/*" \
    ! -path "*/archive/*" \
    ! -path "*/index.md" \
    -print0 | while IFS= read -r -d '' section_file; do
    
    # Skip proposal files (they have their own processor)
    if [[ "$section_file" == *"/proposals/"* ]]; then
        continue
    fi
    
    # Create temp file
    temp_file="${section_file}.tmp"
    
    # Process the file - convert ### headers to bold text to match other sections
    # - Keep ## headers (main section headers)
    # - Convert ### headers to **bold text**
    # - Keep all other lines as-is
    awk '
        # Match headers but preserve the top-level ## headers
        /^##[^#]/ { 
            # This is a section header (##), keep it as is
            print
            next
        }
        /^### / { 
            # Convert ### header to bold text
            # Extract the header text after "### "
            header_text = substr($0, 5)
            print "**" header_text "**"
            next
        }
        # All other lines pass through unchanged
        { print }
    ' "$section_file" > "$temp_file"
    
    # Check if file was actually modified
    if ! cmp -s "$section_file" "$temp_file"; then
        mv "$temp_file" "$section_file"
        echo "  Processed: $section_file"
    else
        rm "$temp_file"
    fi
done

echo "✓ Section preprocessing complete"