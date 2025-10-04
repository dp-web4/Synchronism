#!/bin/bash

# Preprocess proposal files to demote headers by 3 levels
# This ensures consistent formatting regardless of how proposals are submitted

SECTIONS_DIR="sections"

echo "Preprocessing proposal files - demoting headers..."

# Find all proposal markdown files
find "$SECTIONS_DIR" -path "*/meta/proposals/*.md" -type f | while IFS= read -r proposal_file; do
    echo "  Processing: ${proposal_file#$SECTIONS_DIR/}"
    
    # Create temporary file
    temp_file="${proposal_file}.tmp"
    
    # Demote headers by 3 levels: # → ####, ## → #####, ### → ######
    # But for section headers (##), make them ###### to match what scripts expect
    # Use awk to avoid sequential replacement issues
    awk '
        /^###[^#]/ { print "######" substr($0, 4); next }
        /^##[^#]/ { print "######" substr($0, 3); next }
        /^#[^#]/ { print "####" substr($0, 2); next }
        { print }
    ' "$proposal_file" > "$temp_file"
    
    # Replace original file
    mv "$temp_file" "$proposal_file"
done

echo "✓ Proposal preprocessing complete"