#!/bin/bash

# Setup governance structure for all whitepaper sections
# Creates meta directories (proposals, reviews, changelog) and archive directories

echo "Setting up governance structure for Synchronism whitepaper sections..."

for dir in */; do
    if [ -d "$dir" ]; then
        section_name=$(basename "$dir")
        echo "Processing: $section_name"
        
        # Create meta subdirectory structure
        mkdir -p "${dir}meta/proposals"
        mkdir -p "${dir}meta/reviews"
        
        # Create archive subdirectory
        mkdir -p "${dir}archive"
        
        # Create section changelog if it doesn't exist
        if [ ! -f "${dir}meta/CHANGELOG.md" ]; then
            cat > "${dir}meta/CHANGELOG.md" << 'EOF'
# Section Changelog

## Format
Each entry should include:
- Date (ISO 8601)
- Author (LCT ID or name)
- Change type (ADD/MODIFY/DELETE/ARCHIVE)
- Description of change
- Rationale for change

## Entries
<!-- Entries added chronologically below -->

EOF
            echo "  ✓ Created meta/CHANGELOG.md"
        fi
        
        # Create section metadata file if it doesn't exist
        if [ ! -f "${dir}meta/metadata.json" ]; then
            cat > "${dir}meta/metadata.json" << EOF
{
  "section_id": "$section_name",
  "section_type": "fractal",
  "maintainers": [],
  "last_modified": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "status": "active",
  "dependencies": [],
  "tags": [],
  "governance": {
    "review_required": true,
    "min_reviewers": 2,
    "arbitration_threshold": 0.7
  },
  "archive_policy": {
    "trigger": "significant_changes",
    "description": "Archive previous version when making substantial removals or structural changes"
  }
}
EOF
            echo "  ✓ Created meta/metadata.json"
        fi
        
        # Create archive README if it doesn't exist
        if [ ! -f "${dir}archive/README.md" ]; then
            cat > "${dir}archive/README.md" << 'EOF'
# Archive Directory

This directory stores previous versions of the section's fractal file (index.md) when significant changes are made.

## Purpose
- Preserve historical versions before major modifications
- Enable rollback if needed
- Maintain intellectual history of the document

## When to Archive
Archive the current index.md before:
- Substantial content removal (>30% of section)
- Major structural reorganization
- Complete rewriting of core concepts
- Integration of conflicting proposals that override existing content

## Naming Convention
Archived files should be named:
`index_YYYY-MM-DD_description.md`

Example: `index_2025-08-19_pre-r6-framework.md`

## Archive Process
1. Copy current index.md to archive with descriptive name
2. Note the archival in meta/CHANGELOG.md
3. Proceed with modifications to index.md
4. Commit both the archive and the new version

## Restoration
To restore an archived version:
1. Review the archived file
2. Copy it back to index.md if appropriate
3. Document the restoration in meta/CHANGELOG.md
EOF
            echo "  ✓ Created archive/README.md"
        fi
        
        echo "  ✓ Structure complete for $section_name"
        echo ""
    fi
done

# Create a global governance documentation file
cat > ../GOVERNANCE_STRUCTURE.md << 'EOF'
# Synchronism Whitepaper Governance Structure

## Directory Organization

Each section of the whitepaper follows a fractal governance structure:

```
section-name/
├── index.md           # The fractal file (actual content)
├── meta/              # Governance metadata
│   ├── proposals/     # Proposed changes from participants
│   ├── reviews/       # Reviews of proposals
│   ├── CHANGELOG.md   # Section change history
│   └── metadata.json  # Section configuration
└── archive/           # Historical versions
    └── README.md      # Archive usage guide
```

## Fractal Files vs Meta Files

### Fractal Files (index.md)
- The actual content of the whitepaper
- Modified only by arbiters after governance approval
- Represents the "truth" of the document
- Triggers rebuilds when changed

### Meta Files
- Governance artifacts (proposals, reviews, changelog)
- Can be modified by participants during their phase
- Track the evolution and decision process
- Do not trigger rebuilds

### Archive Files
- Historical versions of fractal files
- Created at arbiter's discretion before major changes
- Preserve document history and enable rollback
- Named with date and description

## Governance Workflow

1. **Proposal Phase**: Participants create files in `meta/proposals/`
2. **Review Phase**: Reviewers add files in `meta/reviews/`
3. **Arbitration Phase**: 
   - Arbiter may archive current `index.md` if making major changes
   - Arbiter modifies `index.md` based on accepted proposals
   - Changes logged in `meta/CHANGELOG.md`
4. **Build Phase**: System automatically rebuilds if fractal files changed

## Archive Policy

Archive the current index.md when:
- Removing >30% of content
- Restructuring the section significantly
- Replacing core concepts
- Integrating conflicting proposals

Archive naming: `index_YYYY-MM-DD_description.md`

## Participant Roles

- **Proposers**: Create content in `meta/proposals/`
- **Reviewers**: Evaluate proposals in `meta/reviews/`
- **Arbiters**: 
  - Decide what to integrate
  - Archive if needed
  - Modify fractal files
  - Update changelog

## Important Notes

1. **Only arbiters modify fractal files** - This maintains document integrity
2. **Always document in CHANGELOG.md** - Every change needs a record
3. **Archive before major changes** - Preserve history
4. **Meta files are collaborative** - Multiple participants can contribute
5. **Fractal files are authoritative** - They represent consensus

---

*Last Updated: 2025-08-19*
*Governance System: Rev_0*
EOF

echo "✅ Governance structure setup complete!"
echo ""
echo "Created for all sections:"
echo "  - meta/proposals/ (for participant proposals)"
echo "  - meta/reviews/ (for review feedback)"
echo "  - meta/CHANGELOG.md (section history)"
echo "  - meta/metadata.json (section configuration)"
echo "  - archive/ (for historical versions)"
echo "  - archive/README.md (archive usage guide)"
echo ""
echo "📖 See GOVERNANCE_STRUCTURE.md for complete documentation"