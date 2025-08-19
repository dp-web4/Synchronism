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
