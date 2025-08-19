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
