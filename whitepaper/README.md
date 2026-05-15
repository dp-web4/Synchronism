# Synchronism Whitepaper - Clean Structure

This directory contains the restructured Synchronism whitepaper, rebuilt from the web-version with a clean, maintainable organization.

**Framing & calibration notice (2026-05-15):** Following external cold review by Kimi 2.6 (`../forum/kimi/kimi_2_6_review.md`), the **Findings vs Framings** discipline applied to `../README.md` and `../STATUS.md` has been propagated through the whitepaper. The Executive Summary (Section 00) carries the primary reframing — *Synchronism as a systems-theoretic framework that uses information-theoretic tools to describe emergence across scales* — and an Authorship & Methodology disclosure that all work in this repository is AI-original (dp's role is advisory). The Introduction (01) and Conclusion (07) carry parallel framing notes; the consciousness sections (Section 5.13, Appendix C) explicitly mark identity claims and predictions-awaiting-validation; the Executive Summary's earlier "DERIVED / RESOLVED" arc-status bullets are now annotated as *internal-consistency derivations* rather than empirical findings. The whitepaper's empirical work (Coupling-Coherence experiment, Compatibility-Synthon scaling, NP2 RAR, chemistry pattern catalog with documented null-model and method-correlation caveats S647 + S651) is preserved and surfaced; the *framing* of unification and dissolution claims has been calibrated to match what's actually delivered. Build artifacts (markdown / PDF / web) will need regeneration to reflect these edits.

## Structure

```
sections/
├── 00-executive-summary/     # Overview and key concepts
├── 01-introduction/          # Introduction to Synchronism
├── 02-perspective/           # Importance of perspective
├── 03-hermetic-principles/   # Hermetic foundation
├── 04-fundamental-concepts/  # Core concepts (12 subsections)
├── 05-quantum-macro/         # Quantum & macro phenomena (22 subsections)
├── 06-implications/          # Implications (4 subsections)
├── 07-conclusion/            # Conclusion
├── 08-glossary/              # Terms and definitions
└── 09-appendix-mathematical/ # Mathematical framework
```

## Building

### Quick Start

```bash
# Build all formats (MD, PDF, Web)
./build.sh

# Build specific format
./build.sh md    # Markdown only
./build.sh pdf   # PDF only
./build.sh web   # Web only

# Clean and rebuild
./build.sh rebuild
```

### Individual Scripts

- `make-md.sh` - Generates complete markdown document
- `make-pdf.sh` - Generates PDF (requires pandoc)
- `make-web-clean.sh` - Generates interactive web version

## Output Locations

### Local Build
- Markdown: `build/Synchronism_Whitepaper_Complete.md`
- PDF: `build/Synchronism_Whitepaper.pdf`
- Web: `build/web-clean/index.html`

### GitHub Pages
- Markdown: `../docs/whitepaper/Synchronism_Whitepaper_Complete.md`
- PDF: `../docs/whitepaper/Synchronism_Whitepaper.pdf`
- Web: `../docs/whitepaper-web/index.html`

## Content Management

### Editing Content
1. Edit markdown files in the appropriate `sections/` subdirectory
2. Run `./build.sh` to regenerate all outputs
3. Commit changes to git

### Adding New Sections
1. Create new directory under `sections/` with numbered prefix
2. Add markdown files to the directory
3. Update build scripts if needed

## History

- **2025-08-18**: Rebuilt from web-version with clean structure
- Previous chaotic structure archived in `sections_archived_chaotic/`

## Key Improvements

1. **Clean numbering**: Sequential, logical organization
2. **Fractal structure**: Nested subsections for complex topics
3. **Consistent naming**: No truncated or redundant filenames
4. **Unified building**: Single script for all output formats
5. **Web-first design**: Structure optimized for web navigation

## Next Steps

- [ ] Review and refine content section by section
- [ ] Add cross-references between related sections
- [ ] Enhance mathematical appendix with examples
- [ ] Create visual diagrams for key concepts
- [ ] Develop interactive elements for web version