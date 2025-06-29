# Testing the Synchronism Web Interface

## Issue with File:// Protocol

The HTML-based navigation system uses `fetch()` to load section files dynamically. This requires a web server because browsers block local file access via fetch for security reasons.

## Solution: Local Web Server

### Quick Start (Python)
```bash
# From the web-version directory:
python start-server.py
```

### Alternative (Python built-in)
```bash
# From the web-version directory:
python -m http.server 8000
```

### Alternative (Node.js)
```bash
# Install globally: npm install -g http-server
http-server -p 8000
```

## Testing the Navigation

1. Start the local server (choose one method above)
2. Open http://localhost:8000 in your browser
3. Click on navigation links to test section loading
4. Check browser console (F12) for any errors

## Expected Behavior

- ✅ Introduction loads by default
- ✅ Navigation links load corresponding HTML sections
- ✅ Status indicator shows current section
- ✅ Cross-references work between sections
- ✅ MathJax renders mathematical content

## Current Implementation Status

### Completed Sections (with HTML files):
- ✅ 1. Introduction
- ✅ 2. Importance of Perspective  
- ✅ 3. Hermetic Principles
- ✅ 4. Fundamental Concepts (Chapter overview)
- ✅ 4.1-4.12 All fundamental concept subsections

### Pending Sections:
- ⏳ Chapter 5: Quantum & Macro Phenomena (22 sections)
- ⏳ Chapter 6: Implications & Applications (4 major sections)
- ⏳ Chapter 7: Conclusion
- ⏳ Glossary
- ⏳ Appendix A: Mathematics

## Troubleshooting

### Navigation Not Loading
- Ensure you're accessing via http://localhost:8000, not file://
- Check browser console for fetch errors
- Verify section HTML files exist in correct directories

### Missing Sections
- Sections without HTML files will show "Section Not Yet Available"
- This is expected behavior for sections not yet implemented

### Styling Issues
- Ensure styles.css is in the web-version directory
- Check browser console for CSS loading errors