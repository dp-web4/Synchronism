#!/bin/bash

# make-web.sh - Generate interactive web version of Synchronism whitepaper
# Usage: ./make-web.sh

echo "Building Synchronism whitepaper web version..."

OUTPUT_DIR="build/web"
SECTIONS_DIR="sections"
ASSETS_DIR="$OUTPUT_DIR/assets"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ASSETS_DIR"

echo "  ✓ Created output directories"

# Create main CSS file
cat > "$ASSETS_DIR/style.css" << 'CSS'
/* Synchronism Whitepaper Styles */
:root {
    --primary-color: #374151;
    --secondary-color: #64748b;
    --background: #ffffff;
    --text-color: #1e293b;
    --border-color: #e2e8f0;
    --code-bg: #f8fafc;
    --sidebar-width: 280px;
}

* {
    margin: 0;
    padding: 0;
    box-sizing: border-box;
}

body {
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    line-height: 1.6;
    color: var(--text-color);
    background: var(--background);
}

/* Layout */
.container {
    display: flex;
    min-height: 100vh;
}

/* Sidebar Navigation */
.sidebar {
    width: var(--sidebar-width);
    background: #f8fafc;
    border-right: 1px solid var(--border-color);
    padding: 2rem 1rem;
    position: fixed;
    height: 100vh;
    overflow-y: auto;
}

.sidebar h2 {
    font-size: 1.2rem;
    margin-bottom: 1rem;
    color: var(--primary-color);
}

.nav-list {
    list-style: none;
}

.nav-list li {
    margin-bottom: 0.5rem;
}

.nav-list a {
    color: var(--text-color);
    text-decoration: none;
    padding: 0.5rem 1rem;
    display: block;
    border-radius: 4px;
    transition: all 0.2s;
}

.nav-list a:hover {
    background: white;
    color: var(--primary-color);
}

.nav-list a.active {
    background: var(--primary-color);
    color: white;
}

/* Main Content */
.main-content {
    margin-left: var(--sidebar-width);
    flex: 1;
    padding: 3rem;
    max-width: 900px;
}

/* Section Pages */
.section {
    display: none;
    animation: fadeIn 0.3s;
}

.section.active {
    display: block;
}

@keyframes fadeIn {
    from { opacity: 0; transform: translateY(10px); }
    to { opacity: 1; transform: translateY(0); }
}

/* Typography */
h1 {
    font-size: 2.5rem;
    margin-bottom: 1rem;
    color: var(--primary-color);
    border-bottom: 2px solid var(--border-color);
    padding-bottom: 0.5rem;
}

h2 {
    font-size: 2rem;
    margin: 2rem 0 1rem;
    color: var(--text-color);
}

h3 {
    font-size: 1.5rem;
    margin: 1.5rem 0 0.75rem;
    color: var(--text-color);
}

p {
    margin-bottom: 1rem;
    line-height: 1.8;
}

/* Code */
code {
    background: var(--code-bg);
    padding: 0.2rem 0.4rem;
    border-radius: 3px;
    font-family: 'Courier New', monospace;
    font-size: 0.9em;
}

pre {
    background: var(--code-bg);
    padding: 1rem;
    border-radius: 6px;
    overflow-x: auto;
    margin: 1rem 0;
    border: 1px solid var(--border-color);
}

/* Blockquotes */
blockquote {
    border-left: 4px solid var(--primary-color);
    padding-left: 1rem;
    margin: 1rem 0;
    color: var(--secondary-color);
    font-style: italic;
}

/* Links */
a {
    color: var(--primary-color);
    text-decoration: none;
}

a:hover {
    text-decoration: underline;
}

/* Mobile Responsive */
@media (max-width: 768px) {
    .sidebar {
        width: 100%;
        height: auto;
        position: relative;
        border-right: none;
        border-bottom: 1px solid var(--border-color);
    }
    
    .main-content {
        margin-left: 0;
        padding: 2rem 1rem;
    }
    
    .container {
        flex-direction: column;
    }
}

/* Dark Mode Support */
@media (prefers-color-scheme: dark) {
    :root {
        --background: #0f172a;
        --text-color: #e2e8f0;
        --border-color: #334155;
        --code-bg: #1e293b;
    }
    
    .sidebar {
        background: #1e293b;
    }
    
    .nav-list a:hover {
        background: #334155;
    }
}
CSS

echo "  ✓ Created CSS styles"

# Create JavaScript for navigation
cat > "$ASSETS_DIR/navigation.js" << 'JAVASCRIPT'
// Synchronism Whitepaper Navigation
document.addEventListener('DOMContentLoaded', function() {
    // Get all navigation links and sections
    const navLinks = document.querySelectorAll('.nav-link');
    const sections = document.querySelectorAll('.section');
    
    // Function to show a specific section
    function showSection(sectionId) {
        // Hide all sections
        sections.forEach(section => {
            section.classList.remove('active');
        });
        
        // Remove active class from all nav links
        navLinks.forEach(link => {
            link.classList.remove('active');
        });
        
        // Show the selected section
        const targetSection = document.getElementById(sectionId);
        if (targetSection) {
            targetSection.classList.add('active');
            
            // Mark the corresponding nav link as active
            const targetLink = document.querySelector(`[data-section="${sectionId}"]`);
            if (targetLink) {
                targetLink.classList.add('active');
            }
            
            // Update URL hash
            window.location.hash = sectionId;
            
            // Scroll to top of content
            window.scrollTo(0, 0);
        }
    }
    
    // Add click handlers to navigation links
    navLinks.forEach(link => {
        link.addEventListener('click', function(e) {
            e.preventDefault();
            const sectionId = this.getAttribute('data-section');
            showSection(sectionId);
        });
    });
    
    // Handle direct URL access with hash
    function handleHashChange() {
        const hash = window.location.hash.slice(1);
        if (hash) {
            showSection(hash);
        } else {
            // Show first section by default (introduction)
            showSection('introduction');
        }
    }
    
    // Handle browser back/forward buttons
    window.addEventListener('hashchange', handleHashChange);
    
    // Initialize on page load
    handleHashChange();
});
JAVASCRIPT

echo "  ✓ Created JavaScript navigation"

# Function to convert markdown to basic HTML
md_to_html() {
    local input=$1
    local output=$2
    
    # Use pandoc if available, otherwise basic conversion
    if command -v pandoc &> /dev/null; then
        pandoc "$input" -f markdown -t html -o "$output" 2>/dev/null
    else
        # Basic markdown to HTML conversion
        sed -e 's/^# \(.*\)/<h1>\1<\/h1>/' \
            -e 's/^## \(.*\)/<h2>\1<\/h2>/' \
            -e 's/^### \(.*\)/<h3>\1<\/h3>/' \
            -e 's/\*\*\([^*]*\)\*\*/<strong>\1<\/strong>/g' \
            -e 's/\*\([^*]*\)\*/<em>\1<\/em>/g' \
            -e 's/^- \(.*\)/<li>\1<\/li>/' \
            -e 's/^[0-9]\+\. \(.*\)/<li>\1<\/li>/' \
            "$input" > "$output"
    fi
}

# Create main HTML file
cat > "$OUTPUT_DIR/index.html" << 'HTML'
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Synchronism: A Comprehensive Model of Reality</title>
    <link rel="stylesheet" href="assets/style.css">
</head>
<body>
    <div class="container">
        <!-- Sidebar Navigation -->
        <nav class="sidebar">
            <h2>Synchronism Whitepaper</h2>
            
            <!-- Navigation List -->
            <ul class="nav-list">
                <li><a href="#introduction" class="nav-link active" data-section="introduction">Introduction</a></li>
                <li><a href="#origins" class="nav-link" data-section="origins">Origins</a></li>
                <li><a href="#hermetic" class="nav-link" data-section="hermetic">Hermetic Principles</a></li>
                <li><a href="#concepts" class="nav-link" data-section="concepts">Fundamental Concepts</a></li>
                <li><a href="#quantum" class="nav-link" data-section="quantum">Quantum Perspective</a></li>
                <li><a href="#implications" class="nav-link" data-section="implications">Implications</a></li>
                <li><a href="#mathematical" class="nav-link" data-section="mathematical">Mathematical Framework</a></li>
                <li><a href="#conclusion" class="nav-link" data-section="conclusion">Conclusion</a></li>
                <li><a href="#appendices" class="nav-link" data-section="appendices">Appendices</a></li>
            </ul>
            
            <!-- Download Links -->
            <div style="margin-top: 2rem; padding-top: 1rem; border-top: 1px solid var(--border-color);">
                <h3 style="font-size: 1rem; margin-bottom: 0.5rem;">Downloads</h3>
                <a href="../Synchronism_Whitepaper.pdf" style="display: block; margin-bottom: 0.5rem;">📕 PDF Version</a>
                <a href="../Synchronism_Whitepaper_Complete.md" style="display: block;">📄 Markdown Version</a>
            </div>
        </nav>
        
        <!-- Main Content Area -->
        <main class="main-content">
            <header style="margin-bottom: 3rem;">
                <h1>Synchronism</h1>
                <p style="font-size: 1.2rem; color: var(--secondary-color);">
                    A Comprehensive Model of Reality Through Intent Dynamics
                </p>
            </header>
            
            <!-- Content Sections (dynamically populated) -->
HTML

# Process sections and add to HTML
echo "Converting sections to HTML..."

# Get list of section files
section_files=($(ls "$SECTIONS_DIR" | grep -E '^[0-9]{2}-' | sort))

# Map sections to IDs
declare -A section_ids=(
    ["00-introduction.md"]="introduction"
    ["02-section-2.md"]="origins"
    ["04-section-4.md"]="hermetic"
    ["05-section-5.md"]="concepts"
    ["07-quantum-cosmic-bridge.md"]="quantum"
    ["14-section-14.md"]="implications"
    ["08-mathematical-framework.md"]="mathematical"
    ["11-conclusion.md"]="conclusion"
    ["13-appendices.md"]="appendices"
)

# Process each section
for file in "${section_files[@]}"; do
    section_id="${section_ids[$file]:-section-${file%.md}}"
    
    if [ -f "$SECTIONS_DIR/$file" ]; then
        # Add appropriate class for first section
        if [ "$section_id" = "introduction" ]; then
            echo "            <section id=\"$section_id\" class=\"section active\">" >> "$OUTPUT_DIR/index.html"
        else
            echo "            <section id=\"$section_id\" class=\"section\">" >> "$OUTPUT_DIR/index.html"
        fi
        
        # Convert markdown to HTML and append
        md_to_html "$SECTIONS_DIR/$file" "$OUTPUT_DIR/temp_section.html"
        cat "$OUTPUT_DIR/temp_section.html" >> "$OUTPUT_DIR/index.html"
        rm "$OUTPUT_DIR/temp_section.html"
        
        echo "            </section>" >> "$OUTPUT_DIR/index.html"
        echo "  ✓ Converted $file"
    fi
done

# Close HTML
cat >> "$OUTPUT_DIR/index.html" << 'HTML'
        </main>
    </div>
    
    <script src="assets/navigation.js"></script>
</body>
</html>
HTML

echo "  ✓ Created main HTML file"

# Create README for web directory
cat > "$OUTPUT_DIR/README.md" << 'README'
# Synchronism Whitepaper - Web Version

This is the interactive web version of the Synchronism whitepaper.

## Features
- Interactive navigation
- Responsive design
- Dark mode support
- Print-friendly CSS

## Viewing
Open `index.html` in any modern web browser.

## Source
Generated from: Synchronism_0.pdf
Latest version: https://dpcars.net/Synchronism_0.pdf
README

echo "  ✓ Created README for web directory"

echo ""
echo "✅ Web version created successfully!"
echo ""
echo "📊 Output Structure:"
echo "   build/web/"
echo "   ├── index.html (main file)"
echo "   ├── assets/"
echo "   │   ├── style.css"
echo "   │   └── navigation.js"
echo "   └── README.md"
echo ""

# Copy to docs for GitHub Pages
echo "📋 Copying to GitHub Pages location..."
DOCS_DIR="../docs/whitepaper"
if [ ! -d "$DOCS_DIR" ]; then
    mkdir -p "$DOCS_DIR"
    echo "📁 Created docs/whitepaper directory"
fi

# Copy all web files
cp -r "$OUTPUT_DIR/"* "$DOCS_DIR/"
echo "🌐 Copied web files to: $DOCS_DIR/"

echo ""
echo "✅ GitHub Pages deployment ready at: https://dp-web4.github.io/Synchronism/whitepaper/"