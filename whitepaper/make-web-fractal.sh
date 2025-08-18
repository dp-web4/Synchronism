#!/bin/bash

# make-web-fractal.sh - Generate interactive web version from fractal structure
# Usage: ./make-web-fractal.sh

echo "Building Synchronism whitepaper web version from fractal structure..."

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
    --sidebar-width: 320px;
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
    margin-bottom: 0.25rem;
}

.nav-list a {
    color: var(--text-color);
    text-decoration: none;
    padding: 0.4rem 0.8rem;
    display: block;
    border-radius: 4px;
    transition: all 0.2s;
    font-size: 0.95rem;
}

.nav-list a:hover {
    background: white;
    color: var(--primary-color);
}

.nav-list a.active {
    background: var(--primary-color);
    color: white;
}

/* Nested navigation */
.nav-list .nav-list {
    margin-left: 1rem;
    margin-top: 0.25rem;
    font-size: 0.9rem;
}

.nav-section {
    margin-bottom: 1rem;
}

.nav-section-title {
    font-weight: 600;
    color: var(--primary-color);
    padding: 0.5rem 0;
    border-bottom: 1px solid var(--border-color);
    margin-bottom: 0.5rem;
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

h4 {
    font-size: 1.25rem;
    margin: 1rem 0 0.5rem;
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
            // Show executive summary by default
            showSection('executive-summary');
        }
    }
    
    // Handle browser back/forward buttons
    window.addEventListener('hashchange', handleHashChange);
    
    // Initialize on page load
    handleHashChange();
});
JAVASCRIPT

echo "  ✓ Created JavaScript navigation"

# Function to convert markdown to HTML
md_to_html() {
    local input=$1
    local output=$2
    
    if command -v pandoc &> /dev/null; then
        pandoc "$input" -f markdown -t html -o "$output" 2>/dev/null
    else
        # Basic markdown to HTML conversion
        sed -e 's/^# \(.*\)/<h1>\1<\/h1>/' \
            -e 's/^## \(.*\)/<h2>\1<\/h2>/' \
            -e 's/^### \(.*\)/<h3>\1<\/h3>/' \
            "$input" > "$output"
    fi
}

# Start HTML file
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
            
            <!-- Navigation List (will be generated) -->
            <ul class="nav-list">
HTML

# Function to process directory for navigation
process_nav_directory() {
    local dir=$1
    local level=$2
    local parent_id=$3
    
    for item in "$dir"/*; do
        if [ -d "$item" ]; then
            dirname=$(basename "$item")
            # Clean name for display
            clean_name=${dirname#*-}
            clean_name=${clean_name//-/ }
            clean_name=$(echo "$clean_name" | sed 's/\b\(.\)/\u\1/g')
            
            # Create section ID
            section_id="${parent_id}-${dirname}"
            section_id=${section_id#-}
            
            echo "                <li class=\"nav-section\">" >> "$OUTPUT_DIR/index.html"
            echo "                    <div class=\"nav-section-title\">$clean_name</div>" >> "$OUTPUT_DIR/index.html"
            echo "                    <ul class=\"nav-list\">" >> "$OUTPUT_DIR/index.html"
            
            # Process subdirectory
            process_nav_directory "$item" $((level+1)) "$section_id"
            
            echo "                    </ul>" >> "$OUTPUT_DIR/index.html"
            echo "                </li>" >> "$OUTPUT_DIR/index.html"
            
        elif [ -f "$item" ] && [[ "$item" == *.md ]] && [[ "$(basename "$item")" != "index.md" ]]; then
            filename=$(basename "$item")
            # Clean name for display
            clean_name=${filename%.md}
            clean_name=${clean_name#*-}
            clean_name=${clean_name//-/ }
            clean_name=$(echo "$clean_name" | sed 's/\b\(.\)/\u\1/g')
            
            # Create section ID
            section_id="${parent_id}-${filename%.md}"
            section_id=${section_id#-}
            section_id=${section_id//\//-}
            
            echo "                <li><a href=\"#$section_id\" class=\"nav-link\" data-section=\"$section_id\">$clean_name</a></li>" >> "$OUTPUT_DIR/index.html"
        fi
    done
}

# Generate navigation
echo "Generating navigation structure..."
process_nav_directory "$SECTIONS_DIR" 0 ""

# Continue HTML
cat >> "$OUTPUT_DIR/index.html" << 'HTML'
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

# Function to process directory for content sections
process_content_directory() {
    local dir=$1
    local parent_id=$2
    
    for item in "$dir"/*; do
        if [ -d "$item" ]; then
            # Process subdirectory recursively
            dirname=$(basename "$item")
            section_id="${parent_id}-${dirname}"
            section_id=${section_id#-}
            process_content_directory "$item" "$section_id"
            
        elif [ -f "$item" ] && [[ "$item" == *.md ]] && [[ "$(basename "$item")" != "index.md" ]]; then
            filename=$(basename "$item")
            section_id="${parent_id}-${filename%.md}"
            section_id=${section_id#-}
            section_id=${section_id//\//-}
            
            # Determine if this is the first section (executive summary)
            if [[ "$filename" == "00-executive-summary.md" ]]; then
                echo "            <section id=\"$section_id\" class=\"section active\">" >> "$OUTPUT_DIR/index.html"
            else
                echo "            <section id=\"$section_id\" class=\"section\">" >> "$OUTPUT_DIR/index.html"
            fi
            
            # Convert markdown to HTML and append
            md_to_html "$item" "$OUTPUT_DIR/temp_section.html"
            cat "$OUTPUT_DIR/temp_section.html" >> "$OUTPUT_DIR/index.html"
            rm "$OUTPUT_DIR/temp_section.html"
            
            echo "            </section>" >> "$OUTPUT_DIR/index.html"
            echo "  ✓ Processed $(basename "$item")"
        fi
    done
}

# Generate content sections
echo "Converting sections to HTML..."
process_content_directory "$SECTIONS_DIR" ""

# Close HTML
cat >> "$OUTPUT_DIR/index.html" << 'HTML'
        </main>
    </div>
    
    <script src="assets/navigation.js"></script>
</body>
</html>
HTML

echo "  ✓ Created main HTML file"

echo ""
echo "✅ Web version created successfully!"
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

# Also copy PDF and MD if they exist
if [ -f "build/Synchronism_Whitepaper.pdf" ]; then
    cp "build/Synchronism_Whitepaper.pdf" "$DOCS_DIR/"
    echo "📕 Copied PDF to docs"
fi

if [ -f "build/Synchronism_Whitepaper_Complete.md" ]; then
    cp "build/Synchronism_Whitepaper_Complete.md" "$DOCS_DIR/"
    echo "📄 Copied markdown to docs"
fi

echo ""
echo "✅ GitHub Pages deployment ready at: https://dp-web4.github.io/Synchronism/whitepaper/"