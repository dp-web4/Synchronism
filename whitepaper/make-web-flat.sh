#!/bin/bash

# make-web-flat.sh - Generate interactive web version with 2-level navigation
# Each subsection combines all its content into a single scrollable page
# Usage: ./make-web-flat.sh

echo "Building Synchronism whitepaper web version with flattened navigation..."

OUTPUT_DIR="build/web"
SECTIONS_DIR="sections"
ASSETS_DIR="$OUTPUT_DIR/assets"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ASSETS_DIR"

echo "  ✓ Created output directories"

# Create main CSS file
cat > "$ASSETS_DIR/style.css" << 'CSS'
/* Synchronism Whitepaper Styles - Flattened Navigation */
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
    padding: 1rem;
    position: fixed;
    height: 100vh;
    overflow-y: auto;
}

.sidebar h2 {
    font-size: 1.05rem;
    margin-bottom: 0.5rem;
    color: var(--primary-color);
    font-weight: 600;
}

.nav-list {
    list-style: none;
}

.nav-list li {
    margin-bottom: 0;
}

/* Main sections */
.nav-section {
    margin-bottom: 0.5rem;
}

.nav-section-title {
    font-weight: 600;
    color: var(--primary-color);
    padding: 0.3rem 0.6rem;
    border-bottom: 1px solid var(--border-color);
    margin-bottom: 0.2rem;
    font-size: 0.95rem;
    cursor: pointer;
    user-select: none;
    display: flex;
    justify-content: space-between;
    align-items: center;
}

.nav-section-title:hover {
    background: rgba(55, 65, 81, 0.05);
}

.nav-section-title::after {
    content: '▶';
    font-size: 0.7rem;
    transition: transform 0.2s;
}

.nav-section.expanded .nav-section-title::after {
    transform: rotate(90deg);
}

/* Subsection links */
.nav-subsection {
    display: none;
    margin-left: 0.5rem;
}

.nav-section.expanded .nav-subsection {
    display: block;
}

.nav-subsection a {
    color: var(--text-color);
    text-decoration: none;
    padding: 0.2rem 0.6rem;
    display: block;
    border-radius: 4px;
    transition: all 0.2s;
    font-size: 0.85rem;
    line-height: 1.3;
}

.nav-subsection a:hover {
    background: white;
    color: var(--primary-color);
}

.nav-subsection a.active {
    background: var(--primary-color);
    color: white;
}

/* Main Content */
.main-content {
    margin-left: var(--sidebar-width);
    flex: 1;
    padding: 1.5rem 2rem;
    max-width: 1000px;
}

/* Section Pages */
.section {
    display: none;
}

.section.active {
    display: block;
}

/* Typography */
h1 {
    font-size: 2rem;
    margin-bottom: 0.75rem;
    color: var(--primary-color);
    border-bottom: 2px solid var(--border-color);
    padding-bottom: 0.5rem;
}

h2 {
    font-size: 1.6rem;
    margin: 1.5rem 0 0.75rem;
    color: var(--text-color);
    padding-top: 1rem;
}

h3 {
    font-size: 1.3rem;
    margin: 1.25rem 0 0.5rem;
    color: var(--text-color);
}

h4 {
    font-size: 1.15rem;
    margin: 1rem 0 0.5rem;
    color: var(--text-color);
}

p {
    margin-bottom: 0.75rem;
    line-height: 1.65;
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

/* Smooth scrolling */
html {
    scroll-behavior: smooth;
}

/* Section divider */
.section-divider {
    margin: 2rem 0;
    padding: 1rem 0;
    border-top: 2px solid var(--border-color);
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
    
    .nav-subsection a:hover {
        background: #334155;
    }
}
CSS

echo "  ✓ Created CSS styles"

# Create JavaScript for navigation
cat > "$ASSETS_DIR/navigation.js" << 'JAVASCRIPT'
// Synchronism Whitepaper Navigation - Flattened Version
document.addEventListener('DOMContentLoaded', function() {
    const navSections = document.querySelectorAll('.nav-section');
    const navLinks = document.querySelectorAll('.nav-subsection a');
    const sections = document.querySelectorAll('.section');
    const sectionTitles = document.querySelectorAll('.nav-section-title');
    
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
                
                // Expand parent section if needed
                const parentSection = targetLink.closest('.nav-section');
                if (parentSection && !parentSection.classList.contains('expanded')) {
                    // Collapse others first
                    navSections.forEach(section => {
                        if (section !== parentSection) {
                            section.classList.remove('expanded');
                        }
                    });
                    parentSection.classList.add('expanded');
                }
            }
            
            // Update URL hash
            window.location.hash = sectionId;
            
            // Scroll to top of content
            window.scrollTo(0, 0);
        }
    }
    
    // Add click handlers to section titles for expand/collapse
    sectionTitles.forEach(title => {
        title.addEventListener('click', function(e) {
            e.stopPropagation();
            const section = this.parentElement;
            
            // Toggle this section
            section.classList.toggle('expanded');
            
            // If expanding, collapse others
            if (section.classList.contains('expanded')) {
                navSections.forEach(otherSection => {
                    if (otherSection !== section) {
                        otherSection.classList.remove('expanded');
                    }
                });
                
                // Click the first subsection link
                const firstLink = section.querySelector('.nav-subsection a');
                if (firstLink) {
                    firstLink.click();
                }
            }
        });
    });
    
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
            const firstSection = document.querySelector('.nav-subsection a');
            if (firstSection) {
                firstSection.click();
            }
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

# Function to clean section names
clean_section_name() {
    local name=$1
    # Remove .md extension first
    name=${name%.md}
    # Remove leading numbers and dashes
    name=${name#[0-9][0-9]-}
    name=${name#[0-9]-}
    # Replace dashes and underscores with spaces
    name=${name//-/ }
    name=${name//_/ }
    # Capitalize first letter of each word
    name=$(echo "$name" | awk '{for(i=1;i<=NF;i++)sub(/./,toupper(substr($i,1,1)),$i)}1')
    # Fix common words
    name=${name/ And / and }
    name=${name/ Of / of }
    name=${name/ In / in }
    name=${name/ As / as }
    name=${name/ To / to }
    name=${name/ The / the }
    name=${name/ A / a }
    name=${name/ For / for }
    name=${name/ With / with }
    # Capitalize first word regardless
    name="$(echo ${name:0:1} | tr '[:lower:]' '[:upper:]')${name:1}"
    echo "$name"
}

# Start HTML file
cat > "$OUTPUT_DIR/index.html" << 'HTML'
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Synchronism: A Computational Framework for Pattern Dynamics</title>
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

# Process navigation - only show top 2 levels
for main_dir in "$SECTIONS_DIR"/*; do
    if [ -d "$main_dir" ]; then
        main_dirname=$(basename "$main_dir")
        main_clean_name=$(clean_section_name "$main_dirname")
        
        echo "                <li class=\"nav-section\">" >> "$OUTPUT_DIR/index.html"
        echo "                    <div class=\"nav-section-title\">$main_clean_name</div>" >> "$OUTPUT_DIR/index.html"
        echo "                    <div class=\"nav-subsection\">" >> "$OUTPUT_DIR/index.html"
        
        # Check if this directory has subdirectories or just files
        has_subdirs=false
        for item in "$main_dir"/*; do
            if [ -d "$item" ]; then
                has_subdirs=true
                break
            fi
        done
        
        if [ "$has_subdirs" = true ]; then
            # Process subdirectories as clickable items
            for sub_dir in "$main_dir"/*; do
                if [ -d "$sub_dir" ]; then
                    sub_dirname=$(basename "$sub_dir")
                    sub_clean_name=$(clean_section_name "$sub_dirname")
                    section_id="${main_dirname}-${sub_dirname}"
                    
                    echo "                        <a href=\"#$section_id\" data-section=\"$section_id\">$sub_clean_name</a>" >> "$OUTPUT_DIR/index.html"
                fi
            done
        else
            # If no subdirectories, this main section is itself clickable
            section_id="${main_dirname}"
            echo "                        <a href=\"#$section_id\" data-section=\"$section_id\">View Section</a>" >> "$OUTPUT_DIR/index.html"
        fi
        
        echo "                    </div>" >> "$OUTPUT_DIR/index.html"
        echo "                </li>" >> "$OUTPUT_DIR/index.html"
    fi
done

# Continue HTML
cat >> "$OUTPUT_DIR/index.html" << 'HTML'
            </ul>
            
            <!-- Download Links -->
            <div style="margin-top: 1rem; padding-top: 0.5rem; border-top: 1px solid var(--border-color);">
                <h3 style="font-size: 0.9rem; margin-bottom: 0.3rem; font-weight: 600;">Downloads</h3>
                <a href="../Synchronism_Whitepaper.pdf" style="display: block; margin-bottom: 0.2rem; font-size: 0.85rem; padding: 0.1rem 0;">📕 PDF Version</a>
                <a href="../Synchronism_Whitepaper_Complete.md" style="display: block; font-size: 0.85rem; padding: 0.1rem 0;">📄 Markdown Version</a>
            </div>
        </nav>
        
        <!-- Main Content Area -->
        <main class="main-content">
            <header style="margin-bottom: 1.5rem; padding: 1rem 0; border-bottom: 2px solid var(--border-color);">
                <h1 style="font-size: 1.8rem; margin-bottom: 0.25rem; padding-bottom: 0; border-bottom: none;">Synchronism: A Computational Framework for Pattern Dynamics</h1>
            </header>
            
            <!-- Content Sections (dynamically populated) -->
HTML

# Function to recursively combine all markdown files in a directory
combine_markdown_content() {
    local dir=$1
    local first_file=true
    
    # Process all markdown files in this directory
    for file in "$dir"/*.md; do
        if [ -f "$file" ] && [ "$(basename "$file")" != "index.md" ]; then
            if [ "$first_file" = false ]; then
                echo "<div class=\"section-divider\"></div>" >> "$OUTPUT_DIR/temp_section.html"
            fi
            md_to_html "$file" "$OUTPUT_DIR/temp_content.html"
            cat "$OUTPUT_DIR/temp_content.html" >> "$OUTPUT_DIR/temp_section.html"
            rm "$OUTPUT_DIR/temp_content.html"
            first_file=false
            echo "  ✓ Added $(basename "$file")"
        fi
    done
    
    # Recursively process subdirectories
    for subdir in "$dir"/*; do
        if [ -d "$subdir" ]; then
            if [ "$first_file" = false ]; then
                echo "<div class=\"section-divider\"></div>" >> "$OUTPUT_DIR/temp_section.html"
            fi
            combine_markdown_content "$subdir"
            first_file=false
        fi
    done
}

# Generate content sections
echo "Generating flattened content sections..."

for main_dir in "$SECTIONS_DIR"/*; do
    if [ -d "$main_dir" ]; then
        main_dirname=$(basename "$main_dir")
        
        # Check if this directory has subdirectories
        has_subdirs=false
        for item in "$main_dir"/*; do
            if [ -d "$item" ]; then
                has_subdirs=true
                break
            fi
        done
        
        if [ "$has_subdirs" = true ]; then
            # Process each subdirectory as a separate section
            for sub_dir in "$main_dir"/*; do
                if [ -d "$sub_dir" ]; then
                    sub_dirname=$(basename "$sub_dir")
                    section_id="${main_dirname}-${sub_dirname}"
                    
                    echo "Processing $main_dirname/$sub_dirname..."
                    
                    # Create section container
                    echo "            <section id=\"$section_id\" class=\"section\">" >> "$OUTPUT_DIR/index.html"
                    
                    # Create temp file for this section
                    > "$OUTPUT_DIR/temp_section.html"
                    
                    # Combine all content from this subdirectory
                    combine_markdown_content "$sub_dir"
                    
                    # Add to main HTML
                    cat "$OUTPUT_DIR/temp_section.html" >> "$OUTPUT_DIR/index.html"
                    rm "$OUTPUT_DIR/temp_section.html"
                    
                    echo "            </section>" >> "$OUTPUT_DIR/index.html"
                fi
            done
        else
            # Process the main directory itself as a section
            section_id="${main_dirname}"
            
            echo "Processing $main_dirname..."
            
            # Create section container
            echo "            <section id=\"$section_id\" class=\"section\">" >> "$OUTPUT_DIR/index.html"
            
            # Create temp file for this section
            > "$OUTPUT_DIR/temp_section.html"
            
            # Combine all content from this directory
            combine_markdown_content "$main_dir"
            
            # Add to main HTML
            cat "$OUTPUT_DIR/temp_section.html" >> "$OUTPUT_DIR/index.html"
            rm "$OUTPUT_DIR/temp_section.html"
            
            echo "            </section>" >> "$OUTPUT_DIR/index.html"
        fi
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

echo ""
echo "✅ Flattened web version created successfully!"
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