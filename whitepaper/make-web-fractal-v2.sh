#!/bin/bash

# make-web-fractal-v2.sh - Generate interactive web version from fractal structure (improved)
# Usage: ./make-web-fractal-v2.sh

echo "Building Synchronism whitepaper web version from fractal structure (v2)..."

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

.nav-list a {
    color: var(--text-color);
    text-decoration: none;
    padding: 0.2rem 0.6rem;
    display: block;
    border-radius: 4px;
    transition: all 0.2s;
    font-size: 0.9rem;
    line-height: 1.3;
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
    margin-top: 0;
    font-size: 0.85rem;
}

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

.nav-section > .nav-list {
    display: none;
}

.nav-section.expanded > .nav-list {
    display: block;
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
    const navSections = document.querySelectorAll('.nav-section');
    const sectionTitles = document.querySelectorAll('.nav-section-title');
    
    // Function to collapse all sections except the one containing the active link
    function updateExpandedSections(activeLink) {
        // Collapse all sections first
        navSections.forEach(section => {
            section.classList.remove('expanded');
        });
        
        // Expand the section containing the active link and its parents
        if (activeLink) {
            let parent = activeLink.closest('.nav-section');
            while (parent) {
                parent.classList.add('expanded');
                parent = parent.parentElement.closest('.nav-section');
            }
        }
    }
    
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
                updateExpandedSections(targetLink);
            }
            
            // Update URL hash
            window.location.hash = sectionId;
            
            // Scroll to top of content
            window.scrollTo(0, 0);
        }
    }
    
    // Add click handlers to section titles for expand/collapse and navigation
    sectionTitles.forEach(title => {
        title.addEventListener('click', function(e) {
            e.stopPropagation();
            const section = this.parentElement;
            
            // Check if section is already expanded
            const wasExpanded = section.classList.contains('expanded');
            
            // If section is collapsed, expand it
            if (!wasExpanded) {
                // Expand this section
                section.classList.add('expanded');
                
                // Collapse siblings at the same level
                const siblings = Array.from(section.parentElement.children)
                    .filter(child => child !== section && child.classList.contains('nav-section'));
                siblings.forEach(sibling => {
                    // Don't collapse if it contains the active link
                    if (!sibling.querySelector('.nav-link.active')) {
                        sibling.classList.remove('expanded');
                    }
                });
            }
            
            // Navigate to the first link in this section
            const firstLink = section.querySelector('.nav-link');
            if (firstLink) {
                const sectionId = firstLink.getAttribute('data-section');
                showSection(sectionId);
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
            showSection('00-executive-summary-00-executive-summary');
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
    # Fix specific issues
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

# Track processed files to avoid duplicates
declare -A processed_files

# Function to process directory for navigation (improved)
process_nav_directory() {
    local dir=$1
    local level=$2
    
    # First, process subdirectories
    for subdir in "$dir"/*/; do
        if [ -d "$subdir" ]; then
            dirname=$(basename "$subdir")
            clean_name=$(clean_section_name "$dirname")
            
            # Skip if this is just a container with a single file
            local md_files=()
            while IFS= read -r -d '' file; do
                md_files+=("$file")
            done < <(find "$subdir" -maxdepth 1 -name "*.md" -not -name "index.md" -print0)
            
            if [ ${#md_files[@]} -eq 1 ] && [ ! -d "$subdir"/*/ ]; then
                # Single file in directory - just show the file
                filename=$(basename "${md_files[0]}")
                file_clean_name=$(clean_section_name "$filename")
                section_id="${dirname}-${filename%.md}"
                
                if [ -z "${processed_files[$section_id]}" ]; then
                    processed_files[$section_id]=1
                    echo "                <li><a href=\"#$section_id\" class=\"nav-link\" data-section=\"$section_id\">$file_clean_name</a></li>" >> "$OUTPUT_DIR/index.html"
                fi
            else
                # Directory with multiple items or subdirectories
                echo "                <li class=\"nav-section\">" >> "$OUTPUT_DIR/index.html"
                echo "                    <div class=\"nav-section-title\">$clean_name</div>" >> "$OUTPUT_DIR/index.html"
                echo "                    <ul class=\"nav-list\">" >> "$OUTPUT_DIR/index.html"
                
                # Process files in this directory
                for file in "$subdir"/*.md; do
                    if [ -f "$file" ] && [ "$(basename "$file")" != "index.md" ]; then
                        filename=$(basename "$file")
                        file_clean_name=$(clean_section_name "$filename")
                        section_id="${dirname}-${filename%.md}"
                        
                        if [ -z "${processed_files[$section_id]}" ]; then
                            processed_files[$section_id]=1
                            echo "                        <li><a href=\"#$section_id\" class=\"nav-link\" data-section=\"$section_id\">$file_clean_name</a></li>" >> "$OUTPUT_DIR/index.html"
                        fi
                    fi
                done
                
                # Process subdirectories
                process_nav_directory "$subdir" $((level+1))
                
                echo "                    </ul>" >> "$OUTPUT_DIR/index.html"
                echo "                </li>" >> "$OUTPUT_DIR/index.html"
            fi
        fi
    done
}

# Generate navigation
echo "Generating navigation structure..."
process_nav_directory "$SECTIONS_DIR" 0

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
                <h1 style="font-size: 1.8rem; margin-bottom: 0.25rem; padding-bottom: 0; border-bottom: none;">Synchronism: A Comprehensive Model of Reality</h1>
            </header>
            
            <!-- Content Sections (dynamically populated) -->
HTML

# Reset processed files for content generation
unset processed_files
declare -A processed_files

# Function to process directory for content sections
process_content_directory() {
    local dir=$1
    local level=$2
    
    for subdir in "$dir"/*/; do
        if [ -d "$subdir" ]; then
            dirname=$(basename "$subdir")
            
            # Process all markdown files in this directory
            for file in "$subdir"/*.md; do
                if [ -f "$file" ] && [ "$(basename "$file")" != "index.md" ]; then
                    filename=$(basename "$file")
                    section_id="${dirname}-${filename%.md}"
                    
                    if [ -z "${processed_files[$section_id]}" ]; then
                        processed_files[$section_id]=1
                        
                        # Determine if this is the first section (executive summary)
                        if [[ "$section_id" == "00-executive-summary-00-executive-summary" ]]; then
                            echo "            <section id=\"$section_id\" class=\"section active\">" >> "$OUTPUT_DIR/index.html"
                        else
                            echo "            <section id=\"$section_id\" class=\"section\">" >> "$OUTPUT_DIR/index.html"
                        fi
                        
                        # Convert markdown to HTML and append
                        md_to_html "$file" "$OUTPUT_DIR/temp_section.html"
                        cat "$OUTPUT_DIR/temp_section.html" >> "$OUTPUT_DIR/index.html"
                        rm "$OUTPUT_DIR/temp_section.html"
                        
                        echo "            </section>" >> "$OUTPUT_DIR/index.html"
                        echo "  ✓ Processed $dirname/$(basename "$file")"
                    fi
                fi
            done
            
            # Process subdirectories recursively
            process_content_directory "$subdir" $((level+1))
        fi
    done
}

# Generate content sections
echo "Converting sections to HTML..."
process_content_directory "$SECTIONS_DIR" 0

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