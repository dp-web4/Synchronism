#!/bin/bash

# make-web-clean.sh - Generate interactive web version from clean structure
# Usage: ./make-web-clean.sh

# Pull latest changes before building to avoid conflicts
echo "Checking for updates..."
git fetch

# Check if we're behind the remote
LOCAL=$(git rev-parse @)
REMOTE=$(git rev-parse @{u})
BASE=$(git merge-base @ @{u})

if [ $LOCAL = $REMOTE ]; then
    echo "Already up to date."
elif [ $LOCAL = $BASE ]; then
    echo "Pulling latest changes..."
    git pull
else
    echo "❌ Error: Your branch has diverged from the remote branch."
    echo "Please resolve conflicts manually before building:"
    echo "  1. Review changes with: git status"
    echo "  2. Either stash your changes: git stash"
    echo "  3. Or commit them: git add . && git commit -m 'your message'"
    echo "  4. Then pull: git pull"
    echo "  5. Run this script again"
    exit 1
fi

echo "Building Synchronism whitepaper web version from clean structure..."

OUTPUT_DIR="build/web-clean"
SECTIONS_DIR="sections"
ASSETS_DIR="$OUTPUT_DIR/assets"
WEB_SECTIONS_DIR="$OUTPUT_DIR/sections"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ASSETS_DIR"
mkdir -p "$WEB_SECTIONS_DIR"

echo "  ✓ Created output directories"

# Create main CSS file (using established dark grey theme)
cat > "$ASSETS_DIR/style.css" << 'CSS'
/* Synchronism Whitepaper Styles - Clean Structure */
:root {
    --primary-color: #374151;
    --secondary-color: #64748b;
    --background: #ffffff;
    --text-color: #1e293b;
    --border-color: #e2e8f0;
    --code-bg: #f8fafc;
    --sidebar-width: 300px;
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
    display: flex;
    flex-direction: column;
}

.sidebar h2 {
    font-size: 1.1rem;
    margin-bottom: 0.8rem;
    color: var(--primary-color);
    font-weight: 600;
}

/* Download Links */
.download-links {
    display: flex;
    flex-direction: column;
    gap: 0.5rem;
    margin-top: auto;
    padding: 0.75rem;
    background: white;
    border-radius: 8px;
    border: 1px solid var(--border-color);
    position: sticky;
    bottom: 1rem;
}

.download-links a {
    display: flex;
    align-items: center;
    padding: 0.5rem 0.75rem;
    text-decoration: none;
    color: var(--text-color);
    border-radius: 4px;
    transition: background-color 0.2s;
    font-size: 0.9rem;
}

.download-links a:hover {
    background-color: #f1f5f9;
    color: var(--primary-color);
}

.nav-list {
    list-style: none;
    flex: 1;
    overflow-y: auto;
    padding-bottom: 1rem;
}

.nav-section {
    margin-bottom: 0.5rem;
}

.nav-section-title {
    font-weight: 600;
    color: var(--primary-color);
    padding: 0.4rem 0.6rem;
    border-radius: 4px;
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

/* Only show arrow for sections with subsections */
.nav-section.has-subsections .nav-section-title::after {
    content: '▼';
    font-size: 0.7rem;
    transition: transform 0.3s ease;
}

.nav-section.has-subsections.collapsed .nav-section-title::after {
    transform: rotate(-90deg);
}

.nav-subsections {
    list-style: none;
    margin-left: 1rem;
    overflow: hidden;
    max-height: 500px;
    transition: max-height 0.3s ease;
}

.nav-section.collapsed .nav-subsections {
    max-height: 0;
}

.nav-link {
    display: block;
    padding: 0.3rem 0.6rem;
    color: var(--secondary-color);
    text-decoration: none;
    font-size: 0.9rem;
    border-radius: 4px;
    transition: all 0.2s ease;
}

.nav-link:hover {
    background: rgba(55, 65, 81, 0.05);
    color: var(--primary-color);
}

.nav-link.active {
    background: var(--primary-color);
    color: white;
}

/* Main Content */
.main-content {
    margin-left: var(--sidebar-width);
    flex: 1;
    padding: 2rem 3rem;
    max-width: 900px;
}

/* Content Sections */
.content-section {
    display: none;
    animation: fadeIn 0.3s ease;
}

.content-section.active {
    display: block;
}

@keyframes fadeIn {
    from { opacity: 0; }
    to { opacity: 1; }
}

/* Typography */
h1 {
    font-size: 2.5rem;
    margin-bottom: 1rem;
    color: var(--primary-color);
}

h2 {
    font-size: 2rem;
    margin-top: 2rem;
    margin-bottom: 1rem;
    color: var(--primary-color);
    border-bottom: 2px solid var(--border-color);
    padding-bottom: 0.5rem;
}

h3 {
    font-size: 1.5rem;
    margin-top: 1.5rem;
    margin-bottom: 0.8rem;
    color: var(--primary-color);
}

h4 {
    font-size: 1.2rem;
    margin-top: 1.2rem;
    margin-bottom: 0.6rem;
    color: var(--primary-color);
}

p {
    margin-bottom: 1rem;
}

/* Lists */
ul, ol {
    margin-left: 2rem;
    margin-bottom: 1rem;
}

li {
    margin-bottom: 0.3rem;
}

/* Code */
code {
    background: var(--code-bg);
    padding: 0.2rem 0.4rem;
    border-radius: 3px;
    font-family: 'Consolas', 'Monaco', monospace;
    font-size: 0.9em;
}

pre {
    background: var(--code-bg);
    padding: 1rem;
    border-radius: 5px;
    overflow-x: auto;
    margin-bottom: 1rem;
}

pre code {
    background: none;
    padding: 0;
}

/* Blockquotes */
blockquote {
    border-left: 4px solid var(--primary-color);
    padding-left: 1rem;
    margin: 1rem 0;
    font-style: italic;
    color: var(--secondary-color);
}

/* Links */
a {
    color: var(--primary-color);
    text-decoration: none;
}

a:hover {
    text-decoration: underline;
}

/* Key Concepts */
.key-concept {
    background: #f0f4f8;
    border-left: 4px solid var(--primary-color);
    padding: 1rem;
    margin: 1rem 0;
    border-radius: 0 5px 5px 0;
}

/* Navigation Hints */
.navigation-hints {
    margin-top: 2rem;
    padding-top: 1rem;
    border-top: 1px solid var(--border-color);
}

/* Responsive */
@media (max-width: 768px) {
    .sidebar {
        width: 100%;
        position: relative;
        height: auto;
        border-right: none;
        border-bottom: 1px solid var(--border-color);
    }
    
    .main-content {
        margin-left: 0;
        padding: 1rem;
    }
    
    .container {
        flex-direction: column;
    }
}
CSS

echo "  ✓ Created CSS file"

# Function to convert markdown to HTML
md_to_html() {
    local input=$1
    local output=$2
    
    if command -v pandoc &> /dev/null; then
        pandoc "$input" -f markdown -t html5 --no-highlight -o "$output" 2>/dev/null
    else
        # Fallback: basic conversion
        echo "<div class='content-section'>" > "$output"
        sed -e 's/^# \(.*\)$/<h1>\1<\/h1>/' \
            -e 's/^## \(.*\)$/<h2>\1<\/h2>/' \
            -e 's/^### \(.*\)$/<h3>\1<\/h3>/' \
            -e 's/^#### \(.*\)$/<h4>\1<\/h4>/' \
            -e 's/^\* \(.*\)$/<li>\1<\/li>/' \
            -e 's/^- \(.*\)$/<li>\1<\/li>/' \
            -e 's/^> \(.*\)$/<blockquote>\1<\/blockquote>/' \
            -e 's/`\([^`]*\)`/<code>\1<\/code>/g' \
            -e 's/\*\*\([^*]*\)\*\*/<strong>\1<\/strong>/g' \
            -e 's/\*\([^*]*\)\*/<em>\1<\/em>/g' \
            -e 's/^$/<p>/' \
            "$input" >> "$output"
        echo "</div>" >> "$output"
    fi
}

# Generate content sections
echo "Generating content sections..."

# Process each section
process_content() {
    local section_dir=$1
    local section_id=$2
    local output_file="$WEB_SECTIONS_DIR/section_${section_id}.html"
    
    echo "<section id='section-$section_id' class='content-section'>" > "$output_file"
    
    # Process all markdown files in the section
    for md_file in "$section_dir"/*.md; do
        if [ -f "$md_file" ] && [ "$(basename "$md_file")" != "index.md" ]; then
            md_to_html "$md_file" "$OUTPUT_DIR/temp.html"
            cat "$OUTPUT_DIR/temp.html" >> "$output_file"
        fi
    done
    
    # Process subdirectories
    for subdir in "$section_dir"/*; do
        if [ -d "$subdir" ]; then
            for md_file in "$subdir"/*.md; do
                if [ -f "$md_file" ] && [ "$(basename "$md_file")" != "index.md" ]; then
                    md_to_html "$md_file" "$OUTPUT_DIR/temp.html"
                    cat "$OUTPUT_DIR/temp.html" >> "$output_file"
                fi
            done
        fi
    done
    
    echo "</section>" >> "$output_file"
}

# Process all sections
section_counter=0
for section in "$SECTIONS_DIR"/*; do
    if [ -d "$section" ]; then
        section_name=$(basename "$section")
        echo "  Processing $section_name..."
        process_content "$section" "$section_counter"
        ((section_counter++))
    fi
done

# Generate Appendix B: Proposals
echo "  Processing Appendix B: Proposals..."
output_file="$WEB_SECTIONS_DIR/section_${section_counter}.html"
echo "<section id='section-$section_counter' class='content-section'>" > "$output_file"
echo "<h1>Appendix B: Current Proposals</h1>" >> "$output_file"
echo "<p><em>This appendix contains all active proposals for improvements to the Synchronism whitepaper. These are suggestions under review and not yet integrated into the main text.</em></p>" >> "$output_file"

# Find and add all proposals
proposal_files=$(find "$SECTIONS_DIR" -path "*/meta/proposals/*.md" -type f 2>/dev/null | sort)
if [ -n "$proposal_files" ]; then
    echo "<h2>Navigation</h2>" >> "$output_file"
    echo "<ul>" >> "$output_file"
    
    # Build navigation with anchor links
    current_section=""
    while IFS= read -r proposal_file; do
        rel_path="${proposal_file#$SECTIONS_DIR/}"
        section_path=$(echo "$rel_path" | sed 's|/meta/proposals/.*||')
        proposal_name=$(basename "$proposal_file" .md)
        # Create anchor ID from section path and proposal name
        anchor_id=$(echo "${section_path}-${proposal_name}" | sed 's/[^a-zA-Z0-9-]/-/g')
        
        if [ "$section_path" != "$current_section" ]; then
            if [ -n "$current_section" ]; then
                echo "</ul>" >> "$output_file"
            fi
            current_section="$section_path"
            echo "<li><strong>$section_path:</strong><ul>" >> "$output_file"
        fi
        echo "<li><a href=\"#$anchor_id\">$proposal_name</a></li>" >> "$output_file"
    done <<< "$proposal_files"
    echo "</ul></ul>" >> "$output_file"
    
    echo "<hr>" >> "$output_file"
    
    # Add actual proposals with anchor IDs
    current_section=""
    while IFS= read -r proposal_file; do
        rel_path="${proposal_file#$SECTIONS_DIR/}"
        section_path=$(echo "$rel_path" | sed 's|/meta/proposals/.*||')
        proposal_name=$(basename "$proposal_file" .md)
        # Create anchor ID matching the navigation
        anchor_id=$(echo "${section_path}-${proposal_name}" | sed 's/[^a-zA-Z0-9-]/-/g')
        
        if [ "$section_path" != "$current_section" ]; then
            current_section="$section_path"
            echo "<h4>Proposals for: $section_path</h4>" >> "$output_file"
        fi
        
        # Extract key info from the proposal file
        proposal_id=$(grep "\*\*ID\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
        title=$(grep -m1 "^#### Proposal" "$proposal_file" | sed 's/^#### Proposal [0-9]*: //')
        # Extract author - handle both formats: "Author" and "Author (LCT: hash)"
        author=$(grep "\*\*Author\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //' | sed 's/ (LCT:.*)$//')
        date=$(grep "\*\*Date\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
        status=$(grep "\*\*Status\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
        type=$(grep "\*\*Type\*\*:" "$proposal_file" | head -1 | sed 's/.*\*\*: //')
        
        # Compressed single-line header with anchor
        echo "<p id=\"$anchor_id\"><strong>${proposal_id}. ${title}</strong> — ${author} | ${date} | ${status} | ${type}</p>" >> "$output_file"
        
        # Extract and convert just the content sections
        {
            awk '
                /^###### Proposed Change/,/^###### Rationale/ { 
                    if (!/^######/) print 
                }
                /^###### Specific Text Changes/,/^###### Impact Assessment/ { 
                    if (!/^###### Impact Assessment/) {
                        gsub(/^######/, "#####", $0)
                        print
                    }
                }
            ' "$proposal_file"
        } > "$OUTPUT_DIR/temp_proposal.md"
        
        md_to_html "$OUTPUT_DIR/temp_proposal.md" "$OUTPUT_DIR/temp.html"
        cat "$OUTPUT_DIR/temp.html" >> "$output_file"
        rm -f "$OUTPUT_DIR/temp_proposal.md"
        echo "<hr>" >> "$output_file"
    done <<< "$proposal_files"
else
    echo "<p><em>No active proposals at this time.</em></p>" >> "$output_file"
fi

echo "</section>" >> "$output_file"
((section_counter++))

# Generate Appendix C: Change Log
echo "  Processing Appendix C: Change Log..."
output_file="$WEB_SECTIONS_DIR/section_${section_counter}.html"
echo "<section id='section-$section_counter' class='content-section'>" > "$output_file"
echo "<h1>Appendix C: Change Log</h1>" >> "$output_file"
echo "<p><em>Version history and evolution of the Synchronism whitepaper.</em></p>" >> "$output_file"

# Collect all changelog entries
for section in "$SECTIONS_DIR"/*; do
    if [ -d "$section" ]; then
        changelog_file="$section/meta/CHANGELOG.md"
        if [ -f "$changelog_file" ]; then
            section_name=$(basename "$section")
            echo "<h2>$section_name</h2>" >> "$output_file"
            
            # Process changelog, removing redundant "Section Changelog" and styling Format/Entries
            awk '
                /^#### Section Changelog/ { next }  # Skip redundant title
                /^###### Format/ { print "<h4><strong>Format</strong></h4>"; next }
                /^###### Entries/ { print "<h4><strong>Entries</strong></h4>"; next }
                { print }
            ' "$changelog_file" > "$OUTPUT_DIR/temp_changelog_clean.md"
            
            # Downgrade remaining headers before converting to HTML
            sed 's/^###/######/g; s/^##/#####/g; s/^#\([^#]\)/####\1/g' "$OUTPUT_DIR/temp_changelog_clean.md" > "$OUTPUT_DIR/temp_changelog.md"
            md_to_html "$OUTPUT_DIR/temp_changelog.md" "$OUTPUT_DIR/temp.html"
            cat "$OUTPUT_DIR/temp.html" >> "$output_file"
            rm -f "$OUTPUT_DIR/temp_changelog.md" "$OUTPUT_DIR/temp_changelog_clean.md"
            echo "<hr>" >> "$output_file"
        fi
    fi
done

echo "</section>" >> "$output_file"
((section_counter++))

# Clean up temp file
rm -f "$OUTPUT_DIR/temp.html"

echo "  ✓ Generated content sections"

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
        <nav class="sidebar">
            <h2>Synchronism</h2>
            <ul class="nav-list" id="nav-list">
                <!-- Navigation will be generated by JavaScript -->
            </ul>
            <div class="download-links">
                <a href="./Synchronism_Whitepaper.pdf" download>📄 Download PDF</a>
                <a href="./Synchronism_Whitepaper_Complete.md" download>📝 Download Markdown</a>
            </div>
        </nav>
        
        <main class="main-content" id="main-content">
            <!-- Content sections will be loaded here -->
        </main>
    </div>
    
    <script>
    // Navigation structure
    const sections = [
        {
            title: "Executive Summary",
            id: "executive-summary",
            file: "section_0.html"
        },
        {
            title: "Introduction",
            id: "introduction",
            file: "section_1.html"
        },
        {
            title: "Perspective",
            id: "perspective",
            file: "section_2.html"
        },
        {
            title: "Hermetic Principles",
            id: "hermetic-principles",
            file: "section_3.html"
        },
        {
            title: "Fundamental Concepts",
            id: "fundamental-concepts",
            file: "section_4.html",
            subsections: [
                "Universe Grid",
                "Time Slices",
                "Intent Transfer",
                "Emergence",
                "Field Effects",
                "Interaction Modes",
                "Coherence",
                "Markov Blankets",
                "Markov Relevancy",
                "Spectral Existence",
                "Abstraction",
                "Entity Interactions"
            ]
        },
        {
            title: "Quantum & Macro Phenomena",
            id: "quantum-macro",
            file: "section_5.html",
            subsections: [
                "CRT Analogy",
                "Superposition",
                "Wave-Particle",
                "Entanglement",
                "Witness Effect",
                "Relativity",
                "Speed Limits",
                "Decoherence",
                "Temperature",
                "Energy",
                "Universal Field",
                "Chemistry",
                "Life & Cognition",
                "Gravity",
                "Dark Matter",
                "Superconductivity",
                "Permeability",
                "Electromagnetic",
                "Energy Refinement",
                "Temperature Refinement",
                "Cognition Refinement",
                "String Theory"
            ]
        },
        {
            title: "Implications",
            id: "implications",
            file: "section_6.html",
            subsections: [
                "Unified Understanding",
                "Scientific Inquiry",
                "Ethical & Philosophical",
                "Open Questions"
            ]
        },
        {
            title: "Conclusion",
            id: "conclusion",
            file: "section_7.html"
        },
        {
            title: "Glossary",
            id: "glossary",
            file: "section_8.html"
        },
        {
            title: "Mathematical Appendix",
            id: "appendix-mathematical",
            file: "section_9.html"
        },
        {
            title: "Appendix B: Proposals",
            id: "appendix-proposals",
            file: "section_10.html"
        },
        {
            title: "Appendix C: Change Log",
            id: "appendix-changelog",
            file: "section_11.html"
        }
    ];
    
    // Content cache
    const contentCache = {};
    
    // Generate navigation
    function generateNavigation() {
        const navList = document.getElementById('nav-list');
        
        sections.forEach((section, index) => {
            const navSection = document.createElement('li');
            navSection.className = 'nav-section';
            
            // Add classes for sections with subsections
            if (section.subsections) {
                navSection.classList.add('has-subsections');
                navSection.classList.add('collapsed');
            }
            
            const title = document.createElement('div');
            title.className = 'nav-section-title';
            title.textContent = section.title;
            title.onclick = () => {
                if (section.subsections) {
                    // Collapse all other sections first
                    document.querySelectorAll('.nav-section').forEach(otherSection => {
                        if (otherSection !== navSection && otherSection.querySelector('.nav-subsections')) {
                            otherSection.classList.add('collapsed');
                        }
                    });
                    // Toggle current section
                    const wasCollapsed = navSection.classList.contains('collapsed');
                    navSection.classList.toggle('collapsed');
                    
                    // If we're expanding (was collapsed, now open), load first subsection
                    if (wasCollapsed) {
                        // Get the first subsection's anchor
                        const firstSub = section.subsections[0];
                        let anchorId = firstSub.toLowerCase().replace(/[&\s]/g, '-').replace(/--+/g, '-');
                        
                        // Apply special case mappings
                        const anchorMappings = {
                            'time-slices': 'time-as-planck-timed-slices',
                            'entity-interactions': 'entity-interaction-effects'
                        };
                        
                        if (anchorMappings[anchorId]) {
                            anchorId = anchorMappings[anchorId];
                        }
                        
                        // Load the section with first subsection anchor
                        loadSection(section.id, section.file, anchorId);
                        
                        // Highlight the first subsection link
                        const firstLink = navSection.querySelector('.nav-subsections .nav-link');
                        if (firstLink) {
                            highlightNav(firstLink);
                        }
                    }
                } else {
                    // Collapse all sections when clicking a non-expandable item
                    document.querySelectorAll('.nav-section').forEach(otherSection => {
                        if (otherSection.querySelector('.nav-subsections')) {
                            otherSection.classList.add('collapsed');
                        }
                    });
                    loadSection(section.id, section.file);
                    highlightNav(title);
                }
            };
            navSection.appendChild(title);
            
            if (section.subsections) {
                const subList = document.createElement('ul');
                subList.className = 'nav-subsections';
                
                section.subsections.forEach((sub, subIndex) => {
                    const subItem = document.createElement('li');
                    const subLink = document.createElement('a');
                    subLink.className = 'nav-link';
                    subLink.href = '#' + section.id;
                    subLink.textContent = sub;
                    // Generate anchor ID from subsection name
                    let anchorId = sub.toLowerCase().replace(/[&\s]/g, '-').replace(/--+/g, '-');
                    
                    // Special case mappings for mismatched IDs
                    const anchorMappings = {
                        'time-slices': 'time-as-planck-timed-slices',
                        'entity-interactions': 'entity-interaction-effects'
                    };
                    
                    if (anchorMappings[anchorId]) {
                        anchorId = anchorMappings[anchorId];
                    }
                    
                    subLink.onclick = (e) => {
                        e.preventDefault();
                        loadSection(section.id, section.file, anchorId);
                        highlightNav(subLink);
                    };
                    subItem.appendChild(subLink);
                    subList.appendChild(subItem);
                });
                
                navSection.appendChild(subList);
            }
            
            navList.appendChild(navSection);
        });
    }
    
    // Load section content
    async function loadSection(sectionId, fileName, anchorId) {
        const mainContent = document.getElementById('main-content');
        
        // Hide all existing sections
        document.querySelectorAll('.content-section').forEach(section => {
            section.style.display = 'none';
        });
        
        // Function to scroll to anchor after content loads
        const scrollToAnchor = () => {
            if (anchorId) {
                // Try to find element with matching ID
                let targetElement = document.getElementById(anchorId);
                
                // If not found, try to find heading with similar text
                if (!targetElement) {
                    const headings = mainContent.querySelectorAll('h1, h2, h3, h4');
                    for (let heading of headings) {
                        const headingId = heading.id || '';
                        const headingText = heading.textContent.toLowerCase().replace(/[&\s]/g, '-').replace(/--+/g, '-');
                        if (headingId.includes(anchorId) || headingText.includes(anchorId)) {
                            targetElement = heading;
                            break;
                        }
                    }
                }
                
                if (targetElement) {
                    setTimeout(() => {
                        targetElement.scrollIntoView({ behavior: 'smooth', block: 'start' });
                    }, 100);
                } else {
                    // Scroll to top if anchor not found
                    window.scrollTo(0, 0);
                }
            } else {
                // Scroll to top for main sections
                window.scrollTo(0, 0);
            }
        };
        
        // Check cache
        if (contentCache[fileName]) {
            mainContent.innerHTML = contentCache[fileName];
            // Make the loaded section visible
            const loadedSection = mainContent.querySelector('.content-section');
            if (loadedSection) {
                loadedSection.style.display = 'block';
            }
            scrollToAnchor();
        } else {
            try {
                const response = await fetch('sections/' + fileName);
                const content = await response.text();
                contentCache[fileName] = content;
                mainContent.innerHTML = content;
                // Make the loaded section visible
                const loadedSection = mainContent.querySelector('.content-section');
                if (loadedSection) {
                    loadedSection.style.display = 'block';
                }
                scrollToAnchor();
            } catch (error) {
                console.error('Error loading section:', error);
                mainContent.innerHTML = '<div class="content-section" style="display: block;"><p>Error loading content. Please run from a web server (not file://).</p></div>';
            }
        }
        
        // Update URL
        window.location.hash = anchorId ? sectionId + '-' + anchorId : sectionId;
    }
    
    // Highlight active navigation
    function highlightNav(activeLink) {
        document.querySelectorAll('.nav-link').forEach(link => {
            link.classList.remove('active');
        });
        if (activeLink) {
            activeLink.classList.add('active');
        }
    }
    
    // Initialize
    generateNavigation();
    
    // Load initial section
    loadSection('executive-summary', 'section_0.html');
    
    // Handle hash changes
    window.addEventListener('hashchange', () => {
        const hash = window.location.hash.slice(1);
        const section = sections.find(s => s.id === hash);
        if (section) {
            loadSection(section.id, section.file);
        }
    });
    </script>
</body>
</html>
HTML

echo "  ✓ Created main HTML file"

# Copy to docs directory for GitHub Pages
DOCS_DIR="../docs/whitepaper"
# Note: Don't create directory as it should already exist with MD and PDF

# Copy all web files including assets
cp -r "$OUTPUT_DIR"/* "$DOCS_DIR/"
echo "  ✓ Copied web files to docs"

# Explicitly ensure assets are copied
if [ -d "$OUTPUT_DIR/assets" ]; then
    cp -r "$OUTPUT_DIR/assets" "$DOCS_DIR/"
    echo "  ✓ Copied assets (navigation.js, style.css) to docs"
fi

# Also copy PDF and Markdown files if they exist
if [ -f "build/Synchronism_Whitepaper.pdf" ]; then
    cp "build/Synchronism_Whitepaper.pdf" "$DOCS_DIR/"
    echo "  ✓ Copied PDF to docs"
fi

if [ -f "build/Synchronism_Whitepaper_Complete.md" ]; then
    cp "build/Synchronism_Whitepaper_Complete.md" "$DOCS_DIR/"
    echo "  ✓ Copied Markdown to docs"
fi

echo "  ✓ Copied to GitHub Pages location"

echo ""
echo "✅ Web version generated successfully!"
echo "   Output: $OUTPUT_DIR/index.html"
echo "   GitHub Pages: ../docs/whitepaper/index.html"
echo ""
echo "To view locally, run:"
echo "   cd $OUTPUT_DIR && python3 -m http.server 8000"
echo "   Then open: http://localhost:8000"