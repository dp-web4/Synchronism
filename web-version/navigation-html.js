// HTML-based Synchronism Navigation
// Loads actual HTML files from the sections directory structure

class HTMLSynchronismNavigation {
    constructor() {
        this.currentSection = null;
        this.mainVersion = 'V0'; // Default version
        this.lastUpdated = null; // Will be loaded from JSON
        this.init();
    }

    init() {
        this.setupEventListeners();
        this.setupIntersectionObserver();
        this.loadInitialContent();
        this.updateActiveNavigation();
        
        // Load version info first, then update display
        this.loadVersionInfo();
    }

    setupEventListeners() {
        const navLinks = document.querySelectorAll('.nav-link');
        console.log(`Navigation: Found ${navLinks.length} nav links`);
        
        navLinks.forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const targetId = link.getAttribute('href').substring(1);
                console.log(`Navigation: Loading section ${targetId}`);
                this.loadSection(targetId);
            });
        });

        this.setupMobileMenu();
    }

    setupMobileMenu() {
        if (window.innerWidth <= 768) {
            const menuButton = document.createElement('button');
            menuButton.className = 'mobile-menu-toggle';
            menuButton.innerHTML = '☰';
            menuButton.addEventListener('click', () => {
                document.querySelector('.sidebar').classList.toggle('open');
            });
            document.querySelector('.content').prepend(menuButton);
        }
    }

    setupIntersectionObserver() {
        const observer = new IntersectionObserver((entries) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    this.currentSection = entry.target.id;
                    this.updateActiveNavigation();
                }
            });
        }, {
            rootMargin: '-20% 0px -60% 0px'
        });

        document.querySelectorAll('.content-section').forEach(section => {
            observer.observe(section);
        });
    }

    updateActiveNavigation() {
        document.querySelectorAll('.nav-link').forEach(link => {
            link.classList.remove('active');
        });

        if (this.currentSection) {
            const activeLink = document.querySelector(`a[href="#${this.currentSection}"]`);
            if (activeLink) {
                activeLink.classList.add('active');
            }
        }
    }

    loadInitialContent() {
        // Load introduction by default
        this.loadSection('introduction');
    }

    async loadSection(sectionId) {
        console.log(`Navigation: loadSection called for ${sectionId}`);
        const dynamicContent = document.getElementById('dynamic-content');
        if (!dynamicContent) {
            console.error('Navigation: dynamic-content element not found!');
            return;
        }

        // Hide the static about section when navigation starts
        const aboutSection = document.getElementById('about');
        if (aboutSection) {
            aboutSection.style.display = 'none';
        }

        // Clear existing content
        dynamicContent.innerHTML = '';

        // No status message needed - direct loading

        // Get the HTML file path for this section
        const htmlPath = this.getSectionHTMLPath(sectionId);
        console.log(`Navigation: HTML path for ${sectionId} = ${htmlPath}`);
        
        if (htmlPath) {
            try {
                const response = await fetch(htmlPath);
                if (response.ok) {
                    const html = await response.text();
                    // Set content to just the section HTML
                    dynamicContent.innerHTML = html;
                    
                    // Scroll to top of page
                    window.scrollTo(0, 0);
                    
                    // Setup math rendering with proper loading check
                    this.setupMathJax(dynamicContent);
                    
                    // Update active navigation
                    this.currentSection = sectionId;
                    this.updateActiveNavigation();
                    
                    console.log(`Navigation: Successfully loaded ${sectionId} (${html.length} chars)`);
                } else {
                    throw new Error(`Failed to load: ${response.status}`);
                }
            } catch (error) {
                console.error(`Navigation: Error loading section ${sectionId}:`, error);
                console.error(`Navigation: Failed path was: ${htmlPath}`);
                console.error(`Navigation: Error details:`, error.message);
                this.showSectionError(sectionId, dynamicContent, error, htmlPath);
            }
        } else {
            this.showSectionNotFound(sectionId, dynamicContent);
        }
    }

    getSectionHTMLPath(sectionId) {
        // Map section IDs to their HTML file paths
        const sectionPaths = {
            'about': 'sections/01-introduction/about.html',
            'introduction': 'sections/01-introduction/index.html',
            'perspective': 'sections/02-perspective/index.html',
            'hermetic-principles': 'sections/03-hermetic-principles/index.html',
            'fundamental-concepts-header': 'sections/04-fundamental-concepts/index.html',
            'universe-grid': 'sections/04-fundamental-concepts/01-universe-grid/index.html',
            'time-slices': 'sections/04-fundamental-concepts/02-time-slices/index.html',
            'intent-transfer': 'sections/04-fundamental-concepts/03-intent-transfer/index.html',
            'emergence-patterns': 'sections/04-fundamental-concepts/04-emergence/index.html',
            'field-effects': 'sections/04-fundamental-concepts/05-field-effects/index.html',
            'interaction-modes': 'sections/04-fundamental-concepts/06-interaction-modes/index.html',
            'coherence-feedback': 'sections/04-fundamental-concepts/07-coherence/index.html',
            'markov-blankets': 'sections/04-fundamental-concepts/08-markov-blankets/index.html',
            'mrh': 'sections/04-fundamental-concepts/09-mrh/index.html',
            'spectral-existence': 'sections/04-fundamental-concepts/10-spectral-existence/index.html',
            'abstraction': 'sections/04-fundamental-concepts/11-abstraction/index.html',
            'entity-interactions': 'sections/04-fundamental-concepts/12-entity-interactions/index.html',
            // Chapter 5: Quantum & Macro Phenomena
            'quantum-macro-header': 'sections/05-quantum-macro/index.html',
            'crt-analogy': 'sections/05-quantum-macro/01-crt-analogy/index.html',
            'quantum-superposition': 'sections/05-quantum-macro/02-superposition/index.html',
            'wave-particle': 'sections/05-quantum-macro/03-wave-particle/index.html',
            'entanglement': 'sections/05-quantum-macro/04-entanglement/index.html',
            'witness-effect': 'sections/05-quantum-macro/05-witness-effect/index.html',
            'relativity-view': 'sections/05-quantum-macro/06-relativity/index.html',
            'speed-limits': 'sections/05-quantum-macro/07-speed-limits/index.html',
            'macro-decoherence': 'sections/05-quantum-macro/08-decoherence/index.html',
            'temperature-phases': 'sections/05-quantum-macro/09-temperature/index.html',
            'energy': 'sections/05-quantum-macro/10-energy/index.html',
            'universal-field': 'sections/05-quantum-macro/11-universal-field/index.html',
            'chemistry': 'sections/05-quantum-macro/12-chemistry/index.html',
            'life-cognition': 'sections/05-quantum-macro/13-life-cognition/index.html',
            'gravity': 'sections/05-quantum-macro/14-gravity/index.html',
            'dark-matter': 'sections/05-quantum-macro/15-dark-matter/index.html',
            'superconductivity': 'sections/05-quantum-macro/16-superconductivity/index.html',
            'permeability': 'sections/05-quantum-macro/17-permeability/index.html',
            'electromagnetic': 'sections/05-quantum-macro/18-electromagnetic/index.html',
            'energy-refinement': 'sections/05-quantum-macro/19-energy-refinement/index.html',
            'temperature-refinement': 'sections/05-quantum-macro/20-temperature-refinement/index.html',
            'cognition-refinement': 'sections/05-quantum-macro/21-cognition-refinement/index.html',
            'string-theory': 'sections/05-quantum-macro/22-string-theory/index.html',
            // Chapter 6: Implications & Applications
            'unified-understanding': 'sections/06-implications/index.html',
            'unified-understanding-detail': 'sections/06-implications/01-unified-understanding/index.html',
            'scientific-inquiry': 'sections/06-implications/02-scientific-inquiry/index.html',
            'ethical-philosophical': 'sections/06-implications/03-ethical-philosophical/index.html',
            'open-questions': 'sections/06-implications/04-open-questions/index.html',
            // Chapter 7: Conclusion
            'conclusion': 'sections/07-conclusion/index.html',
            // Chapter 8: Glossary
            'glossary': 'sections/08-glossary/index.html',
            // Appendix A: Mathematical Foundations
            'appendix-mathematical': 'sections/09-appendix-mathematical/index.html',
            'appendix-a': 'sections/09-appendix-mathematical/index.html',
        };

        return sectionPaths[sectionId] || null;
    }

    getSectionTitle(sectionId) {
        const titles = {
            'about': 'About this Document',
            'introduction': '1. Introduction',
            'perspective': '2. Importance of Perspective', 
            'hermetic-principles': '3. Hermetic Principles',
            'fundamental-concepts-header': '4. Fundamental Concepts',
            'universe-grid': '4.1 Universe as Grid',
            'time-slices': '4.2 Time as Planck Slices',
            'intent-transfer': '4.3 Intent Transfer',
            'emergence-patterns': '4.4 Emergence & Patterns',
            'field-effects': '4.5 Field Effects',
            'interaction-modes': '4.6 Interaction Modes',
            'coherence-feedback': '4.7 Coherence & Feedback',
            'markov-blankets': '4.8 Markov Blankets',
            'mrh': '4.9 Markov Relevancy Horizon',
            'spectral-existence': '4.10 Spectral Existence',
            'abstraction': '4.11 Abstraction',
            'entity-interactions': '4.12 Entity Interactions'
        };
        return titles[sectionId] || sectionId;
    }

    showSectionNotFound(sectionId, container) {
        // Set content to just the not found message
        container.innerHTML = `
            <div style="background: #444; padding: 20px; margin: 20px 0; border-radius: 8px; color: #FFA500; text-align: center;">
                <h3>Section "${sectionId}" Not Yet Available</h3>
                <p>This section HTML file has not been created yet in the directory structure.</p>
                <p>Expected path: <code>${this.getSectionHTMLPath(sectionId) || 'Path not mapped'}</code></p>
                <p><a href="#introduction" style="color: #4CAF50;">← Return to Introduction</a></p>
            </div>
        `;
        // Scroll to top for error messages too
        window.scrollTo(0, 0);
    }

    setupMathJax(container) {
        // Check if MathJax is available and loaded
        if (window.MathJax && window.MathJax.typesetPromise) {
            try {
                MathJax.typesetPromise([container]).catch((err) => {
                    console.warn('MathJax typeset error:', err);
                });
            } catch (error) {
                console.warn('MathJax setup error:', error);
            }
        } else if (window.MathJax) {
            // MathJax exists but not fully loaded, wait for it
            const checkMathJax = () => {
                if (window.MathJax.typesetPromise) {
                    try {
                        MathJax.typesetPromise([container]).catch((err) => {
                            console.warn('MathJax delayed typeset error:', err);
                        });
                    } catch (error) {
                        console.warn('MathJax delayed setup error:', error);
                    }
                } else {
                    // Check again in 100ms
                    setTimeout(checkMathJax, 100);
                }
            };
            setTimeout(checkMathJax, 100);
        }
        // If MathJax doesn't exist at all, that's fine - equations will show as plain text
    }

    showSectionError(sectionId, container, error, htmlPath) {
        // Set content to just the error message
        container.innerHTML = `
            <div style="background: #660000; padding: 20px; margin: 20px 0; border-radius: 8px; color: #ffcccc; text-align: center;">
                <h3>🚨 Error Loading Section "${sectionId}"</h3>
                <p><strong>Error:</strong> ${error.message}</p>
                <p><strong>Path:</strong> <code>${htmlPath}</code></p>
                <p><strong>Possible causes:</strong></p>
                <ul style="text-align: left; display: inline-block;">
                    <li>Server not running (need http://localhost:8000)</li>
                    <li>File permissions issue</li>
                    <li>Network connectivity problem</li>
                    <li>CORS policy blocking request</li>
                </ul>
                <p><a href="#introduction" style="color: #4CAF50;">← Return to Introduction</a></p>
                <p style="font-size: 0.9em; margin-top: 15px;">Check browser console (F12) for more details</p>
            </div>
        `;
        // Scroll to top for error messages too  
        window.scrollTo(0, 0);
    }

    async loadVersionInfo() {
        console.log('Loading version info...');
        
        // Embedded version data to avoid CORS issues with file:// protocol
        // In production, this would fetch from main-version.json
        const versionData = {
            "mainVersion": "V0",
            "lastUpdated": "2025-06-27T16:45:00Z"
        };
        
        // Try to fetch from JSON file first (works in HTTP server mode)
        try {
            const response = await fetch('main-version.json');
            if (response.ok) {
                const jsonData = await response.json();
                this.mainVersion = jsonData.mainVersion || versionData.mainVersion;
                this.lastUpdated = jsonData.lastUpdated || versionData.lastUpdated;
                console.log('Loaded version from JSON file:', this.mainVersion, 'Last updated:', this.lastUpdated);
            } else {
                throw new Error('Fetch failed');
            }
        } catch (error) {
            // Fall back to embedded data (for file:// protocol)
            console.log('Using embedded version data (file:// mode)');
            this.mainVersion = versionData.mainVersion;
            this.lastUpdated = versionData.lastUpdated;
            console.log('Loaded version:', this.mainVersion, 'Last updated:', this.lastUpdated);
        }
        
        // Always update display after loading version info
        this.updateVersionDisplay();
    }

    updateVersionDisplay() {
        console.log('updateVersionDisplay called');
        console.log('this.mainVersion:', this.mainVersion);
        console.log('this.lastUpdated:', this.lastUpdated);
        
        // Format the stored lastUpdated timestamp from the JSON file
        let displayVersion = this.mainVersion;
        
        if (this.lastUpdated) {
            console.log('Processing lastUpdated timestamp...');
            try {
                const updateDate = new Date(this.lastUpdated);
                console.log('Parsed date:', updateDate);
                
                // Check if the date is valid
                if (!isNaN(updateDate.getTime())) {
                    const year = String(updateDate.getFullYear()).slice(-2);
                    const month = String(updateDate.getMonth() + 1).padStart(2, '0');
                    const day = String(updateDate.getDate()).padStart(2, '0');
                    const hours = String(updateDate.getHours()).padStart(2, '0');
                    const minutes = String(updateDate.getMinutes()).padStart(2, '0');
                    
                    displayVersion = `${this.mainVersion}.${year}.${month}.${day}.${hours}:${minutes}`;
                    console.log('Formatted version:', displayVersion);
                } else {
                    console.warn('Invalid date in lastUpdated field');
                    displayVersion = this.mainVersion;
                }
            } catch (error) {
                console.warn('Error parsing lastUpdated timestamp:', error);
                // Fall back to just the main version
                displayVersion = this.mainVersion;
            }
        } else {
            console.log('No lastUpdated timestamp found, using main version only');
        }
        
        // Update all version displays
        const versionElement = document.querySelector('.document-version');
        if (versionElement) {
            versionElement.textContent = displayVersion;
            console.log('Updated version display to:', displayVersion);
        } else {
            // If element not found, retry once after a short delay
            setTimeout(() => {
                const retryElement = document.querySelector('.document-version');
                if (retryElement) {
                    retryElement.textContent = displayVersion;
                    console.log('Updated version display on retry to:', displayVersion);
                }
            }, 100);
        }
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    // Check if we're running from file:// protocol
    if (window.location.protocol === 'file:') {
        console.warn('⚠️  Navigation requires a web server. Please run: python3 -m http.server 8000');
        const dynamicContent = document.getElementById('dynamic-content');
        if (dynamicContent) {
            dynamicContent.innerHTML = `
                <div style="background: #ff6b6b; color: white; padding: 20px; margin: 20px 0; border-radius: 8px; text-align: center;">
                    <h3>⚠️  Web Server Required</h3>
                    <p>The navigation system requires a local web server to function.</p>
                    <h4>To fix this:</h4>
                    <ol style="text-align: left; display: inline-block;">
                        <li>Open terminal/command prompt</li>
                        <li>Navigate to the web-version directory</li>
                        <li>Run: <code style="background: rgba(0,0,0,0.3); padding: 2px 6px; border-radius: 3px;">python3 -m http.server 8000</code></li>
                        <li>Open: <code style="background: rgba(0,0,0,0.3); padding: 2px 6px; border-radius: 3px;">http://localhost:8000</code></li>
                    </ol>
                    <p><strong>Current URL:</strong> ${window.location.href}</p>
                    <p><strong>Required URL:</strong> http://localhost:8000</p>
                </div>
            `;
        }
        return;
    }
    
    window.synchronismNav = new HTMLSynchronismNavigation();
});

// Handle browser navigation
window.addEventListener('popstate', (e) => {
    const hash = window.location.hash;
    if (hash && window.synchronismNav) {
        const sectionId = hash.substring(1);
        window.synchronismNav.loadSection(sectionId);
    }
});