// Debug version of HTML-based Synchronism Navigation
// This version includes extensive console logging to help identify issues

class HTMLSynchronismNavigation {
    constructor() {
        console.log('🔧 DEBUG: Navigation constructor called');
        this.currentSection = null;
        this.init();
    }

    init() {
        console.log('🔧 DEBUG: Navigation init started');
        this.setupEventListeners();
        this.setupIntersectionObserver();
        this.loadInitialContent();
        this.updateActiveNavigation();
        console.log('🔧 DEBUG: Navigation init completed');
    }

    setupEventListeners() {
        console.log('🔧 DEBUG: Setting up event listeners');
        const navLinks = document.querySelectorAll('.nav-link');
        console.log(`🔧 DEBUG: Found ${navLinks.length} navigation links`);
        
        navLinks.forEach((link, index) => {
            const href = link.getAttribute('href');
            console.log(`🔧 DEBUG: Link ${index}: ${href}`);
            
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const targetId = link.getAttribute('href').substring(1);
                console.log(`🔧 DEBUG: Navigation click detected - targetId: ${targetId}`);
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
        console.log('🔧 DEBUG: Loading initial content');
        // Load introduction by default
        this.loadSection('introduction');
    }

    async loadSection(sectionId) {
        console.log(`🔧 DEBUG: loadSection called with sectionId: ${sectionId}`);
        
        const dynamicContent = document.getElementById('dynamic-content');
        if (!dynamicContent) {
            console.error('🔧 DEBUG ERROR: dynamic-content element not found!');
            return;
        }
        console.log('🔧 DEBUG: dynamic-content element found');

        // Clear existing content
        dynamicContent.innerHTML = '';
        console.log('🔧 DEBUG: Cleared existing content');

        // Add status indicator
        const statusDiv = document.createElement('div');
        statusDiv.className = 'modular-status';
        statusDiv.style.cssText = 'background: #2a2a2a; padding: 15px; margin: 20px 0; border-radius: 8px; color: #4CAF50; border-left: 4px solid #4CAF50;';
        statusDiv.innerHTML = `
            <h3 style="margin: 0 0 10px 0; color: #4CAF50;">🔧 DEBUG: HTML-Based Modular Framework</h3>
            <p style="margin: 0; opacity: 0.9;">Loading section from directory structure • <span style="color: #FFA500;">Section: ${this.getSectionTitle(sectionId)}</span></p>
            <p style="margin: 5px 0 0 0; font-size: 0.9em; color: #888;">Section ID: ${sectionId}</p>
        `;
        dynamicContent.appendChild(statusDiv);
        console.log('🔧 DEBUG: Added status indicator');

        // Get the HTML file path for this section
        const htmlPath = this.getSectionHTMLPath(sectionId);
        console.log(`🔧 DEBUG: HTML path resolved to: ${htmlPath}`);
        
        if (htmlPath) {
            try {
                console.log(`🔧 DEBUG: Attempting to fetch: ${htmlPath}`);
                const response = await fetch(htmlPath);
                console.log(`🔧 DEBUG: Fetch response status: ${response.status}`);
                
                if (response.ok) {
                    const html = await response.text();
                    console.log(`🔧 DEBUG: Successfully loaded HTML, length: ${html.length} characters`);
                    console.log(`🔧 DEBUG: HTML preview: ${html.substring(0, 100)}...`);
                    
                    // Replace status with the actual content
                    dynamicContent.innerHTML = statusDiv.outerHTML + html;
                    
                    // Setup math rendering
                    if (window.MathJax) {
                        console.log('🔧 DEBUG: MathJax found, rendering...');
                        MathJax.typesetPromise([dynamicContent]);
                    }
                    
                    // Update active navigation
                    this.currentSection = sectionId;
                    this.updateActiveNavigation();
                    
                    console.log(`🔧 DEBUG: Successfully loaded section: ${sectionId}`);
                } else {
                    throw new Error(`Failed to load: ${response.status} ${response.statusText}`);
                }
            } catch (error) {
                console.error(`🔧 DEBUG ERROR: Error loading section ${sectionId}:`, error);
                this.showSectionNotFound(sectionId, dynamicContent);
            }
        } else {
            console.log(`🔧 DEBUG: No HTML path found for section: ${sectionId}`);
            this.showSectionNotFound(sectionId, dynamicContent);
        }
    }

    getSectionHTMLPath(sectionId) {
        // Map section IDs to their HTML file paths
        const sectionPaths = {
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
        };

        const path = sectionPaths[sectionId] || null;
        console.log(`🔧 DEBUG: getSectionHTMLPath(${sectionId}) = ${path}`);
        return path;
    }

    getSectionTitle(sectionId) {
        const titles = {
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
        console.log(`🔧 DEBUG: Showing section not found for: ${sectionId}`);
        container.innerHTML += `
            <div style="background: #444; padding: 20px; margin: 20px 0; border-radius: 8px; color: #FFA500; text-align: center;">
                <h3>🔧 DEBUG: Section "${sectionId}" Not Yet Available</h3>
                <p>This section HTML file has not been created yet in the directory structure.</p>
                <p>Expected path: <code>${this.getSectionHTMLPath(sectionId) || 'Path not mapped'}</code></p>
                <p><a href="#introduction" style="color: #4CAF50;">← Return to Introduction</a></p>
            </div>
        `;
    }
}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    console.log('🔧 DEBUG: DOMContentLoaded event fired');
    window.synchronismNav = new HTMLSynchronismNavigation();
    console.log('🔧 DEBUG: Navigation instance created and assigned to window.synchronismNav');
});

// Handle browser navigation
window.addEventListener('popstate', (e) => {
    console.log('🔧 DEBUG: popstate event fired');
    const hash = window.location.hash;
    if (hash && window.synchronismNav) {
        const sectionId = hash.substring(1);
        console.log(`🔧 DEBUG: Loading section from URL hash: ${sectionId}`);
        window.synchronismNav.loadSection(sectionId);
    }
});