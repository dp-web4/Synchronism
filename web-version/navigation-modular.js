// Modular Synchronism Navigation - Integration Test Version
// This loads our enhanced modular sections alongside legacy content

class ModularSynchronismNavigation {
    constructor() {
        this.currentSection = null;
        this.init();
    }

    init() {
        this.setupEventListeners();
        this.setupIntersectionObserver();
        this.loadModularContent();
        this.updateActiveNavigation();
    }

    setupEventListeners() {
        document.querySelectorAll('.nav-link').forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const targetId = link.getAttribute('href').substring(1);
                this.navigateToSection(targetId);
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
            const activeLink = document.querySelector(\`a[href="#\${this.currentSection}"]\`);
            if (activeLink) {
                activeLink.classList.add('active');
            }
        }
    }

    async navigateToSection(sectionId) {
        // Load only the selected section
        this.loadSingleSection(sectionId);
    }

    async loadModularContent() {
        // Initialize with default section (introduction)
        this.loadSingleSection('introduction');
        console.log('Modular content system initialized');
    }

    loadSingleSection(sectionId) {
        const dynamicContent = document.getElementById('dynamic-content');
        if (!dynamicContent) return;

        // Clear existing content
        dynamicContent.innerHTML = '';

        // Add modular status indicator
        const statusDiv = document.createElement('div');
        statusDiv.className = 'modular-status';
        statusDiv.style.cssText = 'background: #2a2a2a; padding: 15px; margin: 20px 0; border-radius: 8px; color: #4CAF50; border-left: 4px solid #4CAF50;';
        statusDiv.innerHTML = \`
            <h3 style="margin: 0 0 10px 0; color: #4CAF50;">✅ Enhanced Modular Framework Active</h3>
            <p style="margin: 0; opacity: 0.9;">Single-section view • <span style="color: #FFA500;">Section: \${this.getSectionTitle(sectionId)}</span></p>
        \`;
        dynamicContent.appendChild(statusDiv);

        // Load the specific section
        const content = this.getSectionContent(sectionId);
        if (content) {
            dynamicContent.innerHTML += content;
        } else {
            // Show section not found message
            dynamicContent.innerHTML += \`
                <div style="background: #444; padding: 20px; margin: 20px 0; border-radius: 8px; color: #FFA500; text-align: center;">
                    <h3>Section "\${sectionId}" Not Yet Available</h3>
                    <p>This section is being developed in the modular structure.</p>
                    <p><a href="#introduction" style="color: #4CAF50;">← Return to Introduction</a></p>
                </div>
            \`;
        }

        // Setup math rendering
        if (window.MathJax) {
            MathJax.typesetPromise([dynamicContent]);
        }

        // Update active navigation
        this.updateActiveNavigation();
        this.currentSection = sectionId;

        console.log(\`Loaded single section: \${sectionId}\`);
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

    getSectionContent(sectionId) {
        const sections = {
            'introduction': this.generateIntroductionContent(),
            'perspective': this.generatePerspectiveContent(),
            'hermetic-principles': this.generateHermeticPrinciplesContent(),
            'fundamental-concepts-header': this.generateFundamentalConceptsContent(),
            'universe-grid': this.generateUniverseGridContent(),
            'time-slices': this.generateTimeSlicesContent()
        };

        return sections[sectionId] || null;
    }


    generateIntroductionContent() {
        return \`
            <section id="introduction" class="content-section">
                <h2>1. Introduction</h2>
                <div class="section-content">
                    <p>Synchronism is a comprehensive model of reality that seeks to unify and transcend existing belief systems, including religions and scientific theories. It offers a perspective that aims to encompass all aspects of existence, providing a framework for understanding the universe and its phenomena.</p>
                    
                    <p>Developed through a series of philosophical inquiries and thought experiments, Synchronism attempts to bridge the gap between scientific understanding and spiritual or metaphysical concepts.</p>
                    
                    <p>In the Synchronism model, a key concept is that of "Intent," which serves as a reification of the abstract "greater force" proposed by various belief systems. This reification allows for the quantification and analysis of the underlying dynamics of reality. By representing the fundamental driving force of the universe as measurable "intent," Synchronism provides a framework that bridges scientific, philosophical, and spiritual perspectives, enabling a more unified approach to understanding existence.</p>
                    
                    <div class="key-concept">
                        <h3>Single Observer Model</h3>
                        <p>Synchronism is inherently a single-observer model, meaning that all phenomena are interpreted from the standpoint of a singular, unified observer. This perspective simplifies the complexities associated with multiple reference frames and provides a coherent narrative of intent transfer and emergence.</p>
                    </div>
                    
                    <div class="key-concept">
                        <h3>Mathematical Formalism</h3>
                        <p>For Synchronism to be a useful and relevant model, it is necessary to introduce formal mathematical treatments for its key concepts. The proposed mathematical formalism is introduced separately in <a href="#appendix-a">Appendix A</a>.</p>
                    </div>
                    
                    <div class="navigation-hints">
                        <h4>Continue Reading</h4>
                        <ul>
                            <li><a href="#perspective">Next: Importance of Perspective →</a></li>
                            <li><a href="#hermetic-principles">Understanding the Hermetic Foundation →</a></li>
                            <li><a href="#fundamental-concepts-header">Explore Core Concepts →</a></li>
                        </ul>
                    </div>
                </div>
            </section>\`;
    }

    generatePerspectiveContent() {
        return \`
            <section id="perspective" class="content-section">
                <h2>2. Importance of Perspective</h2>
                <div class="section-content">
                    <p>The significance of perspective in understanding reality is illustrated by the "Six Blind Men and the Elephant" analogy, an ancient parable that highlights the limitations of individual perception and the importance of holistic understanding.</p>
                    
                    <div class="key-concept">
                        <h3>The Six Blind Men and the Elephant</h3>
                        <p>In this story, six blind men encounter an elephant for the first time. Each man touches a different part of the elephant and describes what he believes the elephant to be based on his limited experience:</p>
                        <ul>
                            <li>The man who feels the leg says the elephant is like a pillar</li>
                            <li>The one who touches the tail describes it as a rope</li>
                            <li>The man who feels the trunk thinks it's like a tree branch</li>
                            <li>The one who touches the ear believes it's like a hand fan</li>
                            <li>The man who feels the belly describes it as a wall</li>
                            <li>The one who touches the tusk thinks it's like a solid pipe</li>
                        </ul>
                    </div>
                    
                    <h3>Key Insights from the Analogy</h3>
                    <p>This analogy illustrates several key points:</p>
                    <ul>
                        <li>Different witnesses may experience only parts of a whole, leading to incomplete or inaccurate conclusions</li>
                        <li>Consensus doesn't necessarily lead to truth, as all the men might agree on certain aspects while still missing the full picture</li>
                        <li>A comprehensive understanding requires both broadening one's perspective and gaining detailed knowledge</li>
                    </ul>
                    
                    <div class="key-concept">
                        <h3>Witness and Experience Framework</h3>
                        <p>Synchronism introduces formal concepts of <strong>Witness</strong> and <strong>Experience</strong>, defined as interactions of an entity within its fractal scale and levels of abstraction. This provides a framework for choosing the optimal scale and abstraction for analysis.</p>
                    </div>
                    
                    <div class="navigation-hints">
                        <h4>Continue Reading</h4>
                        <ul>
                            <li><a href="#hermetic-principles">Next: Relation to Hermetic Principles →</a></li>
                            <li><a href="#mrh">Understanding Markov Relevancy Horizon →</a></li>
                        </ul>
                    </div>
                </div>
            </section>\`;
    }

    generateHermeticPrinciplesContent() {
        return \`
            <section id="hermetic-principles" class="content-section">
                <h2>3. Relation to Hermetic Principles</h2>
                <div class="section-content">
                    <p>Synchronism aligns with and expands upon the seven Hermetic principles, offering a scientifically-grounded explanation of their manifestation in the universe.</p>
                    
                    <div class="hermetic-principle">
                        <h3>3.1 Mentalism: "The All is Mind; the Universe is Mental."</h3>
                        <p><strong>Synchronism:</strong> The universe is modeled as a vast array of interconnected "cells" that contain and transfer "intent," which can be seen as a form of mental energy or information.</p>
                        <p><em>Related:</em> <a href="#universe-grid">Universe as Grid</a>, <a href="#intent-transfer">Intent Transfer</a></p>
                    </div>
                    
                    <div class="hermetic-principle">
                        <h3>3.2 Correspondence: "As above, so below; as below, so above."</h3>
                        <p><strong>Synchronism:</strong> The model proposes a fractal nature of reality, where patterns repeat at different scales.</p>
                        <p><em>Related:</em> <a href="#emergence-patterns">Emergence</a>, <a href="#abstraction">Abstraction</a></p>
                    </div>
                    
                    <div class="hermetic-principle">
                        <h3>3.3 Vibration: "Nothing rests; everything moves and vibrates."</h3>
                        <p><strong>Synchronism:</strong> All phenomena are described as patterns of intent transfer, creating oscillations at the fundamental level.</p>
                        <p><em>Related:</em> <a href="#time-slices">Time as Planck Slices</a></p>
                    </div>
                    
                    <div class="hermetic-principle">
                        <h3>3.4 Polarity: "Everything is dual; everything has poles."</h3>
                        <p><strong>Synchronism:</strong> Resonance and dissonance reflect polarity, where entities can reinforce or weaken each other.</p>
                        <p><em>Related:</em> <a href="#interaction-modes">Interaction Modes</a></p>
                    </div>
                    
                    <div class="hermetic-principle">
                        <h3>3.5 Rhythm: "Everything flows, out and in; everything has its tides."</h3>
                        <p><strong>Synchronism:</strong> The universe evolves through discrete "ticks" of time, creating rhythmic progression.</p>
                    </div>
                    
                    <div class="hermetic-principle">
                        <h3>3.6 Cause and Effect: "Every cause has its effect; every effect has its cause."</h3>
                        <p><strong>Synchronism:</strong> Intent transfer creates causal chains, with each state informing the next.</p>
                    </div>
                    
                    <div class="hermetic-principle">
                        <h3>3.7 Gender: "Gender is in everything; everything has its masculine and feminine principles."</h3>
                        <p><strong>Synchronism:</strong> Fundamental duality within entities, similar to Generative Adversarial Networks (GANs) - a generative principle that proposes new patterns and a discriminative principle that evaluates them.</p>
                    </div>
                    
                    <div class="navigation-hints">
                        <h4>Continue Reading</h4>
                        <ul>
                            <li><a href="#fundamental-concepts-header">Next: Fundamental Concepts →</a></li>
                            <li><a href="#universe-grid">Explore the Universe as Grid →</a></li>
                        </ul>
                    </div>
                </div>
            </section>\`;
    }

    generateFundamentalConceptsContent() {
        return \`
            <section id="fundamental-concepts-header" class="content-section">
                <h2>4. Fundamental Concepts of Synchronism</h2>
                <div class="section-content">
                    <p>This chapter introduces the core concepts that form the foundation of the Synchronism model. These concepts build upon each other to create a comprehensive framework for understanding reality from the Planck scale to cosmic phenomena.</p>
                    
                    <div class="concept-map">
                        <h3>Core Structure (4.1-4.2)</h3>
                        <ul>
                            <li><a href="#universe-grid">4.1 Universe as a Grid of Planck Cells</a> - The spatial foundation</li>
                            <li><a href="#time-slices">4.2 Time as Planck-Timed Slices</a> - The temporal foundation</li>
                        </ul>
                        
                        <h3>Dynamic Mechanisms (4.3-4.5)</h3>
                        <ul>
                            <li><a href="#intent-transfer">4.3 Intent Transfer and Tension</a> - The fundamental driving force</li>
                            <li><a href="#emergence-patterns">4.4 Emergence and Patterns</a> - How complexity arises</li>
                            <li><a href="#field-effects">4.5 Emergent Properties and Field Effects</a> - Higher-order phenomena</li>
                        </ul>
                        
                        <h3>Organizational Principles (4.6+)</h3>
                        <ul>
                            <li><a href="#interaction-modes">4.6 Interaction Modes</a> - Resonance, dissonance, and indifference</li>
                            <li><a href="#coherence-feedback">4.7 Coherence and Feedback</a> - Stability mechanisms</li>
                            <li><a href="#markov-blankets">4.8 Markov Blankets</a> - Entity boundaries</li>
                            <li><a href="#mrh">4.9 Markov Relevancy Horizon</a> - Optimal analysis scope</li>
                            <li><a href="#spectral-existence">4.10 Spectral Existence</a> - Degrees of existence</li>
                            <li><a href="#abstraction">4.11 Abstraction</a> - Managing complexity</li>
                            <li><a href="#entity-interactions">4.12 Entity Interaction Effects</a> - How entities influence each other</li>
                        </ul>
                    </div>
                    
                    <div class="navigation-hints">
                        <h4>Begin Exploration</h4>
                        <ul>
                            <li><a href="#universe-grid">Start with 4.1: Universe as Grid →</a></li>
                            <li><a href="#intent-transfer">Jump to 4.3: Intent Transfer →</a></li>
                            <li><a href="#coherence-feedback">Explore 4.7: Coherence →</a></li>
                        </ul>
                    </div>
                </div>
            </section>\`;
    }

    generateUniverseGridContent() {
        return \`
            <section id="universe-grid" class="content-section">
                <h2>4.1 Universe as a Grid of Planck Cells</h2>
                <div class="section-content">
                    <p>Synchronism proposes that the universe can be modeled as an infinite three-dimensional grid of discrete cells. This concept provides a fundamental structure for understanding the nature of space and the interactions that occur within it.</p>
                    
                    <div class="key-concept">
                        <h3>Grid Structure</h3>
                        <p>Key aspects of this grid model include:</p>
                        <ul>
                            <li>Each cell is the size of a Planck length (≈ 1.616 × 10⁻³⁵ meters) in each dimension</li>
                            <li>The grid extends infinitely in all directions</li>
                            <li>Each cell contains a quantized amount of "<a href="#intent-transfer">intent</a>"</li>
                            <li>Intent in each cell is limited by a saturation maximum</li>
                        </ul>
                    </div>
                    
                    <div class="mathematical-note">
                        <h3>Mathematical Foundation</h3>
                        <p>This discrete spatial structure enables:</p>
                        <ul>
                            <li><strong>Precise Location Definition:</strong> Grid coordinates for every point</li>
                            <li><strong>Quantized Interactions:</strong> All phenomena in discrete units</li>
                            <li><strong>Intent Conservation:</strong> Total intent precisely tracked</li>
                            <li><strong>Transfer Mechanics:</strong> Intent moves only between adjacent cells</li>
                        </ul>
                        <p><em>See <a href="#basic-intent">Appendix A.1: Basic Intent Transfer</a> for mathematical details.</em></p>
                    </div>
                    
                    <div class="key-concept">
                        <h3>Connection to Hermetic Principles</h3>
                        <p>The grid embodies several <a href="#hermetic-principles">Hermetic principles</a>:</p>
                        <ul>
                            <li><strong>Mentalism:</strong> Grid as vast neural network</li>
                            <li><strong>Correspondence:</strong> Same structure across all scales</li>
                            <li><strong>Vibration:</strong> Intent patterns create oscillations</li>
                        </ul>
                    </div>
                    
                    <div class="practical-applications">
                        <h3>Understanding Through Analogy</h3>
                        <ul>
                            <li><strong>3D Cellular Automaton:</strong> Like Conway's Game of Life in 3D</li>
                            <li><strong>Computer Memory:</strong> Each cell stores and processes information</li>
                            <li><strong>Neural Network:</strong> Cells communicate with neighbors</li>
                            <li><strong>Spacetime Fabric:</strong> Grid provides the fabric of reality</li>
                        </ul>
                    </div>
                    
                    <div class="navigation-hints">
                        <h4>Continue Exploring</h4>
                        <ul>
                            <li><a href="#time-slices">Next: 4.2 Time as Planck Slices →</a></li>
                            <li><a href="#intent-transfer">Jump to: 4.3 Intent Transfer →</a></li>
                            <li><a href="#speed-limits">Application: Speed Limits →</a></li>
                        </ul>
                    </div>
                </div>
            </section>\`;
    }

    generateTimeSlicesContent() {
        return \`
            <section id="time-slices" class="content-section">
                <h2>4.2 Time as Planck-Timed Slices</h2>
                <div class="section-content">
                    <p>In the Synchronism model, time is the fundamental substrate of reality itself. Time progresses as discrete moments or "ticks," with each tick bringing forth a new slice of reality.</p>
                    
                    <div class="key-concept">
                        <h3>Time as Fundamental Substrate</h3>
                        <p>Time is the driving force behind all existence. The cessation of time implies cessation of all existence, as nothing can manifest without the passage of time. Time is the universal "Mind" that governs and sustains the universe's evolution.</p>
                    </div>
                    
                    <div class="time-structure">
                        <h3>Discrete Time Model</h3>
                        <ul>
                            <li><strong>Quantized Progression:</strong> Time advances in discrete "ticks" of Planck time (≈ 5.39 × 10⁻⁴⁴ seconds)</li>
                            <li><strong>Universal Slices:</strong> Each tick represents a complete snapshot of intent distribution across all cells</li>
                            <li><strong>Static Slices:</strong> Each slice is fixed and unchanging</li>
                            <li><strong>Causal Chain:</strong> Each slice is informed by all preceding states</li>
                        </ul>
                    </div>
                    
                    <div class="mathematical-note">
                        <h3>Mathematical Representation</h3>
                        <ul>
                            <li><strong>State Functions:</strong> Universe state U(t) at time t</li>
                            <li><strong>Transition Rules:</strong> U(t+1) = F(U(t))</li>
                            <li><strong>Deterministic Evolution:</strong> Future states determined by current state</li>
                            <li><strong>Conservation Laws:</strong> Total intent preserved across transitions</li>
                        </ul>
                    </div>
                    
                    <div class="philosophical-implications">
                        <h3>Philosophical Implications</h3>
                        <h4>Determinism and Free Will</h4>
                        <p>While each slice is determined by the preceding one, the universe is deterministic in structure but probabilistic in experience, since the universe cannot predict its next state until it "lives" it.</p>
                    </div>
                    
                    <div class="practical-understanding">
                        <h3>Understanding Through Analogies</h3>
                        <ul>
                            <li><strong>Film Frames:</strong> Static slices create illusion of motion</li>
                            <li><strong>Computer Clock Cycles:</strong> Synchronized component updates</li>
                            <li><strong>Universal Heartbeat:</strong> Each tick drives the universe forward</li>
                        </ul>
                    </div>
                    
                    <div class="navigation-hints">
                        <h4>Continue Exploring</h4>
                        <ul>
                            <li><a href="#intent-transfer">Next: 4.3 Intent Transfer →</a></li>
                            <li><a href="#speed-limits">Application: Speed Limits →</a></li>
                            <li><a href="#free-will-determinism">Philosophy: Free Will Questions →</a></li>
                        </ul>
                    </div>
                </div>
            </section>\`;
    }

}

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', () => {
    window.synchronismNav = new ModularSynchronismNavigation();
});

// Handle browser navigation
window.addEventListener('popstate', (e) => {
    const hash = window.location.hash;
    if (hash && window.synchronismNav) {
        const sectionId = hash.substring(1);
        window.synchronismNav.navigateToSection(sectionId);
    }
});