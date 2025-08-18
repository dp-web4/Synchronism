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
