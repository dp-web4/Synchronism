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
