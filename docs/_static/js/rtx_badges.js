/**
 * RTX Navigation Badge System
 * Adds "RTX" badges to sidebar navigation links for pages containing RTX content
 */

(function() {
    'use strict';
    
    function addRtxBadges() {
        // Get RTX pages from the global variable injected by Sphinx
        const rtxPages = window.rtxPages || [];
        
        console.log('RTX Badge System: Found', rtxPages.length, 'RTX pages:', rtxPages);
        
        if (rtxPages.length === 0) {
            console.log('RTX Badge System: No RTX pages detected');
            return;
        }
        
        // Select all navigation links in the sidebar and main navigation
        const navLinks = document.querySelectorAll(
            '.toctree-l1 > a, .toctree-l2 > a, .toctree-l3 > a, ' +
            '.sidebar a, .sphinxsidebar a, nav a, .bd-sidebar a'
        );
        
        console.log('RTX Badge System: Found', navLinks.length, 'navigation links');
        
        let badgesAdded = 0;
        
        navLinks.forEach(function(link) {
            const href = link.getAttribute('href');
            if (!href) return;
            
            // Extract the page name from the href
            // Remove leading slashes, .html extension, and ../ prefixes
            let pageName = href
                .replace(/^\.\.\//, '')  // Remove leading ../
                .replace(/^\//, '')      // Remove leading /
                .replace(/\.html.*$/, '') // Remove .html and any query strings
                .trim();
            
            // Check if this page is in the RTX pages list
            const isRtxPage = rtxPages.some(function(rtxPage) {
                // Direct match
                if (pageName === rtxPage) return true;
                // Match last component (e.g., 'modular_architecture')
                if (pageName.split('/').pop() === rtxPage.split('/').pop()) return true;
                // Match if pageName ends with rtxPage
                if (pageName.endsWith(rtxPage)) return true;
                // Match if rtxPage ends with pageName
                if (rtxPage.endsWith(pageName)) return true;
                return false;
            });
            
            if (isRtxPage) {
                // Check if badge already exists (in case script runs multiple times)
                if (link.querySelector('.rtx-nav-pill')) {
                    return;
                }
                
                // Create and append the RTX badge
                const badge = document.createElement('span');
                badge.className = 'rtx-nav-pill';
                badge.textContent = 'RTX';
                badge.setAttribute('aria-label', 'Uses TARDIS-RTX Architecture');
                badge.setAttribute('title', 'This page covers TARDIS-RTX architecture concepts');
                
                link.appendChild(document.createTextNode(' '));
                link.appendChild(badge);
                badgesAdded++;
                
                console.log('RTX Badge System: Added badge to', pageName, '(href:', href + ')');
            }
        });
        
        console.log('RTX Badge System: Added', badgesAdded, 'badges');
    }
    
    // Run when DOM is ready
    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', addRtxBadges);
    } else {
        addRtxBadges();
    }
})();
