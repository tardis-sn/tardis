/**
 * TARDIS Version Badge System
 * Adds version-specific badges to sidebar navigation links:
 * - "RTX" (purple) for RTX architecture pages
 * - "Classic" (blue) for Classic architecture pages
 */

(function() {
    'use strict';
    
    const BADGE_CLASS_RTX = 'rtx-nav-pill';
    const BADGE_CLASS_CLASSIC = 'classic-nav-pill';
    let updateTimeout = null;

    function getBadgeHtml(type) {
        const isRtx = type === 'rtx';
        const className = isRtx ? BADGE_CLASS_RTX : BADGE_CLASS_CLASSIC;
        const text = isRtx ? ' RTX' : ' Classic';
        const label = isRtx ? 'Uses TARDIS-RTX Architecture' : 'Uses TARDIS Classic Architecture';
        
        const span = document.createElement('span');
        span.className = className;
        span.textContent = text;
        span.setAttribute('aria-label', label);
        span.setAttribute('title', label);
        return span;
    }

    function isPageMatch(pageName, pagesList) {
        if (!pagesList) return false;
        return pagesList.some(function(targetPage) {
            // Direct match
            if (pageName === targetPage) return true;
            // Match last component (e.g., 'modular_architecture')
            if (pageName.split('/').pop() === targetPage.split('/').pop()) return true;
            // Match if pageName ends with targetPage
            if (pageName.endsWith(targetPage)) return true;
            // Match if targetPage ends with pageName
            if (targetPage.endsWith(pageName)) return true;
            return false;
        });
    }

    function getCurrentPageName() {
        let path = window.location.pathname;
        // Remove leading slash
        path = path.replace(/^\//, '');
        // Remove .html
        path = path.replace(/\.html$/, '');
        return path;
    }

    function addVersionBadges() {
        // Get RTX and Classic pages from the global variables injected by Sphinx
        const rtxPages = window.rtxPages || [];
        const classicPages = window.classicPages || [];
        
        if (rtxPages.length === 0 && classicPages.length === 0) return;
        
        // Select all navigation links
        // Expanded selectors to catch more potential sidebar structures
        const navLinks = document.querySelectorAll(
            '.wy-nav-side a, .bd-sidebar a, .sphinxsidebar a, .sidebar a, ' +
            '.toctree-l1 > a, .toctree-l2 > a, .toctree-l3 > a'
        );
        
        const currentPageName = getCurrentPageName();

        navLinks.forEach(function(link) {
            // Skip if we already processed this link and it still has the badge
            if (link.querySelector('.' + BADGE_CLASS_RTX) || link.querySelector('.' + BADGE_CLASS_CLASSIC)) {
                return;
            }

            let href = link.getAttribute('href');
            let pageName = '';

            // Handle current page link which might be '#' or empty
            if (!href || href === '#') {
                // Check if this link represents the current page
                // Usually indicated by 'current' class on parent LI
                if (link.classList.contains('current') || 
                    (link.parentElement && link.parentElement.tagName === 'LI' && link.parentElement.classList.contains('current'))) {
                    pageName = currentPageName;
                } else {
                    return;
                }
            } else if (href.startsWith('javascript:')) {
                return;
            } else {
                // Extract the page name from the href
                pageName = href
                    .replace(/^\.\.\//, '')  // Remove leading ../
                    .replace(/^\//, '')      // Remove leading /
                    .replace(/\.html.*$/, '') // Remove .html and any query strings
                    .trim();
            }
            
            const isRtxPage = isPageMatch(pageName, rtxPages);
            const isClassicPage = isPageMatch(pageName, classicPages);
            
            if (isRtxPage) {
                link.appendChild(getBadgeHtml('rtx'));
            } else if (isClassicPage) {
                link.appendChild(getBadgeHtml('classic'));
            }
        });
    }
    
    function init() {
        addVersionBadges();
        
        // Observer for sidebar changes
        const sidebar = document.querySelector('.wy-nav-side') || 
                        document.querySelector('.bd-sidebar') || 
                        document.querySelector('.sphinxsidebar') ||
                        document.querySelector('.sidebar');

        if (sidebar) {
            const observer = new MutationObserver(function(mutations) {
                // Simple debounce to avoid performance issues
                if (updateTimeout) clearTimeout(updateTimeout);
                updateTimeout = setTimeout(addVersionBadges, 100);
            });
            
            observer.observe(sidebar, {
                childList: true,
                subtree: true,
                attributes: true,
                attributeFilter: ['class', 'style', 'href']
            });
        }
        
        // Hook into jQuery if available (common in Sphinx themes for AJAX nav)
        if (window.jQuery) {
            window.jQuery(document).ajaxComplete(function() {
                setTimeout(addVersionBadges, 100);
            });
        }
    }

    if (document.readyState === 'loading') {
        document.addEventListener('DOMContentLoaded', init);
    } else {
        init();
    }
    
    // Re-run on load just in case
    window.addEventListener('load', init);

})();