"""
Version Badge System for TARDIS Documentation
Detects pages with version-specific admonitions (RTX/Classic) and injects badges into navigation.
"""
import json


def detect_version_pages(app, docname, source):
    """Detect pages containing RTX or Classic admonitions and store them."""
    # Initialize sets if they don't exist
    if not hasattr(app.env, 'rtx_pages'):
        app.env.rtx_pages = set()
    if not hasattr(app.env, 'classic_pages'):
        app.env.classic_pages = set()
    
    # Check if the source contains RTX or Classic admonitions
    source_text = source[0] if source else ""
    
    # Detect RTX pages
    if (':class: rtx' in source_text or 
        ':class: rtx\n' in source_text or
        '.. rtx-admonition::' in source_text):
        app.env.rtx_pages.add(docname)
    
    # Detect Classic pages
    if (':class: classic' in source_text or 
        ':class: classic\n' in source_text or
        '.. classic-admonition::' in source_text):
        app.env.classic_pages.add(docname)


def inject_version_pages(app, pagename, templatename, context, doctree):
    """Inject RTX and Classic pages lists into HTML context as JavaScript variables."""
    if 'metatags' not in context:
        context['metatags'] = ''
    
    if hasattr(app.env, 'rtx_pages'):
        # Convert set to list and add to context
        rtx_pages_list = sorted(list(app.env.rtx_pages))
        # Create inline script to inject the data
        rtx_script = f'<script>window.rtxPages = {json.dumps(rtx_pages_list)};</script>'
        context['metatags'] += rtx_script
    
    if hasattr(app.env, 'classic_pages'):
        # Convert set to list and add to context
        classic_pages_list = sorted(list(app.env.classic_pages))
        # Create inline script to inject the data
        classic_script = f'<script>window.classicPages = {json.dumps(classic_pages_list)};</script>'
        context['metatags'] += classic_script


def setup(app):
    """Setup the version badge system."""
    app.connect('source-read', detect_version_pages)
    app.connect('html-page-context', inject_version_pages)
    
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
