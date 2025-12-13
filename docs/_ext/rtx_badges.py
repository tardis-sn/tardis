"""
RTX Badge System for TARDIS Documentation
Detects pages with RTX admonitions and injects badges into navigation.
"""
import json


def detect_rtx_pages(app, docname, source):
    """Detect pages containing RTX admonitions and store them."""
    # Initialize rtx_pages set if it doesn't exist
    if not hasattr(app.env, 'rtx_pages'):
        app.env.rtx_pages = set()
    
    # Check if the source contains RTX admonitions
    source_text = source[0] if source else ""
    if (':class: rtx' in source_text or 
        ':class: rtx\n' in source_text or
        '.. rtx-admonition::' in source_text):
        app.env.rtx_pages.add(docname)


def inject_rtx_pages(app, pagename, templatename, context, doctree):
    """Inject RTX pages list into HTML context as JavaScript variable."""
    if hasattr(app.env, 'rtx_pages'):
        # Convert set to list and add to context
        rtx_pages_list = sorted(list(app.env.rtx_pages))
        # Create inline script to inject the data
        rtx_script = f'<script>window.rtxPages = {json.dumps(rtx_pages_list)};</script>'
        
        # Add to metatags or body
        if 'metatags' not in context:
            context['metatags'] = ''
        context['metatags'] += rtx_script


def setup(app):
    """Setup the RTX badges extension."""
    app.connect("source-read", detect_rtx_pages)
    app.connect("html-page-context", inject_rtx_pages)
    
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
