"""
Page redirect system for TARDIS documentation
Creates redirect HTML files for moved pages.
"""
from pathlib import Path
from shutil import copyfile
import os


def to_html_ext(path):
    """Convert extension in the file path to .html"""
    return Path(path).with_suffix(".html")


def create_redirect_files(app, docname):
    """Create redirect html files at old paths specified in conf.py redirects list."""
    template_html_path = Path(app.srcdir) / "_templates/redirect_file.html"
    
    # Get redirects from config
    redirects = getattr(app.config, 'redirects', [])

    if app.builder.name == "html":
        for (old_fpath, new_fpath) in redirects:
            # Create a page redirection html file for old_fpath
            old_html_fpath = to_html_ext(Path(app.outdir) / old_fpath)
            old_html_fpath.parent.mkdir(parents=True, exist_ok=True)
            copyfile(template_html_path, old_html_fpath)

            # Replace url placeholders i.e. "#" in this file with the new url
            try:
                # Try using pathlib's relative_to for paths that share a common parent
                new_url = Path(to_html_ext(new_fpath)).relative_to(
                    Path(old_fpath).parent
                )
            except ValueError:
                # Fall back to os.path.relpath for arbitrary relative paths
                new_url = os.path.relpath(
                    to_html_ext(new_fpath), str(Path(old_fpath).parent)
                )
            # urls in a html file are relative to the dir containing it
            with open(old_html_fpath) as f:
                new_content = f.read().replace("#", str(new_url))
            with open(old_html_fpath, "w") as f:
                f.write(new_content)


def autodoc_skip_member(app, what, name, obj, skip, options):
    """Exclude specific functions/methods from the documentation"""
    exclusions = ("yaml_constructors", "yaml_implicit_resolvers")
    exclude = name in exclusions
    return skip or exclude


def setup(app):
    """Setup the redirects extension."""
    app.connect("autodoc-skip-member", autodoc_skip_member)
    app.connect("build-finished", create_redirect_files)
    
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
