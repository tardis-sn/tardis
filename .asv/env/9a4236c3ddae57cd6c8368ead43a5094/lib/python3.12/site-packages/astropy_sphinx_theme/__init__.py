""" Astropy Sphinx Theme """
import os

__version__ = "1.1"


def get_html_theme_path():
    """Return list of HTML theme paths."""
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    return [cur_dir]


def setup(app):
    app.add_html_theme('bootstrap-astropy', os.path.abspath(os.path.join(os.path.dirname(__file__), 'bootstrap-astropy')))
