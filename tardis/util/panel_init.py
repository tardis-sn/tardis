import panel as pn
from tardis.util.environment import Environment
import logging
import os

logger = logging.getLogger(__name__)

VALID_MODES = ["ssh_jh", "notebook", "vscode", "vscode_noipy"]
preferred_mode = None

# note: pn.extension(comms="ipywidgets") and pn.extension("ipywidgets") behave differently! 

def ssh_jh():
    """Initialize panel for JupyterHub (colab comms)"""
    pn.extension("tabulator", "plotly", comms="ipywidgets")

def notebook():
    """Initialize panel for standard Jupyter notebook (default comms)"""
    pn.extension("tabulator", "plotly", comms="ipywidgets")

def vscode():
    """Initialize panel for VSCode (ipywidgets comms)"""
    pn.extension("tabulator", "plotly", comms="ipywidgets")

def vscode_noipy():
    """Initialize panel for VSCode without ipywidgets (vscode comms)"""
    pn.extension("tabulator", "plotly", comms="vscode")

def sphinx():
    """Initialize panel for Sphinx documentation builds"""
    pn.extension("tabulator", "plotly", "ipywidgets")

def auto():
    """Auto-detect environment and initialize panel"""
    # Use preferred mode if set
    if preferred_mode:
        if preferred_mode not in VALID_MODES:
            raise ValueError(f"Invalid mode: {preferred_mode}. Valid modes are: {VALID_MODES}")
        print(f"Using preferred mode: {preferred_mode}")
        modes = {"ssh_jh": ssh_jh, "notebook": notebook, "vscode": vscode, "vscode_noipy": vscode_noipy}
        modes[preferred_mode]()
        return
    
    # Otherwise auto-detect
    if Environment.is_sshjh():
        ssh_jh()
    elif Environment.is_vscode():
        vscode()
    elif Environment.is_notebook():
        notebook()
    elif Environment.is_sphinx():
        sphinx()
    elif not Environment.is_terminal():
        pn.extension()

