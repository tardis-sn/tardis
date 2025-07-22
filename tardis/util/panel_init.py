import panel as pn
from tardis.util.environment import Environment
import logging

logger = logging.getLogger(__name__)

VALID_MODES = ["ssh_jh", "notebook", "vscode", "vscode_noipy"]
preferred_mode = None

def ssh_jh():
    """Initialize panel for JupyterHub (colab comms)"""
    print("Initializing panel with colab comms for JupyterHub")
    pn.extension(comms="colab")

def notebook():
    """Initialize panel for standard Jupyter notebook (default comms)"""
    print("Initializing panel with ipywidgets comms for Jupyter notebook")
    pn.extension(comms="ipywidgets")

def vscode():
    """Initialize panel for VSCode (ipywidgets comms)"""
    print("Initializing panel with ipywidgets comms for VSCode")
    pn.extension(comms="ipywidgets")

def vscode_noipy():
    """Initialize panel for VSCode without ipywidgets (vscode comms)"""
    print("Initializing panel with vscode comms for VSCode (no ipywidgets)")
    pn.extension(comms="vscode")

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
        print("Auto-detected JupyterHub environment")
        ssh_jh()
    elif Environment.is_vscode():
        print("Auto-detected VSCode environment")
        vscode()
    elif Environment.is_notebook():
        print("Auto-detected Jupyter notebook environment")
        notebook()
    else:
        print("Defaulting to default panel comms.")
        notebook()

