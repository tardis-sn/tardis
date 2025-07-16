import panel as pn
import os
import logging
from IPython import get_ipython

logger = logging.getLogger(__name__)

VALID_MODES = ["moria", "notebook", "vscode", "vscode_noipy"]
preferred_mode = None


def get_environment():
    """Determine the execution environment."""
    if any(x for x in ('VSCODE_PID', 'VSCODE') if x in os.environ):
        return 'vscode'
    return 'jupyter'


def is_notebook():
    """Check if running in a Jupyter notebook environment."""
    try:
        from ipykernel.zmqshell import ZMQInteractiveShell
        from IPython.core.interactiveshell import InteractiveShell

        shell = get_ipython()
        if shell is None:
            return False

        if isinstance(shell, ZMQInteractiveShell):
            return True
        elif isinstance(shell, InteractiveShell):
            return False
        else:
            return False
    except (ImportError, NameError):
        return False


def is_moria():
    """Check if running in Moria/JupyterHub environment."""
    return is_notebook() and any(key.startswith('JUPYTERHUB') for key in os.environ.keys())


def is_vscode():
    """Check if running in VSCode environment."""
    return any(x for x in ('VSCODE_PID', 'VSCODE', 'VSCODE_CWD') if x in os.environ)

def moria():
    """Initialize panel for Moria/JupyterHub (colab comms)"""
    print("Initializing panel with colab comms for Moria/JupyterHub")
    pn.extension(comms="colab")

def notebook():
    """Initialize panel for standard Jupyter notebook (default comms)"""
    print("Initializing panel with default comms for Jupyter notebook")
    pn.extension(comms="ipywidgets")

def vscode():
    """Initialize panel for VSCode (ipywidgets comms)"""
    print("Initializing panel with ipywidgets comms for VSCode")
    pn.extension(comms="ipywidgets")

def vscode_noipy():
    """Initialize panel for VSCode without ipywidgets (vscode comms)"""
    print("Initializing panel with vscode comms for VSCode (no ipywidgets)")
    pn.extension(comms="vscode")

def simple():
    """Initialize panel with minimal configuration"""
    print("Initializing panel with minimal configuration")
    pn.extension()

def auto():
    """Autodetect environment and initialize panel"""
    # Use preferred mode if set
    if preferred_mode:
        if preferred_mode not in VALID_MODES:
            raise ValueError(f"Invalid mode: {preferred_mode}. Valid modes are: {VALID_MODES}")
        print(f"Using preferred mode: {preferred_mode}")
        modes = {"moria": moria, "notebook": notebook, "vscode": vscode, "vscode_noipy": vscode_noipy}
        try:
            modes[preferred_mode]()
        except Exception:
            print("Preferred mode failed, falling back to simple mode")
            simple()
        return

    # Otherwise autodetect
    try:
        if is_moria():
            print("Autodetected Moria/JupyterHub environment")
            moria()
        elif is_vscode():
            print("Autodetected VSCode environment")
            vscode()
        elif is_notebook():
            print("Autodetected Jupyter notebook environment")
            notebook()
        else:
            print("Defaulting to default panel comms.")
            notebook()
    except Exception as e:
        print(f"Auto-detection failed ({e}), using simple panel initialization")
        simple()
