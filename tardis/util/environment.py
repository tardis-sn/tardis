import os
from enum import StrEnum
import logging
from IPython import get_ipython
logger = logging.getLogger(__name__)

class Environment(StrEnum):    
    VSCODE = 'vscode'
    JUPYTER = 'jupyter'
    TERMINAL = 'terminal'
    
    @classmethod
    def get_current(cls) -> 'Environment':
        """Get the current execution environment.
        
        Returns
        -------
        Environment
            The current environment enum member.
        """
        if any(x for x in ('VSCODE_PID', 'VSCODE') if x in os.environ):
            return cls.VSCODE
        elif cls.is_notebook():
            return cls.JUPYTER
        else:
            return cls.TERMINAL
    
    @staticmethod
    def is_notebook() -> bool:
        """
        Checking the shell environment where the simulation is run is Jupyter based

        Returns
        -------
        True : if the shell environment is Jupyter Based
        False : if the shell environment is Terminal or anything else
        """
        try:
          from ipykernel.zmqshell import ZMQInteractiveShell
          from IPython.core.interactiveshell import InteractiveShell
        except ImportError:
            logger.debug("Cannot import IPython/Jupyter modules. Not in IPython environment")
            return False
        
        try:
            # Trying to get the value of the shell via the get_ipython() method
            shell = get_ipython()
        except NameError:
            logger.debug("Cannot infer Shell Id")
            # Returns False if the shell name cannot be inferred correctly
            return False

        # Checking if the shell instance is Jupyter based & if True, returning True
        if isinstance(shell, ZMQInteractiveShell):
            return True
        # Checking if the shell instance is Terminal IPython based & if True, returning False
        elif isinstance(shell, InteractiveShell):
            return False
        # All other shell instances are returned False
        return False



