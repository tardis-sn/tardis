import pandas as pd
import numpy as np
import copy
import tardis

from tardis.io.configuration.config_reader import Configuration
from tardis.io.atom_data import AtomData
from tardis.model import SimulationState


def _set_tardis_config_property(tardis_config, key, value):
    """
    Sets tardis_config[key] = value for a TARDIS config
    dictionary. Recursively searches for the deepest
    dictionary and sets the corresponding value.

    Parameters
    ----------
    tardis_config : tardis.io.config_reader.Configuration
        The config where the key,value pair should be set.
    key : str
        The key in the tardis_config. Should be of the format
        'property.property.property'
    value : object
        The value to assign to the config dict for the
        corresponding key.
    """
    keyitems = key.split(".")
    tmp_dict = getattr(tardis_config, keyitems[0])
    for key in keyitems[1:-1]:
        tmp_dict = getattr(tmp_dict, key)
    setattr(tmp_dict, keyitems[-1], value)
    return


class tardisGrid:
    """
    A class that stores a grid of TARDIS parameters and
    facilitates running large numbers of simulations
    easily.

    Parameters
    ----------
    configFile : str or dict
        path to TARDIS yml file, or a pre-validated config dictionary.
    gridFrame : pandas.core.frame.DataFrame
        dataframe where each row is a set of parameters for
        a TARDIS simulation.

    Attributes
    ----------
    config : tardis.io.config_reader.Configuration
        The validated config dict read from the user provided
        configFile. This provides the base properties for the
        TARDIS simulation which the rows of the grid modify.
    grid : pandas.core.frame.DataFrame
        Dataframe where each row is a set of parameters for
        a TARDIS simulation.
    """

    def __init__(self, configFile, gridFrame):
        try:
            tardis_config = Configuration.from_yaml(configFile)
        except TypeError:
            tardis_config = Configuration.from_config_dict(configFile)

        self.config = tardis_config
        self.grid = gridFrame

        return

    def grid_row_to_config(self, row_index):
        """
        Converts a grid row to a TARDIS config dict.
        Modifies the base self.config according to the
        row in self.grid accessed at the provided row_index.
        Returns a deep copy so that the base config is
        not changed.

        Parameters
        ----------
        row_index : int
            Row index in grid.

        Returns
        -------
        tmp_config : tardis.io.config_reader.Configuration
            Deep copy of the base self.config with modified
            properties according to the selected row in the grid.
        """
        tmp_config = copy.deepcopy(self.config)
        grid_row = self.grid.iloc[row_index]
        for colname, value in zip(self.grid.columns, grid_row.values):
            _set_tardis_config_property(tmp_config, colname, value)
        return tmp_config

    def grid_row_to_simulation_state(self, row_index, atomic_data):
        """
        Generates a TARDIS SimulationState object using the base
        self.config modified by the specified grid row.

        Parameters
        ----------
        row_index : int
            Row index in grid.

        Returns
        -------
        model : tardis.model.base.SimulationState
        """
        rowconfig = self.grid_row_to_config(row_index)
        simulation_state = SimulationState.from_config(
            rowconfig, atom_data=atomic_data
        )
        return simulation_state

    def run_sim_from_grid(self, row_index, **tardiskwargs):
        """
        Runs a full TARDIS simulation using the base self.config
        modified by the user specified row_index.

        Parameters
        ----------
        row_index : int
            Row index in grid.

        Returns
        -------
        sim : tardis.simulation.base.Simulation
            Completed TARDIS simulation object.
        """
        tardis_config = self.grid_row_to_config(row_index)
        sim = tardis.run_tardis(tardis_config, **tardiskwargs)
        return sim

    def save_grid(self, filename):
        """
        Saves the parameter grid. Does not save the base
        self.config in any way.

        Parameters
        ----------
        filename : str
            File name to save grid.
        """
        self.grid.to_csv(filename, index=False)
        return

    @classmethod
    def from_axes(cls, configFile, axesdict):
        """
        Creates a grid from a set of axes. The axes are provided
        as a dictionary, where each key is a valid tardis config
        key, and the value is an iterable of values.

        Parameters
        ----------
        configFile : str
            Path to TARDIS yml file.
        axesdict : dict()
            Dictionary containing tardis config keys and the
            corresponding values to define a grid of tardis
            parameters.
        """
        axes = []
        dim = 1
        for key in axesdict:
            ax = axesdict[key]
            axes.append(ax)
            dim = dim * len(ax)
        axesmesh = np.meshgrid(*axes)
        tmp = np.dstack(axesmesh)
        gridpoints = tmp.reshape((dim, len(axes)), order="F")
        df = pd.DataFrame(data=gridpoints, columns=axesdict.keys())
        return cls(configFile=configFile, gridFrame=df)
