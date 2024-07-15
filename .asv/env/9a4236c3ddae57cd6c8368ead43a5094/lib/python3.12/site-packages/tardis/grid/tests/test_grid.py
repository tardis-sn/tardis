from pathlib import Path

import numpy as np
import pandas as pd

import tardis

from pathlib import Path
import tardis.grid as grid


DATA_PATH = Path(tardis.__path__[0]) / "grid" / "tests" / "data"


def test_grid(atomic_dataset):
    """Tests the basic functionality of the TARDIS grid module."""
    dfpath = DATA_PATH / "example_grid.txt"
    ymlpath = DATA_PATH / "example.yml"
    axesdict = {
        "model.structure.velocity.start": np.arange(10000, 15000, 1000),
        "model.abundances.He": np.arange(0, 1, 0.1),
        "model.abundances.H": np.arange(0, 1, 0.25),
    }

    df = pd.read_csv(dfpath)
    g = grid.tardisGrid(configFile=ymlpath, gridFrame=df)
    g2 = grid.tardisGrid.from_axes(configFile=ymlpath, axesdict=axesdict)

    # Check that grid attribute has the right shape
    assert g.grid.shape == df.shape
    ax_len = 1
    for key in axesdict:
        ax_len *= len(axesdict[key])
    assert g2.grid.shape[0] == ax_len

    newconf = g.grid_row_to_config(row_index=0)
    # Verify that grid_row_to_config modifies the base config attribute
    assert g.config != newconf

    # Verify that a model can be returned.
    simulation_state = g.grid_row_to_simulation_state(
        row_index=0, atomic_data=atomic_dataset
    )
    assert (
        simulation_state.velocity[0].to("km/s").value
        == df.iloc[0]["model.structure.velocity.start"]
    )
