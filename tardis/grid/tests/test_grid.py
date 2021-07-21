import pandas as pd
import numpy as np
import tardis
import os
import pytest
import tardis.grid as grid


DATA_PATH = os.path.join(tardis.__path__[0], "grid", "tests", "data")


def test_grid():
    """Tests the basic functionality of the TARDIS grid module."""
    dfpath = os.path.join(DATA_PATH, "example_grid.txt")
    ymlpath = os.path.join(DATA_PATH, "example.yml")
    axesdict = {
        "model.structure.velocity.start": np.arange(10000, 15000, 1000),
        "model.abundances.He": np.arange(0, 1, 0.1),
        "model.abundances.H": np.arange(0, 1, 0.25),
    }

    df = pd.read_csv(os.path.join(dfpath))
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
    model = g.grid_row_to_model(row_index=0)
    assert (
        model.velocity[0].to("km/s").value
        == df.iloc[0]["model.structure.velocity.start"]
    )
