from tardis.io.model.model_reader import simulation_state_to_dict


import h5py


def store_simulation_state_to_hdf(simulation_state, fname):
    """
    Stores data from SimulationState object into a hdf file.

    Parameters
    ----------
    model : tardis.model.SimulationState
    filename : str
    """
    with h5py.File(fname, "a") as f:
        simulation_state_group = f.require_group("simulation_state")
        simulation_state_group.clear()

        simulation_state_dict = simulation_state_to_dict(simulation_state)

        for key, value in simulation_state_dict.items():
            if key.endswith("_cgs"):
                simulation_state_group.create_dataset(key, data=value[0])
                simulation_state_group.create_dataset(
                    key + "_unit", data=value[1]
                )
            else:
                simulation_state_group.create_dataset(key, data=value)
