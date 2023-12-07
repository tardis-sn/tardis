from tardis.io.model.model_reader import simulation_state_to_dict


import h5py


def store_model_to_hdf(model, fname):
    """
    Stores data from SimulationState object into a hdf file.

    Parameters
    ----------
    model : tardis.model.SimulationState
    filename : str
    """
    with h5py.File(fname, "a") as f:
        model_group = f.require_group("model")
        model_group.clear()

        model_dict = simulation_state_to_dict(model)

        for key, value in model_dict.items():
            if key.endswith("_cgs"):
                model_group.create_dataset(key, data=value[0])
                model_group.create_dataset(key + "_unit", data=value[1])
            else:
                model_group.create_dataset(key, data=value)
