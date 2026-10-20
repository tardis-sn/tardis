import logging

import numpy as np

# Energy: keV, exported as eV for SF solver
# distance: cm
# mass: g
# time: s
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def calculate_ejecta_velocity_volume(model):
    outer_velocities = model.v_outer.to("cm/s").value
    inner_velocities = model.v_inner.to("cm/s").value
    ejecta_velocity_volume = (
        4 * np.pi / 3 * (outer_velocities**3.0 - inner_velocities**3.0)
    )

    return ejecta_velocity_volume


def calculate_average_power_per_mass(energy_per_mass, time_delta):
    # Time averaged energy per mass for constant packet count
    average_power_per_mass = energy_per_mass / (time_delta)

    return average_power_per_mass


def iron_group_fraction_per_shell(model):
    # Taking iron group to be elements 21-30
    # Used as part of the approximations for photoabsorption and pair creation
    # Dependent on atomic data
    return model.abundance.loc[(21):(30)].sum(axis=0)
