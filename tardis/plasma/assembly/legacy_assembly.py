from tardis.plasma.assembly.base import PlasmaSolverFactory
from tardis.plasma.radiation_field import DilutePlanckianRadiationField


def assemble_plasma(config, simulation_state, atom_data=None):
    """
    Create a BasePlasma instance from a Configuration object
    and a SimulationState.

    Parameters
    ----------
    config : io.config_reader.Configuration
    simulation_state : model.SimulationState
    atom_data : atomic.AtomData
        If None, an attempt will be made to read the atomic data
        from config.

    Returns
    -------
    : plasma.BasePlasma

    """
    atomic_numbers = simulation_state.abundance.index
    plasma_solver_factory = PlasmaSolverFactory(
        atom_data,
        atomic_numbers,
        config,
    )
    plasma_solver_factory.setup_factory(config)
    dilute_planckian_radiation_field = DilutePlanckianRadiationField(
        simulation_state.t_radiative, simulation_state.dilution_factor
    )

    return plasma_solver_factory.assemble(
        simulation_state.elemental_number_density,
        dilute_planckian_radiation_field,
        simulation_state.time_explosion,
        simulation_state._electron_densities,
    )
