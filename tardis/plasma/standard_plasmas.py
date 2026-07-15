import pandas as pd

from tardis.plasma import BasePlasma
from tardis.plasma.base import PlasmaSolverSettings
from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.properties import (
    HeliumNumericalNLTE,
    IonNumberDensity,
    IonNumberDensityHeNLTE,
    LevelBoltzmannFactorNLTE,
    RadiationFieldCorrection,
    StimulatedEmissionFactor,
)
from tardis.plasma.properties.legacy_property_collections import (
    basic_inputs,
    basic_properties,
    dilute_lte_excitation_properties,
    helium_lte_properties,
    helium_nlte_properties,
    helium_numerical_nlte_properties,
    lte_excitation_properties,
    lte_ionization_properties,
    macro_atom_properties,
    nebular_ionization_properties,
    nlte_properties,
    non_nlte_properties,
)
from tardis.plasma.properties.level_population import LevelNumberDensity
from tardis.plasma.properties.nlte_rate_equation_solver import (
    NLTEPopulationSolverLU,
    NLTEPopulationSolverRoot,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.util.base import species_string_to_tuple


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
    if config.plasma.continuum_interaction.species:
        raise PlasmaConfigError(
            "Continuum interactions are supported only by the IIP workflow."
        )

    # Convert the nlte species list to a proper format.
    nlte_species = [
        species_string_to_tuple(species)
        for species in config.plasma.nlte.species
    ]

    # Convert the continuum interaction species list to a proper format.
    continuum_interaction_species = pd.MultiIndex.from_tuples(
        [], names=["atomic_number", "ion_number"]
    )

    atom_data.prepare_atom_data(
        simulation_state.abundance.index,
        line_interaction_type=config.plasma.line_interaction_type,
        continuum_interaction_species=continuum_interaction_species,
        nlte_species=nlte_species,
    )

    nlte_ionization_species = [
        species_string_to_tuple(species)
        for species in config.plasma.nlte_ionization_species
    ]
    nlte_excitation_species = [
        species_string_to_tuple(species)
        for species in config.plasma.nlte_excitation_species
    ]

    dilute_planckian_radiation_field = DilutePlanckianRadiationField(
        simulation_state.t_radiative, simulation_state.dilution_factor
    )
    kwargs = dict(
        dilute_planckian_radiation_field=dilute_planckian_radiation_field,
        abundance=simulation_state.abundance,
        number_density=simulation_state.calculate_elemental_number_density(
            atom_data.atom_data.mass
        ),
        atomic_data=atom_data,
        time_explosion=simulation_state.time_explosion,
        link_t_rad_t_electron=config.plasma.link_t_rad_t_electron,
        continuum_interaction_species=continuum_interaction_species,
        nlte_ionization_species=nlte_ionization_species,
        nlte_excitation_species=nlte_excitation_species,
    )

    plasma_modules = basic_inputs + basic_properties
    property_kwargs = {}

    ##### RADIATIVE RATES SETUP

    plasma_solver_settings = PlasmaSolverSettings(
        RADIATIVE_RATES_TYPE=config.plasma.radiative_rates_type
    )

    if (plasma_solver_settings.RADIATIVE_RATES_TYPE == "dilute-blackbody") or (
        plasma_solver_settings.RADIATIVE_RATES_TYPE == "detailed"
    ):
        kwargs["j_blues"] = pd.DataFrame(
            dilute_planckian_radiation_field.calculate_mean_intensity(
                atom_data.lines["nu"].values
            ),
            index=atom_data.lines.index,
        )

    elif plasma_solver_settings.RADIATIVE_RATES_TYPE == "blackbody":
        planckian_rad_field = (
            dilute_planckian_radiation_field.to_planckian_radiation_field()
        )
        kwargs["j_blues"] = pd.DataFrame(
            planckian_rad_field.calculate_mean_intensity(
                atom_data.lines["nu"].values
            ),
            index=atom_data.lines.index,
        )

    else:
        raise ValueError(
            f"radiative_rates_type type unknown - {plasma_solver_settings.RADIATIVE_RATES_TYPE}"
        )

    if config.plasma.excitation == "lte":
        plasma_modules += lte_excitation_properties
    elif config.plasma.excitation == "dilute-lte":
        plasma_modules += dilute_lte_excitation_properties

    if config.plasma.ionization == "lte":
        plasma_modules += lte_ionization_properties
    elif config.plasma.ionization == "nebular":
        plasma_modules += nebular_ionization_properties

    if nlte_species:
        plasma_modules += nlte_properties
        nlte_conf = config.plasma.nlte
        plasma_modules.append(LevelBoltzmannFactorNLTE.from_config(nlte_conf))
        property_kwargs[StimulatedEmissionFactor] = dict(
            nlte_species=nlte_species
        )
    else:
        plasma_modules += non_nlte_properties

    if config.plasma.line_interaction_type in ("downbranch", "macroatom"):
        plasma_modules += macro_atom_properties

    if "delta_treatment" in config.plasma:
        property_kwargs[RadiationFieldCorrection] = dict(
            delta_treatment=config.plasma.delta_treatment
        )

    if (
        config.plasma.helium_treatment == "recomb-nlte"
        or config.plasma.helium_treatment == "numerical-nlte"
    ) and (
        config.plasma.nlte_ionization_species
        or config.plasma.nlte_excitation_species
    ):
        # Prevent the user from using helium NLTE treatment with
        # NLTE ionization and excitation treatment. This is because
        # the helium_nlte_properties could overwrite the NLTE ionization
        # and excitation ion number and electron densities.
        # helium_numerical_nlte_properties is also included here because
        # it is currently in the same if else block, and thus may block
        # the addition of the components from the else block.
        raise PlasmaConfigError(
            "Helium NLTE treatment is incompatible with the NLTE eonization and excitation treatment."
        )

    # TODO: Disentangle these if else block such that compatible components
    # can be added independently.
    if config.plasma.helium_treatment == "recomb-nlte":
        plasma_modules += helium_nlte_properties
    elif config.plasma.helium_treatment == "numerical-nlte":
        plasma_modules += helium_numerical_nlte_properties
        # TODO: See issue #633
        if config.plasma.heating_rate_data_file in ["none", None]:
            raise PlasmaConfigError("Heating rate data file not specified")
        else:
            property_kwargs[HeliumNumericalNLTE] = dict(
                heating_rate_data_file=config.plasma.heating_rate_data_file
            )
    else:
        # If nlte ionization species are present, we don't want to add the
        # IonNumberDensity from helium_lte_properties, since we want
        # to use the IonNumberDensity provided by the NLTE solver.
        if (
            config.plasma.nlte_ionization_species
            or config.plasma.nlte_excitation_species
        ):
            plasma_modules += [LevelNumberDensity]
        else:
            plasma_modules += helium_lte_properties

    if simulation_state._electron_densities is not None:
        electron_densities = pd.Series(
            simulation_state._electron_densities.cgs.value
        )
        if config.plasma.helium_treatment == "numerical-nlte":
            property_kwargs[IonNumberDensityHeNLTE] = dict(
                electron_densities=electron_densities
            )
        elif (
            config.plasma.nlte_ionization_species
            or config.plasma.nlte_excitation_species
        ) and config.plasma.nlte_solver == "root":
            property_kwargs[NLTEPopulationSolverRoot] = dict(
                electron_densities=electron_densities
            )
        elif (
            config.plasma.nlte_ionization_species
            or config.plasma.nlte_excitation_species
        ) and config.plasma.nlte_solver == "lu":
            property_kwargs[NLTEPopulationSolverLU] = dict(
                electron_densities=electron_densities
            )
        else:
            property_kwargs[IonNumberDensity] = dict(
                electron_densities=electron_densities
            )

    kwargs["helium_treatment"] = config.plasma.helium_treatment

    plasma = BasePlasma(
        plasma_properties=plasma_modules,
        property_kwargs=property_kwargs,
        plasma_solver_settings=plasma_solver_settings,
        **kwargs,
    )

    return plasma
