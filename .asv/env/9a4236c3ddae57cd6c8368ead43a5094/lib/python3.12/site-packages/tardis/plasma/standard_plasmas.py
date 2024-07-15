import logging

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.plasma import BasePlasma
from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.properties import (
    HeliumNumericalNLTE,
    IonNumberDensity,
    IonNumberDensityHeNLTE,
    JBluesBlackBody,
    JBluesDetailed,
    JBluesDiluteBlackBody,
    LevelBoltzmannFactorNLTE,
    MarkovChainTransProbsCollector,
    RadiationFieldCorrection,
    StimulatedEmissionFactor,
)
from tardis.plasma.properties.base import TransitionProbabilitiesProperty
from tardis.plasma.properties.level_population import LevelNumberDensity
from tardis.plasma.properties.nlte_rate_equation_solver import (
    NLTEPopulationSolverLU,
    NLTEPopulationSolverRoot,
)
from tardis.plasma.properties.property_collections import (
    adiabatic_cooling_properties,
    basic_inputs,
    basic_properties,
    continuum_interaction_inputs,
    continuum_interaction_properties,
    detailed_j_blues_inputs,
    detailed_j_blues_properties,
    dilute_lte_excitation_properties,
    helium_lte_properties,
    helium_nlte_properties,
    helium_numerical_nlte_properties,
    lte_excitation_properties,
    lte_ionization_properties,
    macro_atom_properties,
    nebular_ionization_properties,
    nlte_lu_solver_properties,
    nlte_properties,
    nlte_root_solver_properties,
    non_nlte_properties,
    two_photon_properties,
)
from tardis.plasma.properties.rate_matrix_index import NLTEIndexHelper
from tardis.util.base import species_string_to_tuple

logger = logging.getLogger(__name__)


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
    # Convert the nlte species list to a proper format.
    nlte_species = [
        species_string_to_tuple(s) for s in config.plasma.nlte.species
    ]

    # Convert the continuum interaction species list to a proper format.
    continuum_interaction_species = [
        species_string_to_tuple(s)
        for s in config.plasma.continuum_interaction.species
    ]
    continuum_interaction_species = pd.MultiIndex.from_tuples(
        continuum_interaction_species, names=["atomic_number", "ion_number"]
    )

    atom_data.prepare_atom_data(
        simulation_state.abundance.index,
        line_interaction_type=config.plasma.line_interaction_type,
        continuum_interaction_species=continuum_interaction_species,
        nlte_species=nlte_species,
    )

    # Check if continuum interaction species are in selected_atoms
    continuum_atoms = continuum_interaction_species.get_level_values(
        "atomic_number"
    )
    continuum_atoms_in_selected_atoms = np.all(
        continuum_atoms.isin(atom_data.selected_atomic_numbers)
    )
    if not continuum_atoms_in_selected_atoms:
        raise PlasmaConfigError(
            "Not all continuum interaction species "
            "belong to atoms that have been specified "
            "in the configuration."
        )

    nlte_ionization_species = [
        species_string_to_tuple(s)
        for s in config.plasma.nlte_ionization_species
    ]
    nlte_excitation_species = [
        species_string_to_tuple(s)
        for s in config.plasma.nlte_excitation_species
    ]

    kwargs = dict(
        t_rad=simulation_state.t_radiative,
        abundance=simulation_state.abundance,
        number_density=simulation_state.elemental_number_density,
        atomic_data=atom_data,
        time_explosion=simulation_state.time_explosion,
        w=simulation_state.dilution_factor,
        link_t_rad_t_electron=config.plasma.link_t_rad_t_electron,
        continuum_interaction_species=continuum_interaction_species,
        nlte_ionization_species=nlte_ionization_species,
        nlte_excitation_species=nlte_excitation_species,
    )

    plasma_modules = basic_inputs + basic_properties
    property_kwargs = {}
    if len(config.plasma.continuum_interaction.species) > 0:
        line_interaction_type = config.plasma.line_interaction_type
        if line_interaction_type != "macroatom":
            raise PlasmaConfigError(
                "Continuum interactions require line_interaction_type "
                f"macroatom (instead of {line_interaction_type})."
            )

        plasma_modules += continuum_interaction_properties
        plasma_modules += continuum_interaction_inputs

        if config.plasma.continuum_interaction.enable_adiabatic_cooling:
            plasma_modules += adiabatic_cooling_properties

        if config.plasma.continuum_interaction.enable_two_photon_decay:
            plasma_modules += two_photon_properties

        transition_probabilities_outputs = [
            plasma_property.transition_probabilities_outputs
            for plasma_property in plasma_modules
            if issubclass(plasma_property, TransitionProbabilitiesProperty)
        ]
        transition_probabilities_outputs = [
            item
            for sublist in transition_probabilities_outputs
            for item in sublist
        ]

        property_kwargs[MarkovChainTransProbsCollector] = {
            "inputs": transition_probabilities_outputs
        }
        if (
            config.plasma.nlte_ionization_species
            or config.plasma.nlte_excitation_species
        ):
            if config.plasma.nlte_ionization_species:
                nlte_ionization_species = config.plasma.nlte_ionization_species
                for species in nlte_ionization_species:
                    if (
                        species
                        not in config.plasma.continuum_interaction.species
                    ):
                        raise PlasmaConfigError(
                            f"NLTE ionization species {species} not in continuum species."
                        )
            if config.plasma.nlte_excitation_species:
                nlte_excitation_species = config.plasma.nlte_excitation_species
                for species in nlte_excitation_species:
                    if (
                        species
                        not in config.plasma.continuum_interaction.species
                    ):
                        raise PlasmaConfigError(
                            f"NLTE excitation species {species} not in continuum species."
                        )
            property_kwargs[NLTEIndexHelper] = {
                "nlte_ionization_species": config.plasma.nlte_ionization_species,
                "nlte_excitation_species": config.plasma.nlte_excitation_species,
            }
            if config.plasma.nlte_solver == "lu":
                plasma_modules += nlte_lu_solver_properties
                logger.warning(
                    "LU solver will be inaccurate for NLTE excitation, proceed with caution."
                )
            elif config.plasma.nlte_solver == "root":
                plasma_modules += nlte_root_solver_properties
            else:
                raise PlasmaConfigError(
                    f"NLTE solver type unknown - {config.plasma.nlte_solver}"
                )

        kwargs.update(
            gamma_estimator=None,
            bf_heating_coeff_estimator=None,
            stim_recomb_cooling_coeff_estimator=None,
            alpha_stim_estimator=None,
            volume=simulation_state.volume,
            r_inner=simulation_state.r_inner.to(u.cm),
            t_inner=simulation_state.t_inner,
        )
    if config.plasma.radiative_rates_type == "blackbody":
        plasma_modules.append(JBluesBlackBody)
    elif config.plasma.radiative_rates_type == "dilute-blackbody":
        plasma_modules.append(JBluesDiluteBlackBody)
    elif config.plasma.radiative_rates_type == "detailed":
        plasma_modules += detailed_j_blues_properties + detailed_j_blues_inputs
        kwargs.update(
            r_inner=simulation_state.r_inner.to(u.cm),
            t_inner=simulation_state.t_inner,
            volume=simulation_state.volume,
            j_blue_estimator=None,
        )
        property_kwargs[JBluesDetailed] = {"w_epsilon": config.plasma.w_epsilon}
    else:
        raise ValueError(
            f"radiative_rates_type type unknown - {config.plasma.radiative_rates_type}"
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
        if not config.plasma.continuum_interaction.species:
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
        **kwargs,
    )

    return plasma
