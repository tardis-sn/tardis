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
from tardis.transport.montecarlo.estimators.continuum_radfield_properties import (
    DiluteBlackBodyContinuumPropertiesSolver,
)
from tardis.util.base import species_string_to_tuple

logger = logging.getLogger(__name__)

def map_species_from_string(species):
    return [species_string_to_tuple(spec) for spec in species]


class PlasmaSolverFactory:

    ## Analytical Approximations
    excitation_analytical_approximation: str = "lte"
    ionization_analytical_approximation: str = "lte"
    nebular_ionization_delta_treatment: tuple  # species to use for the delta_treatment in nebular ionization ML93

    link_t_rad_t_electron: float = 1.0

    radiative_rates_type: str = "dilute-blackbody"

    delta_treatment = None

    ## Statistical Balance Solver
    legacy_nlte_species: list = []

    nlte_excitation_species: list = []
    nlte_ionization_species: list = []
    nlte_solver: str = "lu"

    ## Helium Treatment options
    helium_treatment: str = "none"
    heating_rate_data_file: str = "none"

    ## Continuum Interaction
    continuum_interaction_species: pd.MultiIndex

    ## Opacities
    line_interaction_type: str = "scatter"

    ## Assembly properties
    plasma_modules: list = []
    kwargs: dict = {}
    property_kwargs: dict = {}

    def __init__(self, config, atom_data, atomic_numbers) -> None:

        self.set_nlte_species_from_string(config.plasma.nlte.species)
        self.set_continuum_interaction_species_from_string(
            config.plasma.continuum_interaction.species
        )
        self.line_interaction_type = config.plasma.line_interaction_type

        self.atom_data = atom_data
        self.atom_data.prepare_atom_data(
            atomic_numbers,
            line_interaction_type=config.plasma.line_interaction_type,
            continuum_interaction_species=self.continuum_interaction_species,
            nlte_species=self.legacy_nlte_species,
        )

        #### THIS IS VERY BAD BUT FOR NOW IS CHICKEN/EGG
        # Check if continuum interaction species are in selected_atoms
        continuum_atoms = self.continuum_interaction_species.get_level_values(
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
        ##### ----------------------------

        self.plasma_modules = basic_inputs + basic_properties
        self.link_t_rad_t_electron = config.plasma.link_t_rad_t_electron

        self.setup_analytical_approximations(config.plasma)
        if "delta_treatment" in config.plasma:
            self.property_kwargs[RadiationFieldCorrection] = dict(
                delta_treatment=config.plasma.delta_treatment
            )

        self.setup_legacy_nlte(config.plasma.nlte)
        if self.line_interaction_type in ("downbranch", "macroatom") and (
            len(self.continuum_interaction_species) == 0
        ):
            self.setup_legacy_macro_atom()
        self.delta_treatment = config.plasma.get("delta_treatment", None)
        if self.delta_treatment is not None:
            self.property_kwargs[RadiationFieldCorrection] = dict(
                delta_treatment=config.plasma.delta_treatment
            )

        self.helium_treatment = config.plasma.helium_treatment
        self.heating_rate_data_file = config.plasma.heating_rate_data_file
        self.setup_helium_treatment()

        self.nlte_solver = config.plasma.nlte_solver
        self.nlte_ionization_species = config.plasma.nlte_ionization_species
        self.nlte_excitation_species = config.plasma.nlte_excitation_species
        if len(config.plasma.continuum_interaction.species) > 0:
            self.setup_continuum_interactions(
                config.plasma.continuum_interaction
            )
        self.radiative_rates_type = config.plasma.radiative_rates_type
    def setup_helium_treatment(self):
        """
        Set up the helium treatment for the plasma assembly.

        Parameters
        ----------
        helium_treatment : str
            The type of helium treatment to be used. Possible values are:
            - "recomb-nlte": Use recombination NLTE treatment for helium.
            - "numerical-nlte": Use numerical NLTE treatment for helium.

        heating_rate_data_file : str or None
            The path to the heating rate data file. Required when using
            "numerical-nlte" helium treatment.

        Raises
        ------
        PlasmaConfigError
            If the helium NLTE treatment is incompatible with the NLTE ionization
            and excitation treatment.

            If the heating rate data file is not specified when using
            "numerical-nlte" helium treatment.
        """
        if (
            self.helium_treatment == "recomb-nlte"
            or self.helium_treatment == "numerical-nlte"
        ) and (
            len(self.nlte_ionization_species + self.nlte_excitation_species) > 0
        ):
            # Prevent the user from using helium NLTE treatment with
            # NLTE ionization and excitation treatment. This is because
            # the helium_nlte_properties could overwrite the NLTE ionization
            # and excitation ion number and electron densities.
            # helium_numerical_nlte_properties is also included here because
            # it is currently in the same if else block, and thus may block
            # the addition of the components from the else block.
            raise PlasmaConfigError(
                "Helium NLTE treatment is incompatible with the NLTE ionization and excitation treatment."
            )

        # TODO: Disentangle these if else block such that compatible components
        # can be added independently.
        if self.helium_treatment == "recomb-nlte":
            self.plasma_modules += helium_nlte_properties
        elif self.helium_treatment == "numerical-nlte":
            self.plasma_modules += helium_numerical_nlte_properties
            if heating_rate_data_file in ["none", None]:
                raise PlasmaConfigError("Heating rate data file not specified")
            self.property_kwargs[HeliumNumericalNLTE] = dict(
                heating_rate_data_file=heating_rate_data_file
            )
        else:
            # If nlte ionization species are present, we don't want to add the
            # IonNumberDensity from helium_lte_properties, since we want
            # to use the IonNumberDensity provided by the NLTE solver.
            if (
                len(self.nlte_ionization_species + self.nlte_excitation_species)
                > 0
            ):
                self.plasma_modules.append(LevelNumberDensity)
            else:
                self.plasma_modules += helium_lte_properties

    def setup_legacy_macro_atom(self):
        self.plasma_modules += macro_atom_properties

    def setup_legacy_nlte(self, nlte_config):
        """
        Set up the non-LTE (NLTE) properties for the legacy species.

        Parameters:
        -----------
        nlte_config : dict
            A dictionary containing the NLTE configuration.

        Notes:
        ------
        This method adds the NLTE properties for the legacy species to the plasma modules.
        If there are no legacy NLTE species, it adds the non-NLTE properties instead.
        """
        if len(self.legacy_nlte_species) > 0:
            self.plasma_modules += nlte_properties
            self.plasma_modules.append(
                LevelBoltzmannFactorNLTE.from_config(nlte_config)
            )
            self.property_kwargs[StimulatedEmissionFactor] = dict(
                nlte_species=self.legacy_nlte_species
            )
        else:
            self.plasma_modules += non_nlte_properties

    def setup_analytical_approximations(self, plasma_config):
        """
        Setup the analytical approximations for excitation and ionization.

        Returns
        -------
        None
        """
        self.excitation_analytical_approximation = plasma_config.excitation
        self.ionization_analytical_approximation = plasma_config.ionization

        if self.excitation_analytical_approximation == "lte":
            self.plasma_modules += lte_excitation_properties
        elif self.excitation_analytical_approximation == "dilute-lte":
            self.plasma_modules += dilute_lte_excitation_properties
        else:
            raise PlasmaConfigError(
                f'Invalid excitation analytical approximation. Configured as {self.excitation_analytical_approximation} but needs to be either "lte" or "dilute-lte"'
            )

        if self.ionization_analytical_approximation == "lte":
            self.plasma_modules += lte_ionization_properties
        elif self.ionization_analytical_approximation == "nebular":
            self.plasma_modules += nebular_ionization_properties
        else:
            raise PlasmaConfigError(
                f'Invalid excitation analytical approximation. Configured as {self.ionization_analytical_approximation} but needs to be either "lte" or "nebular"'
            )

    def initialize_j_blues(self, dilute_planckian_radiation_field, lines_df):
        if (self.radiative_rates_type == "dilute-blackbody") or (
            self.radiative_rates_type == "detailed"
        ):
            j_blues = pd.DataFrame(
                dilute_planckian_radiation_field.calculate_mean_intensity(
                    lines_df.nu.values
                ),
                index=lines_df.index,
            )

        elif self.radiative_rates_type == "blackbody":
            planckian_rad_field = (
                dilute_planckian_radiation_field.to_planckian_radiation_field()
            )
            j_blues = pd.DataFrame(
                planckian_rad_field.calculate_mean_intensity(
                    lines_df.nu.values
                ),
                index=lines_df.index,
            )

        else:
            raise ValueError(
                f"radiative_rates_type type unknown - {self.radiative_rates_type}"
            )

        return j_blues

    def set_continuum_interaction_species_from_string(
        self, continuum_interaction_species
    ):
        """
        Set the continuum interaction species from a list of species strings.

        Parameters
        ----------
        continuum_interaction_species : list of str
            List of species strings representing the continuum interaction species.

        Returns
        -------
        None
        """
        continuum_interaction_species = [
            species_string_to_tuple(species)
            for species in continuum_interaction_species
        ]

        self.continuum_interaction_species = pd.MultiIndex.from_tuples(
            continuum_interaction_species, names=["atomic_number", "ion_number"]
        )

    def set_nlte_species_from_string(self, nlte_species):
        """
        Sets the non-LTE species from a string representation.

        Parameters
        ----------
        nlte_species : str
            A string representation of the non-LTE species.

        Returns
        -------
        None
            This method does not return anything.
        """
        self.legacy_nlte_species = map_species_from_string(nlte_species)

    def setup_continuum_interactions(self, config_continuum_interaction):

        if self.line_interaction_type != "macroatom":
            raise PlasmaConfigError(
                "Continuum interactions require line_interaction_type "
                f"macroatom (instead of {self.line_interaction_type})."
            )

        self.plasma_modules += continuum_interaction_properties
        self.plasma_modules += continuum_interaction_inputs

        if config_continuum_interaction.enable_adiabatic_cooling:
            self.plasma_modules += adiabatic_cooling_properties

        if config_continuum_interaction.enable_two_photon_decay:
            self.plasma_modules += two_photon_properties

        transition_probabilities_outputs = [
            plasma_property.transition_probabilities_outputs
            for plasma_property in self.plasma_modules
            if issubclass(plasma_property, TransitionProbabilitiesProperty)
        ]
        transition_probabilities_outputs = [
            item
            for sublist in transition_probabilities_outputs
            for item in sublist
        ]

        self.property_kwargs[MarkovChainTransProbsCollector] = {
            "inputs": transition_probabilities_outputs
        }
        if len(self.nlte_ionization_species + self.nlte_excitation_species) > 0:
            if self.nlte_ionization_species:
                nlte_ionization_species = self.nlte_ionization_species
                for species in nlte_ionization_species:
                    if species not in config_continuum_interaction.species:
                        raise PlasmaConfigError(
                            f"NLTE ionization species {species} not in continuum species."
                        )
            if self.nlte_excitation_species:
                nlte_excitation_species = self.nlte_excitation_species
                for species in nlte_excitation_species:
                    if species not in config_continuum_interaction.species:
                        raise PlasmaConfigError(
                            f"NLTE excitation species {species} not in continuum species."
                        )
            self.property_kwargs[NLTEIndexHelper] = {
                "nlte_ionization_species": self.nlte_ionization_species,
                "nlte_excitation_species": self.nlte_excitation_species,
            }
            if self.nlte_solver == "lu":
                self.plasma_modules += nlte_lu_solver_properties
                logger.warning(
                    "LU solver will be inaccurate for NLTE excitation, proceed with caution."
                )
            elif self.nlte_solver == "root":
                self.plasma_modules += nlte_root_solver_properties
            else:
                raise PlasmaConfigError(
                    f"NLTE solver type unknown - {self.nlte_solver}"
                )

    def setup_electron_densities(self, electron_densities):
        if self.helium_treatment == "numerical-nlte":
            self.property_kwargs[IonNumberDensityHeNLTE] = dict(
                electron_densities=electron_densities
            )
        elif (
            len(self.nlte_ionization_species + self.nlte_excitation_species) > 0
        ) and self.nlte_solver == "root":
            self.property_kwargs[NLTEPopulationSolverRoot] = dict(
                electron_densities=electron_densities
            )
        elif (
            len(self.nlte_ionization_species + self.nlte_excitation_species) > 0
        ) and self.nlte_solver == "lu":
            self.property_kwargs[NLTEPopulationSolverLU] = dict(
                electron_densities=electron_densities
            )
        else:
            self.property_kwargs[IonNumberDensity] = dict(
                electron_densities=electron_densities
            )

    def initialize_continuum_properties(self, dilute_planckian_radiation_field):
        """
        Initialize the continuum properties of the plasma.

        Parameters
        ----------
        dilute_planckian_radiation_field : DilutePlanckianRadiationField
            The dilute Planckian radiation field.

        Returns
        -------
        initial_continuum_properties : `~tardis.plasma.properties.ContinuumProperties`
            The initial continuum properties of the plasma.
        """
        t_electrons = dilute_planckian_radiation_field.temperature.to(u.K).value

        initial_continuum_solver = DiluteBlackBodyContinuumPropertiesSolver(
            self.atom_data
        )
        initial_continuum_properties = initial_continuum_solver.solve(
            dilute_planckian_radiation_field, t_electrons
        )
        return initial_continuum_properties

    def assemble(self, dilute_planckian_radiation_field, simulation_state):
        j_blues = self.initialize_j_blues(
            dilute_planckian_radiation_field, self.atom_data.lines
        )

        kwargs = dict(
            time_explosion=simulation_state.time_explosion,
            dilute_planckian_radiation_field=dilute_planckian_radiation_field,
            abundance=simulation_state.abundance,
            number_density=simulation_state.elemental_number_density,
            link_t_rad_t_electron=self.link_t_rad_t_electron,
            atomic_data=self.atom_data,
            j_blues=j_blues,
            continuum_interaction_species=self.continuum_interaction_species,
            nlte_ionization_species=self.nlte_ionization_species,
            nlte_excitation_species=self.nlte_excitation_species,
        )

        if len(self.continuum_interaction_species) > 0:
            initial_continuum_properties = self.initialize_continuum_properties(
                dilute_planckian_radiation_field
            )
            kwargs.update(
                gamma=initial_continuum_properties.photo_ionization_rate_coefficient,
                bf_heating_coeff_estimator=None,
                stim_recomb_cooling_coeff_estimator=None,
                alpha_stim_factor=initial_continuum_properties.stimulated_recombination_rate_factor,
            )

        if simulation_state._electron_densities is not None:
            electron_densities = pd.Series(
                simulation_state._electron_densities.cgs.value
            )
            self.setup_electron_densities(electron_densities)
        kwargs["helium_treatment"] = self.helium_treatment
        return BasePlasma(
            plasma_properties=self.plasma_modules,
            property_kwargs=self.property_kwargs,
            **kwargs,
        )
