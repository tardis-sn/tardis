import importlib
import logging

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.plasma import BasePlasma
from tardis.plasma.base import PlasmaSolverSettings
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
from tardis.plasma.properties.rate_matrix_index import NLTEIndexHelper
from tardis.transport.montecarlo.estimators.continuum_radfield_properties import (
    DiluteBlackBodyContinuumPropertiesSolver,
)
from tardis.util.base import species_string_to_tuple

logger = logging.getLogger(__name__)


def map_species_from_string(species):
    return [species_string_to_tuple(spec) for spec in species]


def convert_species_to_multi_index(species_strs):
    return pd.MultiIndex.from_tuples(
        map_species_from_string(species_strs),
        names=["atomic_number", "ion_number"],
    )


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
    continuum_interaction_species: list = []
    enable_adiabatic_cooling: bool = False
    enable_two_photon_decay: bool = False

    ## Opacities
    line_interaction_type: str = "scatter"

    ## Assembly properties
    plasma_modules: list = []
    kwargs: dict = {}
    property_kwargs: dict = {}
    plasma_collection = importlib.import_module(
        "tardis.plasma.properties.legacy_property_collections"
    )

    def __init__(
        self,
        atom_data,
        config=None,
    ) -> None:
        if config is not None:
            self.parse_plasma_config(config.plasma)
        self.atom_data = atom_data

    @property
    def continuum_interaction_species_multi_index(self):
        return convert_species_to_multi_index(
            self.continuum_interaction_species
        )

    def parse_plasma_config(self, plasma_config):
        """
        Parse the plasma configuration.

        Parameters
        ----------
        plasma_config : PlasmaConfig
            The plasma configuration object containing the plasma parameters.

        Returns
        -------
        None
        """
        self.continuum_interaction_species = (
            plasma_config.continuum_interaction.species
        )
        self.set_nlte_species_from_string(plasma_config.nlte.species)
        self.line_interaction_type = plasma_config.line_interaction_type
        self.link_t_rad_t_electron = plasma_config.link_t_rad_t_electron

        self.excitation_analytical_approximation = plasma_config.excitation
        self.ionization_analytical_approximation = plasma_config.ionization
        self.delta_treatment = plasma_config.get("delta_treatment", None)

        self.helium_treatment = plasma_config.helium_treatment
        self.heating_rate_data_file = plasma_config.heating_rate_data_file

        self.nlte_ionization_species = plasma_config.nlte_ionization_species
        self.nlte_excitation_species = plasma_config.nlte_excitation_species

        self.nlte_solver = plasma_config.nlte_solver

        self.radiative_rates_type = plasma_config.radiative_rates_type

        self.enable_adiabatic_cooling = (
            plasma_config.continuum_interaction.enable_adiabatic_cooling
        )
        self.enable_two_photon_decay = (
            plasma_config.continuum_interaction.enable_two_photon_decay
        )

    def prepare_factory(
        self, selected_atomic_numbers, property_collections, config=None
    ):
        """
        Set up the plasma factory.

        Parameters
        ----------
        selected_atomic_numbers : list of int
            Selected atomic numbers in the simulation.
        property_collections : str
            The property collection module to be used in the plasma assembly.
        config : object, optional
            Configuration object containing plasma settings (default: None).
        """
        self.plasma_collection = importlib.import_module(property_collections)

        self.atom_data.prepare_atom_data(
            selected_atomic_numbers,
            line_interaction_type=self.line_interaction_type,
            continuum_interaction_species=self.continuum_interaction_species_multi_index,
            nlte_species=self.legacy_nlte_species,
        )

        self.check_continuum_interaction_species()

        self.plasma_modules = (
            self.plasma_collection.basic_inputs
            + self.plasma_collection.basic_properties
        )

        self.setup_analytical_approximations()
        self.property_kwargs[RadiationFieldCorrection] = dict(
            delta_treatment=self.delta_treatment
        )
        if (config is not None) and len(self.legacy_nlte_species) > 0:
            self.setup_legacy_nlte(config.plasma.nlte)
        else:
            self.plasma_modules += self.plasma_collection.non_nlte_properties

        if self.line_interaction_type in ("downbranch", "macroatom") and (
            len(self.continuum_interaction_species) == 0
        ):
            self.plasma_modules += self.plasma_collection.macro_atom_properties

        self.setup_helium_treatment()


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
        
        # nlte_exists determines if NLTE ionization and excitation ion 
        # number and electron densities exists so that it can prevent 
        # user from using Helium NLTE treatments 
        nlte_exists = len(self.nlte_ionization_species + self.nlte_excitation_species) > 0

        if self.helium_treatment == "recomb-nlte":
            self.set_recomb_nlte_helium_treatment(nlte_exists)
        elif self.helium_treatment == "numerical-nlte":
            self.set_numerical_nlte_helium_treatment(nlte_exists)
        else:
            self.set_default_helium_treatment(nlte_exists)


    def set_recomb_nlte_helium_treatment(self, default_nlte):
        if default_nlte:
            raise PlasmaConfigError(
                "Helium NLTE treatment is incompatible with the NLTE ionization and excitation treatment."
            )

        self.plasma_modules += self.plasma_collection.helium_nlte_properties


    def set_numerical_nlte_helium_treatment(self, default_nlte):
        if default_nlte:
            raise PlasmaConfigError(
                "Helium NLTE treatment is incompatible with the NLTE ionization and excitation treatment."
            )
        
        self.plasma_modules += (
                self.plasma_collection.helium_numerical_nlte_properties
            )
        
        if self.heating_rate_data_file in ["none", None]:
                raise PlasmaConfigError("Heating rate data file not specified")
        self.property_kwargs[HeliumNumericalNLTE] = dict(
                heating_rate_data_file=self.heating_rate_data_file
            )
        

    def set_default_helium_treatment(self, default_nlte):
        # If nlte ionization species are present, we don't want to add the
        # IonNumberDensity from helium_lte_properties, since we want
        # to use the IonNumberDensity provided by the NLTE solver.
        self.plasma_modules += [LevelNumberDensity]
        if not default_nlte:
            self.plasma_modules += [IonNumberDensity]
            

    def setup_legacy_nlte(self, nlte_config):
        """
        Set up the non-LTE (NLTE) properties for the legacy species.

        Parameters
        ----------
        nlte_config : dict
            A dictionary containing the NLTE configuration.
        """
        self.plasma_modules += self.plasma_collection.nlte_properties
        self.plasma_modules.append(
            LevelBoltzmannFactorNLTE.from_config(nlte_config)
        )
        self.property_kwargs[StimulatedEmissionFactor] = dict(
            nlte_species=self.legacy_nlte_species
        )

    def setup_analytical_approximations(self):
        """
        Setup the analytical approximations for excitation and ionization.

        Returns
        -------
        None
        """
        if self.excitation_analytical_approximation == "lte":
            self.plasma_modules += (
                self.plasma_collection.lte_excitation_properties
            )
        elif self.excitation_analytical_approximation == "dilute-lte":
            self.plasma_modules += (
                self.plasma_collection.dilute_lte_excitation_properties
            )
        else:
            raise PlasmaConfigError(
                f'Invalid excitation analytical approximation. Configured as {self.excitation_analytical_approximation} but needs to be either "lte" or "dilute-lte"'
            )

        if self.ionization_analytical_approximation == "lte":
            self.plasma_modules += (
                self.plasma_collection.lte_ionization_properties
            )
        elif self.ionization_analytical_approximation == "nebular":
            self.plasma_modules += (
                self.plasma_collection.nebular_ionization_properties
            )
        else:
            raise PlasmaConfigError(
                f'Invalid excitation analytical approximation. Configured as {self.ionization_analytical_approximation} but needs to be either "lte" or "nebular"'
            )

    def initialize_j_blues(self, dilute_planckian_radiation_field, lines_df):
        """
        Initialize the j_blues DataFrame based on the radiative_rates_type and the dilute_planckian_radiation_field.

        Parameters
        ----------
        dilute_planckian_radiation_field : object
            The dilute Planckian radiation field object.
        lines_df : pandas.DataFrame
            The DataFrame containing lines information.

        Returns
        -------
        pandas.DataFrame
            The DataFrame with calculated mean intensity values.

        Raises
        ------
        ValueError
            If the radiative_rates_type is unknown.
        """
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
        self.continuum_interaction_species = [
            species_string_to_tuple(species)
            for species in continuum_interaction_species
        ]

    def check_continuum_interaction_species(self):
        """
        Check if all continuum interaction species belong to atoms that have been specified in the configuration.

        Raises
        ------
            PlasmaConfigError: If not all continuum interaction species belong to specified atoms.
        """
        continuum_atoms = (
            self.continuum_interaction_species_multi_index.get_level_values(
                "atomic_number"
            )
        )

        continuum_atoms_in_selected_atoms = np.all(
            continuum_atoms.isin(self.atom_data.selected_atomic_numbers)
        )

        if not continuum_atoms_in_selected_atoms:
            raise PlasmaConfigError(
                "Not all continuum interaction species "
                "belong to atoms that have been specified "
                "in the configuration."
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

    def setup_continuum_interactions(self):
        """
        Set up continuum interactions for the plasma assembly.

        Raises
        ------
            PlasmaConfigError: If the line_interaction_type is not "macroatom".
            PlasmaConfigError: If an NLTE ionization species is not in the continuum species.
            PlasmaConfigError: If an NLTE excitation species is not in the continuum species.
            PlasmaConfigError: If the NLTE solver type is unknown.
        """
        if self.line_interaction_type != "macroatom":
            raise PlasmaConfigError(
                "Continuum interactions require line_interaction_type "
                f"macroatom (instead of {self.line_interaction_type})."
            )

        self.plasma_modules += (
            self.plasma_collection.continuum_interaction_properties
        )
        self.plasma_modules += (
            self.plasma_collection.continuum_interaction_inputs
        )

        if self.enable_adiabatic_cooling:
            self.plasma_modules += (
                self.plasma_collection.adiabatic_cooling_properties
            )

        if self.enable_two_photon_decay:
            self.plasma_modules += self.plasma_collection.two_photon_properties

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
                    if species not in self.continuum_interaction_species:
                        raise PlasmaConfigError(
                            f"NLTE ionization species {species} not in continuum species."
                        )
            if self.nlte_excitation_species:
                nlte_excitation_species = self.nlte_excitation_species
                for species in nlte_excitation_species:
                    if species not in self.continuum_interaction_species:
                        raise PlasmaConfigError(
                            f"NLTE excitation species {species} not in continuum species."
                        )
            self.property_kwargs[NLTEIndexHelper] = {
                "nlte_ionization_species": self.nlte_ionization_species,
                "nlte_excitation_species": self.nlte_excitation_species,
            }
            if self.nlte_solver == "lu":
                self.plasma_modules += (
                    self.plasma_collection.nlte_lu_solver_properties
                )
                logger.warning(
                    "LU solver will be inaccurate for NLTE excitation, proceed with caution."
                )
            elif self.nlte_solver == "root":
                self.plasma_modules += (
                    self.plasma_collection.nlte_root_solver_properties
                )
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

    def assemble(
        self,
        number_densities,
        dilute_planckian_radiation_field,
        time_explosion,
        electron_densities=None,
        **kwargs,
    ):
        """
        Assemble the plasma based on the provided parameters and settings.

        Parameters
        ----------
        number_densities : dict
            Dictionary of number densities for different species.
        dilute_planckian_radiation_field : object
            The dilute Planckian radiation field object.
        time_explosion : float
            The time of explosion.
        electron_densities : array-like, optional
            Optional electron densities.

        Returns
        -------
        BasePlasma
            The assembled plasma object.

        Raises
        ------
        ValueError
            If an error occurs during assembly.
        """
        j_blues = self.initialize_j_blues(
            dilute_planckian_radiation_field, self.atom_data.lines
        )
        plasma_solver_settings = PlasmaSolverSettings(
            RADIATIVE_RATES_TYPE=self.radiative_rates_type
        )

        plasma_assemble_kwargs = dict(
            time_explosion=time_explosion,
            dilute_planckian_radiation_field=dilute_planckian_radiation_field,
            number_density=number_densities,
            link_t_rad_t_electron=self.link_t_rad_t_electron,
            atomic_data=self.atom_data,
            j_blues=j_blues,
            continuum_interaction_species=self.continuum_interaction_species_multi_index,
            nlte_ionization_species=self.nlte_ionization_species,
            nlte_excitation_species=self.nlte_excitation_species,
        )
        if len(self.continuum_interaction_species) > 0:
            initial_continuum_properties = self.initialize_continuum_properties(
                dilute_planckian_radiation_field
            )
            plasma_assemble_kwargs.update(
                gamma=initial_continuum_properties.photo_ionization_rate_coefficient,
                bf_heating_coeff_estimator=None,
                stim_recomb_cooling_coeff_estimator=None,
                alpha_stim_factor=initial_continuum_properties.stimulated_recombination_rate_factor,
            )

        if electron_densities is not None:
            electron_densities = pd.Series(electron_densities.cgs.value)

        self.setup_electron_densities(electron_densities)
        plasma_assemble_kwargs["helium_treatment"] = self.helium_treatment
        plasma_assemble_kwargs.update(kwargs)
        return BasePlasma(
            plasma_properties=self.plasma_modules,
            property_kwargs=self.property_kwargs,
            plasma_solver_settings=plasma_solver_settings,
            **plasma_assemble_kwargs,
        )
