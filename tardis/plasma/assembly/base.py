import importlib

import numpy as np
import pandas as pd
from astropy import units as u

from tardis.plasma import BasePlasma
from tardis.plasma.base import PlasmaSolverSettings
from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.properties import (
    HeliumNumericalNLTE,
    HydrogenContinuumFractionalHeating,
    IonNumberDensity,
    IonNumberDensityHeNLTE,
    LevelBoltzmannFactorNLTE,
    RadiationFieldCorrection,
    StimulatedEmissionFactor,
)
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.util.base import species_string_to_tuple


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

    ## Helium Treatment options
    helium_treatment: str = "none"
    heating_rate_data_file: str = "none"

    ## Continuum Interaction
    continuum_interaction_species: list = []

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

        self.radiative_rates_type = plasma_config.radiative_rates_type

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

        if self.continuum_interaction_species:
            raise PlasmaConfigError(
                "Continuum interactions are supported only by the IIP workflow."
            )

        self.atom_data.prepare_atom_data(
            selected_atomic_numbers,
            line_interaction_type=self.line_interaction_type,
            continuum_interaction_species=self.continuum_interaction_species_multi_index,
            nlte_species=self.legacy_nlte_species,
        )

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

        if self.line_interaction_type in ("downbranch", "macroatom"):
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

        """
        if self.helium_treatment == "recomb-nlte":
            self.plasma_modules += self.plasma_collection.helium_nlte_properties
        elif self.helium_treatment == "numerical-nlte":
            self.plasma_modules += (
                self.plasma_collection.helium_numerical_nlte_properties
            )
            if self.heating_rate_data_file in ["none", None]:
                raise PlasmaConfigError("Heating rate data file not specified")
            self.property_kwargs[HeliumNumericalNLTE] = dict(
                heating_rate_data_file=self.heating_rate_data_file
            )
        else:
            self.plasma_modules += self.plasma_collection.helium_lte_properties

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

    def setup_electron_densities(self, electron_densities):
        if self.helium_treatment == "numerical-nlte":
            self.property_kwargs[IonNumberDensityHeNLTE] = dict(
                electron_densities=electron_densities
            )
        else:
            self.property_kwargs[IonNumberDensity] = dict(
                electron_densities=electron_densities
            )

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


class IIPPlasmaSolverFactory:
    """Simplified plasma solver factory that produces only IIP plasmas."""

    excitation_analytical_approximation = "dilute-lte"
    ionization_analytical_approximation = "nebular"

    link_t_rad_t_electron: float = 1.0
    radiative_rates_type: str = "dilute-blackbody"
    line_interaction_type: str = "macroatom"

    ## Continuum Interaction
    continuum_interaction_species: list = []

    plasma_modules: list = []
    property_kwargs: dict = {}
    plasma_collection = None

    def __init__(self, atom_data, config=None) -> None:
        self.plasma_modules = []
        self.property_kwargs = {}
        self.atom_data = atom_data
        self.plasma_collection = importlib.import_module(
            "tardis.plasma.properties.iip_property_collections"
        )
        if config is not None:
            self.link_t_rad_t_electron = config.plasma.link_t_rad_t_electron
            self.radiative_rates_type = config.plasma.radiative_rates_type
            self.line_interaction_type = config.plasma.line_interaction_type
            self.continuum_interaction_species = (
                config.plasma.continuum_interaction.species
            )
            self.nlte_species = config.plasma.nlte.species

    @property
    def continuum_interaction_species_multi_index(self):
        return convert_species_to_multi_index(
            self.continuum_interaction_species
        )

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
        if len(self.continuum_interaction_species) == 0:
            return

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

    def prepare_factory(
        self,
        selected_atomic_numbers,
        property_collections,
        nlte_species,
        config=None,
    ):
        """Set up the IIP plasma factory with minimal configuration."""
        self.plasma_collection = importlib.import_module(property_collections)

        self.atom_data.prepare_atom_data(
            selected_atomic_numbers,
            line_interaction_type=self.line_interaction_type,
            continuum_interaction_species=self.continuum_interaction_species_multi_index,
            nlte_species=nlte_species,
        )

        self.check_continuum_interaction_species()

        self.plasma_modules = (
            self.plasma_collection.basic_inputs
            + self.plasma_collection.basic_properties
        )

        self.plasma_modules += (
            self.plasma_collection.dilute_lte_excitation_properties
        )
        self.plasma_modules += (
            self.plasma_collection.nebular_ionization_properties
        )

        self.plasma_modules += self.plasma_collection.non_nlte_properties

        self.plasma_modules += self.plasma_collection.macro_atom_properties

        # Hydrogen continuum properties
        self.plasma_modules += self.plasma_collection.hydrogen_continuum_inputs
        self.plasma_modules += (
            self.plasma_collection.hydrogen_continuum_properties
        )

        # Set up property kwargs for hydrogen continuum
        hydrogen_photoionization_kwargs = dict(
            photo_ion_cross_sections=self.atom_data.photoionization_data
        )
        self.property_kwargs[HydrogenContinuumFractionalHeating] = (
            hydrogen_photoionization_kwargs
        )
        self.property_kwargs[StimulatedEmissionFactor] = dict(
            nlte_species=set(nlte_species)
        )

    def initialize_j_blues(self, dilute_planckian_radiation_field, lines_df):
        """Initialize j_blues for the given radiation field."""
        j_blues = pd.DataFrame(
            dilute_planckian_radiation_field.calculate_mean_intensity(
                lines_df.nu.values
            ),
            index=lines_df.index,
        )
        return j_blues

    def assemble(
        self,
        number_densities,
        dilute_planckian_radiation_field,
        time_explosion,
        electron_densities=None,
        **kwargs,
    ):
        """Assemble a IIP plasma."""
        plasma_solver_settings = PlasmaSolverSettings(
            RADIATIVE_RATES_TYPE=self.radiative_rates_type
        )

        plasma_assemble_kwargs = dict(
            time_explosion=time_explosion,
            dilute_planckian_radiation_field=dilute_planckian_radiation_field,
            number_density=number_densities,
            link_t_rad_t_electron=self.link_t_rad_t_electron,
            atomic_data=self.atom_data,
            j_blues=None,
            continuum_interaction_species=self.continuum_interaction_species_multi_index,
        )

        if electron_densities is not None:
            electron_densities = pd.Series(electron_densities.cgs.value)
            self.property_kwargs[IonNumberDensity] = dict(
                electron_densities=electron_densities
            )

        plasma_assemble_kwargs.update(kwargs)
        return IIPPlasma(
            plasma_properties=self.plasma_modules,
            property_kwargs=self.property_kwargs,
            plasma_solver_settings=plasma_solver_settings,
            **plasma_assemble_kwargs,
        )


class IIPPlasma(BasePlasma):
    """Standard ``BasePlasma`` with the IIP iteration update contract."""

    def update_radiationfield(
        self,
        t_rad,
        dilution_factor,
        j_blues,
        photoionization_rate_estimator=None,
        stimulated_recombination_rate_estimator=None,
        bound_free_heating_estimator=None,
        stimulated_recombination_cooling_estimator=None,
        free_free_heating_estimator=None,
        initialize=False,
    ):
        """Update the radiation field and estimator-backed equilibrium state."""
        if not isinstance(t_rad, u.Quantity):
            t_rad = t_rad * u.K
        if isinstance(j_blues, pd.DataFrame):
            j_blues = j_blues.copy()
            j_blues.columns = np.arange(j_blues.shape[1])
        iteration = self.get_value("iteration") or 0
        next_iteration = 0 if initialize else iteration + 1
        update_kwargs = dict(
            dilute_planckian_radiation_field=DilutePlanckianRadiationField(
                t_rad, dilution_factor
            ),
            j_blues=j_blues,
            iteration=next_iteration,
            previous_beta_sobolev=self.get_value("beta_sobolev").copy(),
            previous_electron_densities=self.get_value("electron_densities"),
            previous_ion_number_density=self.get_value(
                "ion_number_density"
            ).copy(),
            previous_level_number_density=self.get_value(
                "level_number_density"
            ).copy(),
        )
        update_kwargs.update(
            photoionization_rate_estimator=photoionization_rate_estimator,
            stimulated_recombination_rate_estimator=stimulated_recombination_rate_estimator,
            bound_free_heating_estimator=bound_free_heating_estimator,
            stimulated_recombination_cooling_estimator=stimulated_recombination_cooling_estimator,
            free_free_heating_estimator=free_free_heating_estimator,
        )
        self.update(**update_kwargs)
