import numpy as np
import pandas as pd

from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.properties.property_collections import (
    basic_inputs,
    basic_properties,
    dilute_lte_excitation_properties,
    lte_excitation_properties,
    lte_ionization_properties,
    nebular_ionization_properties,
    nlte_properties,
    non_nlte_properties,
)
from tardis.util.base import species_string_to_tuple


from tardis.plasma.properties import (
    HeliumNumericalNLTE,
    IonNumberDensity,
    IonNumberDensityHeNLTE,
    LevelBoltzmannFactorNLTE,
    MarkovChainTransProbsCollector,
    RadiationFieldCorrection,
    StimulatedEmissionFactor,
)


def map_species_from_string(species):
    return [species_string_to_tuple(spec) for spec in species]


class PlasmaSolverFactory:

    continuum_interaction_species: pd.MultiIndex
    line_interaction_type: str = "scatter"

    legacy_nlte_species: list = []
    nlte_exciation_species: list = []
    nlte_ionization_species: list = []
    plasma_modules: list = []
    kwargs: dict = {}
    property_kwargs: dict = {}

    excitation_analytical_approximation: str = "lte"
    ionization_analytical_approximation: str = "lte"

    radiative_rates_type: str = "dilute-blackbody"

    nebular_ionization_delta_treatment: tuple  # species to use for the delta_treatment in nebular ionization ML93

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

        self.kwargs = dict(
            link_t_rad_t_electron=config.plasma.link_t_rad_t_electron,
            continuum_interaction_species=continuum_interaction_species,
            nlte_ionization_species=nlte_ionization_species,
            nlte_excitation_species=nlte_excitation_species,
        )

        self.setup_analytical_approximations(config)

        self.setup_legacy_nlte(config.plasma.nlte)
        if self.line_interaction_type in ("downbranch", "macroatom") and (
            len(self.continuum_interaction_species) == 0
        ):
            self.setup_legacy_macro_atom(config)

    def setup_legacy_macro_atom(self, macro_atom_config=None):
        self.plasma_modules += macro_atom_properties

        if macro_atom_config is not None:
            self.plasma_modules.append(
                MarkovChainTransProbsCollector.from_config(macro_atom_config)
            )

    def setup_legacy_nlte(self, nlte_config):
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
        plasma_modules = []
        if self.excitation_analytical_approximation == "lte":
            plasma_modules += lte_excitation_properties
        elif self.excitation_analytical_approximation == "dilute-lte":
            plasma_modules += dilute_lte_excitation_properties
        else:
            raise PlasmaConfigError(
                f'Invalid excitation analytical approximation. Configured as {self.excitation_analytical_approximation} but needs to be either "lte" or "dilute-lte"'
            )

        if self.ionization_analytical_approximation == "lte":
            plasma_modules += lte_ionization_properties
        elif self.ionization_analytical_approximation == "nebular":
            plasma_modules += nebular_ionization_properties
            if "delta_treatment" in plasma_config:
                self.property_kwargs[RadiationFieldCorrection] = dict(
                    delta_treatment=plasma_config.delta_treatment
                )
        else:
            raise PlasmaConfigError(
                f'Invalid excitation analytical approximation. Configured as {self.ionization_analytical_approximation} but needs to be either "lte" or "nebular"'
            )

    def setup_radiative_rates(self, radiative_rates_config):
        ##### RADIATIVE RATES SETUP

        if (self.radiative_rates_type == "dilute-blackbody") or (
            self.radiative_rates_type == "detailed"
        ):
            self.kwargs["j_blues"] = pd.DataFrame(
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

            self.radiative_rates_type = radiative_rates_config.type

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

    def assemble(
        self, dilute_planckian_radiation_field, simulation_state, atom_data
    ):

        kwargs = dict(
            time_explosion=simulation_state.time_explosion,
            dilute_planckian_radiation_field=dilute_planckian_radiation_field,
            abundance=simulation_state.abundance,
            number_density=simulation_state.elemental_number_density,
            atomic_data=atom_data,
        )

        return BasePlasma(
            plasma_properties=self.plasma_modules,
            property_kwargs=self.property_kwargs,
            **kwargs,
        )
