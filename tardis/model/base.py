import logging
import os
from pathlib import Path

import numpy as np
from astropy import units as u

from tardis import constants
from tardis.io.configuration.config_reader import Configuration
from tardis.io.configuration.config_validator import validate_dict
from tardis.io.model.readers.csvy import (
    load_csvy,
)
from tardis.io.util import HDFWriterMixin
from tardis.model.matter.composition import Composition
from tardis.model.parse_input import (
    parse_abundance_config,
    parse_csvy_composition,
    parse_csvy_geometry,
    parse_csvy_radiation_field_state,
    parse_radiation_field_state,
    parse_structure_config,
    parse_packet_source,
)
from tardis.montecarlo.packet_source import BlackBodySimpleSource
from tardis.model.radiation_field_state import (
    DiluteBlackBodyRadiationFieldState,
)
from tardis.util.base import is_valid_nuclide_or_elem

logger = logging.getLogger(__name__)


SCHEMA_DIR = (
    Path(__file__).parent / ".." / "io" / "configuration" / "schemas"
).resolve()


class SimulationState(HDFWriterMixin):
    """
    An object that hold information about the individual shells.

    Parameters
    ----------
    velocity : np.ndarray
        An array with n+1 (for n shells) velocities "cut" to the provided
        boundaries

        .. note:: To access the entire, "uncut", velocity array, use `raw_velocity`
    abundance : pd.DataFrame
    time_explosion : astropy.units.Quantity
        Time since explosion
    t_inner : astropy.units.Quantity
    elemental_mass: pd.Series
    luminosity_requested : astropy.units.quantity.Quantity
    t_radiative : astropy.units.Quantity
        Radiative temperature for the shells
    dilution_factor : np.ndarray
        If None, the dilution_factor will be initialized with the geometric
        dilution factor.
    v_boundary_inner : astropy.units.Quantity
    v_boundary_outer : astropy.units.Quantity
    raw_velocity : np.ndarray
        The complete array of the velocities, without being cut by
        `v_boundary_inner` and `v_boundary_outer`
    electron_densities : astropy.units.quantity.Quantity

    Attributes
    ----------
    w : numpy.ndarray
        Shortcut for `dilution_factor`
    t_rad : astropy.units.quantity.Quantity
        Shortcut for `t_radiative`
    radius : astropy.units.quantity.Quantity
    r_inner : astropy.units.quantity.Quantity
    r_outer : astropy.units.quantity.Quantity
    r_middle : astropy.units.quantity.Quantity
    v_inner : astropy.units.quantity.Quantity
    v_outer : astropy.units.quantity.Quantity
    v_middle : astropy.units.quantity.Quantity
    density : astropy.units.quantity.Quantity
    volume : astropy.units.quantity.Quantity
    no_of_shells : int
        The number of shells as formed by `v_boundary_inner` and
        `v_boundary_outer`
    no_of_raw_shells : int
    """

    hdf_properties = [
        "t_inner",
        "dilution_factor",
        "t_radiative",
        "v_inner",
        "v_outer",
        "density",
        "r_inner",
        "time_explosion",
    ]
    hdf_name = "simulation_state"

    def __init__(
        self,
        geometry,
        composition,
        radiation_field_state,
        time_explosion,
        packet_source,
        electron_densities=None,
    ):
        self.geometry = geometry
        self.composition = composition
        self.radiation_field_state = radiation_field_state
        self.time_explosion = time_explosion
        self._electron_densities = electron_densities
        self.packet_source = packet_source

    @property
    def blackbody_packet_source(self):
        return self.packet_source

    @property
    def t_inner(self):
        return self.packet_source.temperature

    @t_inner.setter
    def t_inner(self, value):
        self.packet_source.temperature = value

    @property
    def dilution_factor(self):
        return self.radiation_field_state.dilution_factor[
            self.geometry.v_inner_boundary_index : self.geometry.v_outer_boundary_index
        ]

    @dilution_factor.setter
    def dilution_factor(self, new_dilution_factor):
        if len(new_dilution_factor) == self.no_of_shells:
            self.radiation_field_state.dilution_factor[
                self.geometry.v_inner_boundary_index : self.geometry.v_outer_boundary_index
            ] = new_dilution_factor
        else:
            raise ValueError(
                "Trying to set dilution_factor for unmatching number"
                "of shells."
            )

    @property
    def t_radiative(self):
        return self.radiation_field_state.t_radiative[
            self.geometry.v_inner_boundary_index : self.geometry.v_outer_boundary_index
        ]

    @t_radiative.setter
    def t_radiative(self, new_t_radiative):
        if len(new_t_radiative) == self.no_of_shells:
            self.radiation_field_state.t_radiative[
                self.geometry.v_inner_boundary_index : self.geometry.v_outer_boundary_index
            ] = new_t_radiative
        else:
            raise ValueError(
                "Trying to set t_radiative for different number of shells."
            )

    @property
    def radius(self):
        return self.time_explosion * self.velocity

    @property
    def v_boundary_inner(self):
        return self.geometry.v_inner_boundary

    @property
    def v_boundary_outer(self):
        return self.geometry.v_outer_boundary

    @property
    def r_inner(self):
        return self.geometry.r_inner_active

    @property
    def r_outer(self):
        return self.geometry.r_outer_active

    @property
    def r_middle(self):
        return 0.5 * self.r_inner + 0.5 * self.r_outer

    @property
    def velocity(self):
        velocity = self.geometry.v_outer_active.copy()
        return velocity.insert(0, self.geometry.v_inner_active[0])

    @property
    def v_inner(self):
        return self.geometry.v_inner_active

    @property
    def v_outer(self):
        return self.geometry.v_outer_active

    @property
    def v_middle(self):
        return 0.5 * self.v_inner + 0.5 * self.v_outer

    @property
    def density(self):
        return self.composition.density[
            self.geometry.v_inner_boundary_index : self.geometry.v_outer_boundary_index
        ]

    @property
    def abundance(self):
        elemental_mass_fraction = (
            self.composition.elemental_mass_fraction.copy()
        )
        elemental_mass_fraction = elemental_mass_fraction.iloc[
            :,
            self.geometry.v_inner_boundary_index : self.geometry.v_outer_boundary_index,
        ]
        elemental_mass_fraction.columns = range(
            len(elemental_mass_fraction.columns)
        )
        return elemental_mass_fraction

    @property
    def volume(self):
        return ((4.0 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)).cgs

    @property
    def no_of_shells(self):
        return self.geometry.no_of_shells_active

    @property
    def no_of_raw_shells(self):
        return self.geometry.no_of_shells

    @classmethod
    def from_config(cls, config, atom_data, legacy_mode_enabled=False):
        """
        Create a new SimulationState instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        atom_data : tardis.io.AtomData

        Returns
        -------
        SimulationState
        """
        time_explosion = config.supernova.time_explosion.cgs

        (
            electron_densities,
            t_radiative,
            geometry,
            density,
        ) = parse_structure_config(config, time_explosion)

        nuclide_mass_fraction, raw_isotope_abundance = parse_abundance_config(
            config, geometry, time_explosion
        )

        # using atom_data.mass.copy() to ensure that the original atom_data is not modified
        composition = Composition(
            density,
            nuclide_mass_fraction,
            raw_isotope_abundance,
            atom_data.atom_data.mass.copy(),
        )

        packet_source = parse_packet_source(
            config, geometry, legacy_mode_enabled
        )
        radiation_field_state = parse_radiation_field_state(
            config,
            t_radiative,
            geometry,
            dilution_factor=None,
            packet_source=packet_source,
        )

        return cls(
            geometry=geometry,
            composition=composition,
            radiation_field_state=radiation_field_state,
            time_explosion=time_explosion,
            packet_source=packet_source,
            electron_densities=electron_densities,
        )

    @classmethod
    def from_csvy(cls, config, atom_data=None, legacy_mode_enabled=False):
        """
        Create a new SimulationState instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        atom_data : tardis.io.AtomData

        Returns
        -------
        SimulationState
        """
        CSVY_SUPPORTED_COLUMNS = {
            "velocity",
            "density",
            "t_rad",
            "dilution_factor",
        }

        if os.path.isabs(config.csvy_model):
            csvy_model_fname = config.csvy_model
        else:
            csvy_model_fname = os.path.join(
                config.config_dirname, config.csvy_model
            )
        csvy_model_config, csvy_model_data = load_csvy(csvy_model_fname)
        csvy_schema_fname = SCHEMA_DIR / "csvy_model.yml"
        csvy_model_config = Configuration(
            validate_dict(csvy_model_config, schemapath=csvy_schema_fname)
        )

        if hasattr(csvy_model_data, "columns"):
            abund_names = set(
                [
                    name
                    for name in csvy_model_data.columns
                    if is_valid_nuclide_or_elem(name)
                ]
            )
            unsupported_columns = (
                set(csvy_model_data.columns)
                - abund_names
                - CSVY_SUPPORTED_COLUMNS
            )

            field_names = set(
                [field["name"] for field in csvy_model_config.datatype.fields]
            )
            assert (
                set(csvy_model_data.columns) - field_names == set()
            ), "CSVY columns exist without field descriptions"
            assert (
                field_names - set(csvy_model_data.columns) == set()
            ), "CSVY field descriptions exist without corresponding csv data"
            if unsupported_columns != set():
                logger.warning(
                    "The following columns are "
                    "specified in the csvy model file,"
                    f" but are IGNORED by TARDIS: {str(unsupported_columns)}"
                )

        time_explosion = config.supernova.time_explosion.cgs

        electron_densities = None

        geometry = parse_csvy_geometry(
            config, csvy_model_config, csvy_model_data, time_explosion
        )

        composition = parse_csvy_composition(
            atom_data,
            csvy_model_config,
            csvy_model_data,
            time_explosion,
            geometry,
        )

        packet_source = parse_packet_source(
            config, geometry, legacy_mode_enabled
        )

        radiation_field_state = parse_csvy_radiation_field_state(
            config, csvy_model_config, csvy_model_data, geometry, packet_source
        )

        return cls(
            geometry=geometry,
            composition=composition,
            time_explosion=time_explosion,
            radiation_field_state=radiation_field_state,
            packet_source=packet_source,
            electron_densities=electron_densities,
        )
