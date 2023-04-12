import os
import logging
import numpy as np
import pandas as pd
from astropy import units as u
from tardis import constants
import radioactivedecay as rd
from radioactivedecay.utils import Z_DICT
from tardis.model.geometry.radial1d import Radial1DGeometry

from tardis.util.base import quantity_linspace, is_valid_nuclide_or_elem
from tardis.io.parsers.csvy import load_csvy
from tardis.io.model_reader import (
    read_density_file,
    read_abundances_file,
    read_uniform_abundances,
    parse_csv_abundances,
)
from tardis.io.config_validator import validate_dict
from tardis.io.config_reader import Configuration
from tardis.io.util import HDFWriterMixin
from tardis.io.decay import IsotopeAbundances
from tardis.model.density import HomologousDensity

logger = logging.getLogger(__name__)


class Composition:
    """
    Holds information about model composition

    Parameters
    ----------
    density : astropy.units.quantity.Quantity
        An array of densities for each shell.
    isotopic_mass_fraction : pd.DataFrame
    atomic_mass : pd.DataFrame
    atomic_mass_unit: astropy.units.Unit

    Attributes
    ----------
    atomic_mass : pd.DataFrame
        Atomic mass of elements calculated for each shell.
    elemental_number_density : pd.DataFrame
        Number density of each element in each shell.
    """

    def __init__(
        self,
        density,
        elemental_mass_fraction,
        atomic_mass,
        atomic_mass_unit=u.g,
    ):
        self.density = density
        self.elemental_mass_fraction = elemental_mass_fraction
        self.atomic_mass_unit = atomic_mass_unit
        self._atomic_mass = atomic_mass

    @property
    def atomic_mass(self):
        """Atomic mass of elements in each shell"""
        if self._atomic_mass is None:
            raise AttributeError(
                "ModelState was not provided elemental masses."
            )
        return self._atomic_mass

    @property
    def elemental_number_density(self):
        """Elemental Number Density computed using the formula: (elemental_mass_fraction * density) / atomic mass"""
        if self.atomic_mass is None:
            raise AttributeError(
                "ModelState was not provided elemental masses."
            )
        return (self.elemental_mass_fraction * self.density).divide(
            self.atomic_mass, axis=0
        )


class ModelState:
    """
    Holds information about model geometry for radial 1D models.

    Parameters
    ----------
    composition : tardis.model.Composition
    geometry : tardis.model.geometry.radial1d.Radial1DGeometry
    time_explosion : astropy.units.quantity.Quantity

    Attributes
    ----------
    mass : pd.DataFrame
    number : pd.DataFrame
    """

    def __init__(self, composition, geometry, time_explosion):
        self.time_explosion = time_explosion
        self.composition = composition
        self.geometry = geometry

    @property
    def mass(self):
        """Mass calculated using the formula:
        mass_fraction * density * volume"""
        return (
            self.composition.elemental_mass_fraction
            * self.composition.density
            * self.geometry.volume
        )

    @property
    def number(self):
        """Number calculated using the formula:
        mass / atomic_mass"""
        if self.composition.atomic_mass is None:
            raise AttributeError(
                "ModelState was not provided elemental masses."
            )
        return (self.mass).divide(self.composition.atomic_mass, axis=0)


class Radial1DModel(HDFWriterMixin):
    """
    An object that hold information about the individual shells.

    Parameters
    ----------
    velocity : np.ndarray
        An array with n+1 (for n shells) velocities "cut" to the provided
        boundaries

        .. note:: To access the entire, "uncut", velocity array, use `raw_velocity`
    homologous_density : HomologousDensity
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
        "w",
        "t_radiative",
        "v_inner",
        "v_outer",
        "homologous_density",
        "r_inner",
    ]
    hdf_name = "model"

    def __init__(
        self,
        velocity,
        homologous_density,
        abundance,
        isotope_abundance,
        time_explosion,
        t_inner,
        elemental_mass,
        luminosity_requested=None,
        t_radiative=None,
        dilution_factor=None,
        v_boundary_inner=None,
        v_boundary_outer=None,
        electron_densities=None,
    ):
        self._v_boundary_inner = None
        self._v_boundary_outer = None
        self._velocity = None
        self.raw_velocity = velocity
        self.v_boundary_inner = v_boundary_inner
        self.v_boundary_outer = v_boundary_outer
        self.homologous_density = homologous_density
        self._abundance = abundance
        self.time_explosion = time_explosion
        self._electron_densities = electron_densities
        v_outer = self.velocity[1:]
        v_inner = self.velocity[:-1]
        density = (
            self.homologous_density.calculate_density_at_time_of_simulation(
                self.time_explosion
            )[self.v_boundary_inner_index + 1 : self.v_boundary_outer_index + 1]
        )
        self.raw_abundance = self._abundance
        self.raw_isotope_abundance = isotope_abundance

        atomic_mass = None
        if elemental_mass is not None:
            mass = {}
            stable_atomic_numbers = self.raw_abundance.index.to_list()
            for z in stable_atomic_numbers:
                mass[z] = [
                    elemental_mass[z]
                    for i in range(self.raw_abundance.columns.size)
                ]
            stable_isotope_mass = pd.DataFrame(mass).T

            isotope_mass = {}
            for atomic_number, i in self.raw_isotope_abundance.decay(
                self.time_explosion
            ).groupby(level=0):
                i = i.loc[atomic_number]
                for column in i:
                    mass = {}
                    shell_abundances = i[column]
                    isotopic_masses = [
                        rd.Nuclide(Z_DICT[atomic_number] + str(i)).atomic_mass
                        for i in shell_abundances.index.to_numpy()
                    ]
                    mass[atomic_number] = (
                        shell_abundances * isotopic_masses
                    ).sum()
                    mass[atomic_number] /= shell_abundances.sum()
                    mass[atomic_number] = mass[atomic_number] * u.u.to(u.g)
                    if isotope_mass.get(column) is None:
                        isotope_mass[column] = {}
                    isotope_mass[column][atomic_number] = mass[atomic_number]
            isotope_mass = pd.DataFrame(isotope_mass)

            atomic_mass = pd.concat([stable_isotope_mass, isotope_mass])

        composition = Composition(
            density=density,
            elemental_mass_fraction=self.abundance,
            atomic_mass=atomic_mass,
        )
        geometry = Radial1DGeometry(
            r_inner=self.time_explosion * v_inner,
            r_outer=self.time_explosion * v_outer,
            v_inner=v_inner,
            v_outer=v_outer,
        )
        self.model_state = ModelState(
            composition=composition,
            geometry=geometry,
            time_explosion=self.time_explosion,
        )

        if t_inner is None:
            if luminosity_requested is not None:
                self.t_inner = (
                    (
                        luminosity_requested
                        / (
                            4
                            * np.pi
                            * self.r_inner[0] ** 2
                            * constants.sigma_sb
                        )
                    )
                    ** 0.25
                ).to("K")
            else:
                raise ValueError(
                    "Both t_inner and luminosity_requested cannot " "be None."
                )
        else:
            self.t_inner = t_inner

        if t_radiative is None:
            lambda_wien_inner = constants.b_wien / self.t_inner
            self._t_radiative = constants.b_wien / (
                lambda_wien_inner
                * (1 + (self.v_middle - self.v_boundary_inner) / constants.c)
            )
        else:
            # self._t_radiative = t_radiative[self.v_boundary_inner_index + 1:self.v_boundary_outer_index]
            self._t_radiative = t_radiative

        if dilution_factor is None:
            self._dilution_factor = 0.5 * (
                1
                - np.sqrt(
                    1 - (self.r_inner[0] ** 2 / self.r_middle**2).to(1).value
                )
            )
        else:
            # self.dilution_factor = dilution_factor[self.v_boundary_inner_index + 1:self.v_boundary_outer_index]
            self._dilution_factor = dilution_factor

    @property
    def w(self):
        return self.dilution_factor

    @w.setter
    def w(self, value):
        self.dilution_factor = value

    @property
    def t_rad(self):
        return self.t_radiative

    @t_rad.setter
    def t_rad(self, value):
        self.t_radiative = value

    @property
    def dilution_factor(self):
        if len(self._dilution_factor) == self.no_of_shells:
            return self._dilution_factor

        #        if self.v_boundary_inner in self.raw_velocity:
        #            v_inner_ind = np.argwhere(self.raw_velocity == self.v_boundary_inner)[0][0]
        #        else:
        #            v_inner_ind = np.searchsorted(self.raw_velocity, self.v_boundary_inner) - 1
        #        if self.v_boundary_outer in self.raw_velocity:
        #            v_outer_ind = np.argwhere(self.raw_velocity == self.v_boundary_outer)[0][0]
        #        else:
        #            v_outer_ind = np.searchsorted(self.raw_velocity, self.v_boundary_outer)

        return self._dilution_factor[
            self.v_boundary_inner_index + 1 : self.v_boundary_outer_index + 1
        ]

    @dilution_factor.setter
    def dilution_factor(self, value):
        if len(value) == len(self._dilution_factor):
            self._dilution_factor = value
        elif len(value) == self.no_of_shells:
            #            if self.v_boundary_inner in self.raw_velocity:
            #                v_inner_ind = np.argwhere(self.raw_velocity == self.v_boundary_inner)[0][0]
            #            else:
            #                v_inner_ind = np.searchsorted(self.raw_velocity, self.v_boundary_inner) - 1
            #            if self.v_boundary_outer in self.raw_velocity:
            #                v_outer_ind = np.argwhere(self.raw_velocity == self.v_boundary_outer)[0][0]
            #            else:
            #                v_outer_ind = np.searchsorted(self.raw_velocity, self.v_boundary_outer)
            #            assert v_outer_ind - v_inner_ind == self.no_of_shells, "trad shape different from number of shells"
            self._dilution_factor[
                self.v_boundary_inner_index
                + 1 : self.v_boundary_outer_index
                + 1
            ] = value
        else:
            raise ValueError(
                "Trying to set dilution_factor for unmatching number"
                "of shells."
            )

    @property
    def t_radiative(self):
        if len(self._t_radiative) == self.no_of_shells:
            return self._t_radiative

        #        if self.v_boundary_inner in self.raw_velocity:
        #            v_inner_ind = np.argwhere(self.raw_velocity == self.v_boundary_inner)[0][0]
        #        else:
        #            v_inner_ind = np.searchsorted(self.raw_velocity, self.v_boundary_inner) - 1
        #        if self.v_boundary_outer in self.raw_velocity:
        #            v_outer_ind = np.argwhere(self.raw_velocity == self.v_boundary_outer)[0][0]
        #        else:
        #            v_outer_ind = np.searchsorted(self.raw_velocity, self.v_boundary_outer)

        return self._t_radiative[
            self.v_boundary_inner_index + 1 : self.v_boundary_outer_index + 1
        ]

    @t_radiative.setter
    def t_radiative(self, value):
        if len(value) == len(self._t_radiative):
            self._t_radiative = value
        elif len(value) == self.no_of_shells:
            #            if self.v_boundary_inner in self.raw_velocity:
            #                v_inner_ind = np.argwhere(self.raw_velocity == self.v_boundary_inner)[0][0]
            #            else:
            #                v_inner_ind = np.searchsorted(self.raw_velocity, self.v_boundary_inner) - 1
            #            if self.v_boundary_outer in self.raw_velocity:
            #                v_outer_ind = np.argwhere(self.raw_velocity == self.v_boundary_outer)[0][0]
            #            else:
            #                v_outer_ind = np.searchsorted(self.raw_velocity, self.v_boundary_outer)
            #            assert v_outer_ind - v_inner_ind == self.no_of_shells, "trad shape different from number of shells"
            self._t_radiative[
                self.v_boundary_inner_index
                + 1 : self.v_boundary_outer_index
                + 1
            ] = value
        else:
            raise ValueError(
                "Trying to set t_radiative for unmatching number" "of shells."
            )

    @property
    def radius(self):
        return self.time_explosion * self.velocity

    @property
    def r_inner(self):
        return self.model_state.geometry.r_inner

    @property
    def r_outer(self):
        return self.model_state.geometry.r_outer

    @property
    def r_middle(self):
        return 0.5 * self.r_inner + 0.5 * self.r_outer

    @property
    def velocity(self):

        if self._velocity is None:
            self._velocity = self.raw_velocity[
                self.v_boundary_inner_index : self.v_boundary_outer_index + 1
            ]
            self._velocity[0] = self.v_boundary_inner
            self._velocity[-1] = self.v_boundary_outer
        return self._velocity

    @property
    def v_inner(self):
        return self.model_state.geometry.v_inner

    @property
    def v_outer(self):
        return self.model_state.geometry.v_outer

    @property
    def v_middle(self):
        return 0.5 * self.v_inner + 0.5 * self.v_outer

    @property
    def density(self):
        return self.model_state.composition.density

    @property
    def abundance(self):
        if not self.raw_isotope_abundance.empty:
            self._abundance = self.raw_isotope_abundance.decay(
                self.time_explosion
            ).merge(self.raw_abundance)
        abundance = self._abundance.loc[
            :, self.v_boundary_inner_index : self.v_boundary_outer_index - 1
        ]
        abundance.columns = range(len(abundance.columns))
        return abundance

    @property
    def volume(self):
        return ((4.0 / 3) * np.pi * (self.r_outer**3 - self.r_inner**3)).cgs

    @property
    def no_of_shells(self):
        return len(self.velocity) - 1

    @property
    def no_of_raw_shells(self):
        return len(self.raw_velocity) - 1

    @property
    def v_boundary_inner(self):
        if self._v_boundary_inner is None:
            return self.raw_velocity[0]
        if self._v_boundary_inner < 0 * u.km / u.s:
            return self.raw_velocity[0]
        return self._v_boundary_inner

    @v_boundary_inner.setter
    def v_boundary_inner(self, value):
        if value is not None:
            if value > 0 * u.km / u.s:
                value = u.Quantity(value, self.v_boundary_inner.unit)
                if value > self.v_boundary_outer:
                    raise ValueError(
                        f"v_boundary_inner ({value}) must not be higher than "
                        f"v_boundary_outer ({self.v_boundary_outer})."
                    )
                if value > self.raw_velocity[-1]:
                    raise ValueError(
                        f"v_boundary_inner ({value}) is outside of the model range ({self.raw_velocity[-1]})."
                    )
                if value < self.raw_velocity[0]:
                    raise ValueError(
                        f"v_boundary_inner ({value}) is lower than the lowest shell ({self.raw_velocity[0]}) in the model."
                    )
        self._v_boundary_inner = value
        # Invalidate the cached cut-down velocity array
        self._velocity = None

    @property
    def v_boundary_outer(self):
        if self._v_boundary_outer is None:
            return self.raw_velocity[-1]
        if self._v_boundary_outer < 0 * u.km / u.s:
            return self.raw_velocity[-1]
        return self._v_boundary_outer

    @v_boundary_outer.setter
    def v_boundary_outer(self, value):
        if value is not None:
            if value > 0 * u.km / u.s:
                value = u.Quantity(value, self.v_boundary_outer.unit)
                if value < self.v_boundary_inner:
                    raise ValueError(
                        f"v_boundary_outer ({value}) must not be smaller than v_boundary_inner ({self.v_boundary_inner})."
                    )
                if value < self.raw_velocity[0]:
                    raise ValueError(
                        f"v_boundary_outer ({value}) is outside of the model range ({self.raw_velocity[0]})."
                    )
                if value > self.raw_velocity[-1]:
                    raise ValueError(
                        f"v_boundary_outer ({value}) is larger than the largest shell in the model ({self.raw_velocity[-1]})."
                    )
        self._v_boundary_outer = value
        # Invalidate the cached cut-down velocity array
        self._velocity = None

    @property
    def v_boundary_inner_index(self):
        if self.v_boundary_inner in self.raw_velocity:
            v_inner_ind = np.argwhere(
                self.raw_velocity == self.v_boundary_inner
            )[0][0]
        else:
            v_inner_ind = (
                np.searchsorted(self.raw_velocity, self.v_boundary_inner) - 1
            )
        return v_inner_ind

    @property
    def v_boundary_outer_index(self):
        if self.v_boundary_outer in self.raw_velocity:
            v_outer_ind = np.argwhere(
                self.raw_velocity == self.v_boundary_outer
            )[0][0]
        else:
            v_outer_ind = np.searchsorted(
                self.raw_velocity, self.v_boundary_outer
            )
        return v_outer_ind

    #    @property
    #    def v_boundary_inner_index(self):
    #        if self.v_boundary_inner <= self.raw_velocity[0]:
    #            return 0
    #        else:
    #            idx = max(0,
    #                      self.raw_velocity.searchsorted(self.v_boundary_inner) - 1)
    #            # check for zero volume of designated first cell
    #            if np.isclose(self.v_boundary_inner, self.raw_velocity[idx + 1],
    #                          atol=1e-8 * u.km / u.s) and (self.v_boundary_inner <=
    #                                                           self.raw_velocity[idx + 1]):
    #                idx += 1
    #            return idx
    #
    #    @property
    #    def v_boundary_outer_index(self):
    #        if self.v_boundary_outer >= self.raw_velocity[-1]:
    #            return None
    #        return self.raw_velocity.searchsorted(self.v_boundary_outer) + 1

    @classmethod
    def from_config(cls, config, atom_data=None):
        """
        Create a new Radial1DModel instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        atom_data : tardis.io.AtomData

        Returns
        -------
        Radial1DModel
        """
        time_explosion = config.supernova.time_explosion.cgs

        structure = config.model.structure
        electron_densities = None
        temperature = None
        if structure.type == "specific":
            velocity = quantity_linspace(
                structure.velocity.start,
                structure.velocity.stop,
                structure.velocity.num + 1,
            ).cgs
            homologous_density = HomologousDensity.from_config(config)
        elif structure.type == "file":
            if os.path.isabs(structure.filename):
                structure_fname = structure.filename
            else:
                structure_fname = os.path.join(
                    config.config_dirname, structure.filename
                )

            (
                time_0,
                velocity,
                density_0,
                electron_densities,
                temperature,
            ) = read_density_file(structure_fname, structure.filetype)
            density_0 = density_0.insert(0, 0)
            homologous_density = HomologousDensity(density_0, time_0)
        else:
            raise NotImplementedError
        # Note: This is the number of shells *without* taking in mind the
        #       v boundaries.
        no_of_shells = len(velocity) - 1

        if temperature:
            t_radiative = temperature
        elif config.plasma.initial_t_rad > 0 * u.K:
            t_radiative = (
                np.ones(no_of_shells + 1) * config.plasma.initial_t_rad
            )
        else:
            t_radiative = None

        if config.plasma.initial_t_inner < 0.0 * u.K:
            luminosity_requested = config.supernova.luminosity_requested
            t_inner = None
        else:
            luminosity_requested = None
            t_inner = config.plasma.initial_t_inner

        abundances_section = config.model.abundances
        isotope_abundance = pd.DataFrame()

        if abundances_section.type == "uniform":
            abundance, isotope_abundance = read_uniform_abundances(
                abundances_section, no_of_shells
            )

        elif abundances_section.type == "file":
            if os.path.isabs(abundances_section.filename):
                abundances_fname = abundances_section.filename
            else:
                abundances_fname = os.path.join(
                    config.config_dirname, abundances_section.filename
                )

            index, abundance, isotope_abundance = read_abundances_file(
                abundances_fname, abundances_section.filetype
            )

        abundance = abundance.replace(np.nan, 0.0)
        abundance = abundance[abundance.sum(axis=1) > 0]

        norm_factor = abundance.sum(axis=0) + isotope_abundance.sum(axis=0)

        if np.any(np.abs(norm_factor - 1) > 1e-12):
            logger.warning(
                "Abundances have not been normalized to 1." " - normalizing"
            )
            abundance /= norm_factor
            isotope_abundance /= norm_factor

        isotope_abundance = IsotopeAbundances(isotope_abundance)

        elemental_mass = None
        if atom_data is not None:
            elemental_mass = atom_data.atom_data.mass

        return cls(
            velocity=velocity,
            homologous_density=homologous_density,
            abundance=abundance,
            isotope_abundance=isotope_abundance,
            time_explosion=time_explosion,
            t_radiative=t_radiative,
            t_inner=t_inner,
            elemental_mass=elemental_mass,
            luminosity_requested=luminosity_requested,
            dilution_factor=None,
            v_boundary_inner=structure.get("v_inner_boundary", None),
            v_boundary_outer=structure.get("v_outer_boundary", None),
            electron_densities=electron_densities,
        )

    @classmethod
    def from_csvy(cls, config, atom_data=None):
        """
        Create a new Radial1DModel instance from a Configuration object.

        Parameters
        ----------
        config : tardis.io.config_reader.Configuration
        atom_data : tardis.io.AtomData

        Returns
        -------
        Radial1DModel
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
        base_dir = os.path.abspath(os.path.dirname(__file__))
        schema_dir = os.path.join(base_dir, "..", "io", "schemas")
        csvy_schema_file = os.path.join(schema_dir, "csvy_model.yml")
        csvy_model_config = Configuration(
            validate_dict(csvy_model_config, schemapath=csvy_schema_file)
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
        temperature = None

        # if hasattr(csvy_model_config, 'v_inner_boundary'):
        #    v_boundary_inner = csvy_model_config.v_inner_boundary
        # else:
        #    v_boundary_inner = None

        # if hasattr(csvy_model_config, 'v_outer_boundary'):
        #    v_boundary_outer = csvy_model_config.v_outer_boundary
        # else:
        #    v_boundary_outer = None

        if hasattr(config, "model"):
            if hasattr(config.model, "v_inner_boundary"):
                v_boundary_inner = config.model.v_inner_boundary
            else:
                v_boundary_inner = None

            if hasattr(config.model, "v_outer_boundary"):
                v_boundary_outer = config.model.v_outer_boundary
            else:
                v_boundary_outer = None
        else:
            v_boundary_inner = None
            v_boundary_outer = None

        if hasattr(csvy_model_config, "velocity"):
            velocity = quantity_linspace(
                csvy_model_config.velocity.start,
                csvy_model_config.velocity.stop,
                csvy_model_config.velocity.num + 1,
            ).cgs
        else:
            velocity_field_index = [
                field["name"] for field in csvy_model_config.datatype.fields
            ].index("velocity")
            velocity_unit = u.Unit(
                csvy_model_config.datatype.fields[velocity_field_index]["unit"]
            )
            velocity = csvy_model_data["velocity"].values * velocity_unit
            velocity = velocity.to("cm/s")

        if hasattr(csvy_model_config, "density"):
            homologous_density = HomologousDensity.from_csvy(
                config, csvy_model_config
            )
        else:
            time_0 = csvy_model_config.model_density_time_0
            density_field_index = [
                field["name"] for field in csvy_model_config.datatype.fields
            ].index("density")
            density_unit = u.Unit(
                csvy_model_config.datatype.fields[density_field_index]["unit"]
            )
            density_0 = csvy_model_data["density"].values * density_unit
            density_0 = density_0.to("g/cm^3")[1:]
            density_0 = density_0.insert(0, 0)
            homologous_density = HomologousDensity(density_0, time_0)

        no_of_shells = len(velocity) - 1

        # TODO -- implement t_radiative
        # t_radiative = None
        if temperature:
            t_radiative = temperature
        elif hasattr(csvy_model_data, "columns"):
            if "t_rad" in csvy_model_data.columns:
                t_rad_field_index = [
                    field["name"] for field in csvy_model_config.datatype.fields
                ].index("t_rad")
                t_rad_unit = u.Unit(
                    csvy_model_config.datatype.fields[t_rad_field_index]["unit"]
                )
                t_radiative = (
                    csvy_model_data["t_rad"].iloc[0:].values * t_rad_unit
                )
            else:
                t_radiative = None

        dilution_factor = None
        if hasattr(csvy_model_data, "columns"):
            if "dilution_factor" in csvy_model_data.columns:
                dilution_factor = (
                    csvy_model_data["dilution_factor"].iloc[0:].to_numpy()
                )

        elif config.plasma.initial_t_rad > 0 * u.K:
            t_radiative = np.ones(no_of_shells) * config.plasma.initial_t_rad
        else:
            t_radiative = None

        if config.plasma.initial_t_inner < 0.0 * u.K:
            luminosity_requested = config.supernova.luminosity_requested
            t_inner = None
        else:
            luminosity_requested = None
            t_inner = config.plasma.initial_t_inner

        if hasattr(csvy_model_config, "abundance"):
            abundances_section = csvy_model_config.abundance
            abundance, isotope_abundance = read_uniform_abundances(
                abundances_section, no_of_shells
            )
        else:
            index, abundance, isotope_abundance = parse_csv_abundances(
                csvy_model_data
            )
            abundance = abundance.loc[:, 1:]
            abundance.columns = np.arange(abundance.shape[1])
            isotope_abundance = isotope_abundance.loc[:, 1:]
            isotope_abundance.columns = np.arange(isotope_abundance.shape[1])

        abundance = abundance.replace(np.nan, 0.0)
        abundance = abundance[abundance.sum(axis=1) > 0]
        isotope_abundance = isotope_abundance.replace(np.nan, 0.0)
        isotope_abundance = isotope_abundance[isotope_abundance.sum(axis=1) > 0]
        norm_factor = abundance.sum(axis=0) + isotope_abundance.sum(axis=0)

        if np.any(np.abs(norm_factor - 1) > 1e-12):
            logger.warning(
                "Abundances have not been normalized to 1." " - normalizing"
            )
            abundance /= norm_factor
            isotope_abundance /= norm_factor

        # isotope_abundance = IsotopeAbundances(isotope_abundance)
        isotope_abundance = IsotopeAbundances(
            isotope_abundance, time_0=csvy_model_config.model_isotope_time_0
        )
        # isotope_abundance.time_0 = csvy_model_config.model_isotope_time_0

        elemental_mass = None
        if atom_data is not None:
            elemental_mass = atom_data.atom_data.mass

        return cls(
            velocity=velocity,
            homologous_density=homologous_density,
            abundance=abundance,
            isotope_abundance=isotope_abundance,
            time_explosion=time_explosion,
            t_radiative=t_radiative,
            t_inner=t_inner,
            elemental_mass=elemental_mass,
            luminosity_requested=luminosity_requested,
            dilution_factor=dilution_factor,
            v_boundary_inner=v_boundary_inner,
            v_boundary_outer=v_boundary_outer,
            electron_densities=electron_densities,
        )
