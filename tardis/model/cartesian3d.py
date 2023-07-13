from tardis.model.geometry import Cartesian3DGeometry

# packet sources need an origin


class Cartesian3DModelSolver(HDFWriterMixin):
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
        model_state,
        # velocity,
        # homologous_density,           # should go into modelstate
        # abundance,                    # should go into composition
        # isotope_abundance,            # should go into composition
        time_explosion,
        # t_inner,
        # elemental_mass,               # should go into composition
        # luminosity_requested=None,    # should go into radiation field
        # t_radiative=None,             # should go into radiation field
        # dilution_factor=None,         # should go into radiation field
        # v_boundary_inner=None,
        # v_boundary_outer=None,
        # electron_densities=None,      # should go into plasma
    ):
        self.initial_model_state = model_state

    def evolve(self, time_explosion):
        # TODO: this method is overloaded by users to return the modelstate

        self.model_state = ModelState(
            composition=composition,
            geometry=geometry,
            time_explosion=time_explosion,
        )

        if t_inner is None:
            if luminosity_requested is not None:
                self.t_inner = (
                    (
                        luminosity_requested
                        / (
                            4
                            * np.pi
                            * self.r_inner[0]
                            ** 2  # does not make sense! need to define inner boundary velocity/radius instead
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
            self._t_radiative = t_radiative

        if dilution_factor is None:
            self._dilution_factor = 0.5 * (
                1
                - np.sqrt(
                    1 - (self.r_inner[0] ** 2 / self.r_middle**2).to(1).value
                )
            )
        else:
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
        if len(self._dilution_factor) == self.no_of_cells:
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
        else:
            raise ValueError(
                "Trying to set dilution_factor for unmatching number"
                "of shells."
            )

    @property
    def t_radiative(self):
        if len(self._t_radiative) == self.no_of_cells:
            return self._t_radiative

    @t_radiative.setter
    def t_radiative(self, value):
        if len(value) == len(self._t_radiative):
            self._t_radiative = value
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
        raise NotImplementedError

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
    def no_of_cells(self):
        return len(self.velocity) - 1


# this has to do a lot of work!
# @classmethod
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
        t_radiative = np.ones(no_of_shells + 1) * config.plasma.initial_t_rad
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
