# atomic model


import logging
import numpy as np
import pandas as pd

from scipy import interpolate
from collections import OrderedDict
from astropy import units as u
from tardis import constants as const
from astropy.units import Quantity
from tardis.io.atom_data.util import resolve_atom_data_fname


class AtomDataNotPreparedError(Exception):
    pass


class AtomDataMissingError(Exception):
    pass


logger = logging.getLogger(__name__)


class AtomData(object):
    """
    Class for storing atomic data

    Parameters
    ----------
    atom_data: pandas.DataFrame
        A DataFrame containing the *basic atomic data* with:
            index: atomic_number;
            columns: symbol, name, mass[u].

    ionization_data: pandas.DataFrame
        A DataFrame containing the *ionization data* with:
            index: atomic_number, ion_number;
            columns: ionization_energy[eV].

       It is important to note here is that `ion_number` describes the *final ion state*
            e.g. H I - H II is described with ion=1

    levels: pandas.DataFrame
        A DataFrame containing the *levels data* with:
            index: numerical index;
            columns: atomic_number, ion_number, level_number, energy[eV], g[1], metastable.

    lines: pandas.DataFrame
        A DataFrame containing the *lines data* with:
            index: numerical index;
            columns: line_id, atomic_number, ion_number, level_number_lower, level_number_upper,
                wavelength[angstrom], nu[Hz], f_lu[1], f_ul[1], B_ul[?], B_ul[?], A_ul[1/s].

    macro_atom_data:
        A DataFrame containing the *macro atom data* with:
            index: numerical index;
            columns: atomic_number, ion_number, source_level_number, destination_level_number,
                transition_line_id, transition_type, transition_probability;

    macro_atom_references:
        A DataFrame containing  the *macro atom references* with:
            index: numerical index;
            columns: atomic_number, ion_number, source_level_number, count_down, count_up, count_total.

        Refer to the docs: http://tardis.readthedocs.io/en/latest/physics/plasma/macroatom.html

    collision_data: (pandas.DataFrame, np.array)
        A DataFrame containing the *electron collisions data* with:
            index: atomic_number, ion_number, level_number_lower, level_number_upper;
            columns: e_col_id, delta_e, g_ratio, c_ul;

    collision_data_temperatures: np.array
        An array with the collision temperatures.

    zeta_data:
        A DataFrame containing the *zeta data* for the
        nebular ionization calculation
        (i.e., the fraction of recombinations that go directly to the
        ground state) with:
            index: atomic_number, ion_charge;
            columns: temperatures[K]

    synpp_refs: ?

    photoionization_data: pandas.DataFrame
        A DataFrame containing the *photoionization data* with:
            index: numerical index;
            columns: atomic_number, ion_number, level_number, nu[Hz], x_sect[cm^2]

    Attributes
    -------------
    prepared: bool

    atom_data: pandas.DataFrame
    ionization_data: pandas.DataFrame
    macro_atom_data_all: pandas.DataFrame
    macro_atom_references_all: pandas.DataFrame
    collision_data: pandas.DataFrame
    collision_data_temperatures: numpy.array
    zeta_data: pandas.DataFrame
    synpp_refs: pandas.DataFrame
    symbol2atomic_number: OrderedDict
    atomic_number2symbol OrderedDict
    photoionization_data: pandas.DataFrame

    Methods
    --------
    from_hdf
    prepare_atom_data

    Notes
    ------
    1. The units of some columns are given in the square brackets. They are **NOT** the parts of columns' names!

    """

    hdf_names = [
        "atom_data",
        "ionization_data",
        "levels",
        "lines",
        "macro_atom_data",
        "macro_atom_references",
        "zeta_data",
        "collision_data",
        "collision_data_temperatures",
        "synpp_refs",
        "photoionization_data",
    ]

    # List of tuples of the related dataframes.
    # Either all or none of the related dataframes must be given
    related_groups = [
        ("macro_atom_data_all", "macro_atom_references_all"),
        ("collision_data", "collision_data_temperatures"),
    ]

    @classmethod
    def from_hdf(cls, fname=None):
        """
        Function to read the atom data from a TARDIS atom HDF Store

        Parameters
        ----------

        fname: str, optional
            Path to the HDFStore file or name of known atom data file
            (default: None)
        """

        dataframes = dict()
        nonavailable = list()

        fname = resolve_atom_data_fname(fname)

        with pd.HDFStore(fname, "r") as store:
            for name in cls.hdf_names:
                try:
                    dataframes[name] = store[name]
                except KeyError:
                    nonavailable.append(name)

            atom_data = cls(**dataframes)

            try:
                atom_data.uuid1 = store.root._v_attrs["uuid1"].decode("ascii")
            except KeyError:
                atom_data.uuid1 = None

            try:
                atom_data.md5 = store.root._v_attrs["md5"].decode("ascii")
            except KeyError:
                atom_data.md5 = None

            try:
                atom_data.version = store.root._v_attrs["database_version"]
            except KeyError:
                atom_data.version = None

            # ToDo: strore data sources as attributes in carsus

            logger.info(
                "Read Atom Data with UUID={0} and MD5={1}.".format(
                    atom_data.uuid1, atom_data.md5
                )
            )
            if nonavailable:
                logger.info(
                    "Non provided atomic data: {0}".format(
                        ", ".join(nonavailable)
                    )
                )

        return atom_data

    def __init__(
        self,
        atom_data,
        ionization_data,
        levels=None,
        lines=None,
        macro_atom_data=None,
        macro_atom_references=None,
        zeta_data=None,
        collision_data=None,
        collision_data_temperatures=None,
        synpp_refs=None,
        photoionization_data=None,
    ):

        self.prepared = False

        # CONVERT VALUES TO CGS UNITS

        # Convert atomic masses to CGS
        # We have to use constants.u because astropy uses
        # different values for the unit u and the constant.
        # This is changed in later versions of astropy (
        # the value of constants.u is used in all cases)
        if u.u.cgs == const.u.cgs:
            atom_data.loc[:, "mass"] = Quantity(
                atom_data["mass"].values, "u"
            ).cgs
        else:
            atom_data.loc[:, "mass"] = atom_data["mass"].values * const.u.cgs

        # Convert ionization energies to CGS
        ionization_data = ionization_data.squeeze()
        ionization_data[:] = Quantity(ionization_data[:], "eV").cgs

        # Convert energy to CGS
        levels.loc[:, "energy"] = Quantity(levels["energy"].values, "eV").cgs

        # Create a new columns with wavelengths in the CGS units
        lines.loc[:, "wavelength_cm"] = Quantity(
            lines["wavelength"], "angstrom"
        ).cgs

        # SET ATTRIBUTES

        self.atom_data = atom_data
        self.ionization_data = ionization_data
        self.levels = levels
        self.lines = lines

        # Rename these (drop "_all") when `prepare_atom_data` is removed!
        self.macro_atom_data_all = macro_atom_data
        self.macro_atom_references_all = macro_atom_references

        self.zeta_data = zeta_data

        self.collision_data = collision_data
        self.collision_data_temperatures = collision_data_temperatures

        self.synpp_refs = synpp_refs

        self.photoionization_data = photoionization_data

        self._check_related()

        self.symbol2atomic_number = OrderedDict(
            zip(self.atom_data["symbol"].values, self.atom_data.index)
        )
        self.atomic_number2symbol = OrderedDict(
            zip(self.atom_data.index, self.atom_data["symbol"])
        )

    def _check_related(self):
        """
        Check that either all or none of the related dataframes are given.
        """
        for group in self.related_groups:
            check_list = [name for name in group if getattr(self, name) is None]

            if len(check_list) != 0 and len(check_list) != len(group):
                raise AtomDataMissingError(
                    "The following dataframes from the related group [{0}] "
                    "were not given: {1}".format(
                        ", ".join(group), ", ".join(check_list)
                    )
                )

    def prepare_atom_data(
        self,
        selected_atomic_numbers,
        line_interaction_type="scatter",
        nlte_species=[],
    ):
        """
        Prepares the atom data to set the lines, levels and if requested macro
        atom data.  This function mainly cuts the `levels` and `lines` by
        discarding any data that is not needed (any data for atoms that are not
        needed

        Parameters
        ----------

        selected_atoms : `~set`
            set of selected atom numbers, e.g. set([14, 26])

        line_interaction_type : `~str`
            can be 'scatter', 'downbranch' or 'macroatom'

        """
        if not self.prepared:
            self.prepared = True
        else:
            raise AtomDataNotPreparedError("AtomData was already prepared")
        self.selected_atomic_numbers = selected_atomic_numbers

        self._check_selected_atomic_numbers()

        self.nlte_species = nlte_species

        self.levels = self.levels[
            self.levels.index.isin(
                self.selected_atomic_numbers, level="atomic_number"
            )
        ]

        self.levels_index = pd.Series(
            np.arange(len(self.levels), dtype=int), index=self.levels.index
        )

        # cutting levels_lines
        self.lines = self.lines[
            self.lines.index.isin(
                self.selected_atomic_numbers, level="atomic_number"
            )
        ]

        self.lines.sort_values(by="wavelength", inplace=True)

        self.lines_index = pd.Series(
            np.arange(len(self.lines), dtype=int),
            index=self.lines.set_index("line_id").index,
        )

        tmp_lines_lower2level_idx = self.lines.index.droplevel(
            "level_number_upper"
        )

        self.lines_lower2level_idx = (
            self.levels_index.loc[tmp_lines_lower2level_idx]
            .astype(np.int64)
            .values
        )

        tmp_lines_upper2level_idx = self.lines.index.droplevel(
            "level_number_lower"
        )

        self.lines_upper2level_idx = (
            self.levels_index.loc[tmp_lines_upper2level_idx]
            .astype(np.int64)
            .values
        )

        if (
            self.macro_atom_data_all is not None
            and not line_interaction_type == "scatter"
        ):

            self.macro_atom_data = self.macro_atom_data_all.loc[
                self.macro_atom_data_all["atomic_number"].isin(
                    self.selected_atomic_numbers
                )
            ].copy()

            self.macro_atom_references = self.macro_atom_references_all[
                self.macro_atom_references_all.index.isin(
                    self.selected_atomic_numbers, level="atomic_number"
                )
            ].copy()

            if line_interaction_type == "downbranch":
                self.macro_atom_data = self.macro_atom_data.loc[
                    self.macro_atom_data["transition_type"] == -1
                ]
                self.macro_atom_references = self.macro_atom_references.loc[
                    self.macro_atom_references["count_down"] > 0
                ]
                self.macro_atom_references.loc[
                    :, "count_total"
                ] = self.macro_atom_references["count_down"]
                self.macro_atom_references.loc[
                    :, "block_references"
                ] = np.hstack(
                    (
                        0,
                        np.cumsum(
                            self.macro_atom_references["count_down"].values[:-1]
                        ),
                    )
                )

            elif line_interaction_type == "macroatom":
                self.macro_atom_references.loc[
                    :, "block_references"
                ] = np.hstack(
                    (
                        0,
                        np.cumsum(
                            self.macro_atom_references["count_total"].values[
                                :-1
                            ]
                        ),
                    )
                )

            self.macro_atom_references.loc[:, "references_idx"] = np.arange(
                len(self.macro_atom_references)
            )

            self.macro_atom_data.loc[:, "lines_idx"] = self.lines_index.loc[
                self.macro_atom_data["transition_line_id"]
            ].values

            self.lines_upper2macro_reference_idx = (
                self.macro_atom_references.loc[
                    tmp_lines_upper2level_idx, "references_idx"
                ]
                .astype(np.int64)
                .values
            )

            if line_interaction_type == "macroatom":
                # Sets all
                tmp_macro_destination_level_idx = pd.MultiIndex.from_arrays(
                    [
                        self.macro_atom_data["atomic_number"],
                        self.macro_atom_data["ion_number"],
                        self.macro_atom_data["destination_level_number"],
                    ]
                )

                tmp_macro_source_level_idx = pd.MultiIndex.from_arrays(
                    [
                        self.macro_atom_data["atomic_number"],
                        self.macro_atom_data["ion_number"],
                        self.macro_atom_data["source_level_number"],
                    ]
                )

                self.macro_atom_data.loc[:, "destination_level_idx"] = (
                    self.macro_atom_references.loc[
                        tmp_macro_destination_level_idx, "references_idx"
                    ]
                    .astype(np.int64)
                    .values
                )

                self.macro_atom_data.loc[:, "source_level_idx"] = (
                    self.macro_atom_references.loc[
                        tmp_macro_source_level_idx, "references_idx"
                    ]
                    .astype(np.int64)
                    .values
                )

            elif line_interaction_type == "downbranch":
                # Sets all the destination levels to -1 to indicate that they
                # are not used in downbranch calculations
                self.macro_atom_data.loc[:, "destination_level_idx"] = -1

        self.nlte_data = NLTEData(self, nlte_species)

    def _check_selected_atomic_numbers(self):
        selected_atomic_numbers = self.selected_atomic_numbers
        available_atomic_numbers = np.unique(
            self.ionization_data.index.get_level_values(0)
        )
        atomic_number_check = np.isin(
            selected_atomic_numbers, available_atomic_numbers
        )

        if not all(atomic_number_check):
            missing_atom_mask = np.logical_not(atomic_number_check)
            missing_atomic_numbers = selected_atomic_numbers[missing_atom_mask]
            missing_numbers_str = ",".join(missing_atomic_numbers.astype("str"))
            msg = "For atomic numbers {} there is no atomic data.".format(
                missing_numbers_str
            )
            raise AtomDataMissingError(msg)

    def __repr__(self):
        return "<Atomic Data UUID={} MD5={} Lines={:d} Levels={:d}>".format(
            self.uuid1,
            self.md5,
            self.lines.line_id.count(),
            self.levels.energy.count(),
        )


class NLTEData(object):
    def __init__(self, atom_data, nlte_species):
        self.atom_data = atom_data
        self.lines = atom_data.lines.reset_index()
        self.nlte_species = nlte_species

        if nlte_species:
            logger.info("Preparing the NLTE data")
            self._init_indices()
            if atom_data.collision_data is not None:
                self._create_collision_coefficient_matrix()

    def _init_indices(self):
        self.lines_idx = {}
        self.lines_level_number_lower = {}
        self.lines_level_number_upper = {}
        self.A_uls = {}
        self.B_uls = {}
        self.B_lus = {}

        for species in self.nlte_species:
            lines_idx = np.where(
                (self.lines.atomic_number == species[0])
                & (self.lines.ion_number == species[1])
            )
            self.lines_idx[species] = lines_idx
            self.lines_level_number_lower[
                species
            ] = self.lines.level_number_lower.values[lines_idx].astype(int)
            self.lines_level_number_upper[
                species
            ] = self.lines.level_number_upper.values[lines_idx].astype(int)

            self.A_uls[species] = self.atom_data.lines.A_ul.values[lines_idx]
            self.B_uls[species] = self.atom_data.lines.B_ul.values[lines_idx]
            self.B_lus[species] = self.atom_data.lines.B_lu.values[lines_idx]

    def _create_collision_coefficient_matrix(self):
        self.C_ul_interpolator = {}
        self.delta_E_matrices = {}
        self.g_ratio_matrices = {}
        collision_group = self.atom_data.collision_data.groupby(
            level=["atomic_number", "ion_number"]
        )
        for species in self.nlte_species:
            no_of_levels = self.atom_data.levels.loc[species].energy.count()
            C_ul_matrix = np.zeros(
                (
                    no_of_levels,
                    no_of_levels,
                    len(self.atom_data.collision_data_temperatures),
                )
            )
            delta_E_matrix = np.zeros((no_of_levels, no_of_levels))
            g_ratio_matrix = np.zeros((no_of_levels, no_of_levels))

            for (
                (
                    atomic_number,
                    ion_number,
                    level_number_lower,
                    level_number_upper,
                ),
                line,
            ) in collision_group.get_group(species).iterrows():
                # line.columns : delta_e, g_ratio, temperatures ...
                C_ul_matrix[
                    level_number_lower, level_number_upper, :
                ] = line.values[2:]
                delta_E_matrix[level_number_lower, level_number_upper] = line[
                    "delta_e"
                ]
                # TODO TARDISATOMIC fix change the g_ratio to be the otherway round - I flip them now here.
                g_ratio_matrix[level_number_lower, level_number_upper] = (
                    1 / line["g_ratio"]
                )
            self.C_ul_interpolator[species] = interpolate.interp1d(
                self.atom_data.collision_data_temperatures, C_ul_matrix
            )
            self.delta_E_matrices[species] = delta_E_matrix

            self.g_ratio_matrices[species] = g_ratio_matrix

    def get_collision_matrix(self, species, t_electrons):
        """
        Creat collision matrix by interpolating the C_ul values for
        the desired temperatures.
        """
        c_ul_matrix = self.C_ul_interpolator[species](t_electrons)
        no_of_levels = c_ul_matrix.shape[0]
        c_ul_matrix[np.isnan(c_ul_matrix)] = 0.0

        # TODO in tardisatomic the g_ratio is the other way round - here I'll flip it in prepare_collision matrix

        c_lu_matrix = (
            c_ul_matrix
            * np.exp(
                -self.delta_E_matrices[species].reshape(
                    (no_of_levels, no_of_levels, 1)
                )
                / t_electrons.reshape((1, 1, t_electrons.shape[0]))
            )
            * self.g_ratio_matrices[species].reshape(
                (no_of_levels, no_of_levels, 1)
            )
        )
        return c_ul_matrix + c_lu_matrix.transpose(1, 0, 2)
