import logging

import numpy as np
import pandas as pd

from tardis import constants as const
from tardis.io.atom_data.collision_data import CollisionData
from tardis.io.atom_data.decay_radiation_data import DecayRadiationData
from tardis.io.atom_data.ionization_data import (
    IonizationData,
    PhotoIonizationData,
    ZetaData,
)
from tardis.io.atom_data.levels_data import LevelsData
from tardis.io.atom_data.lines_data import LineList, LinesData
from tardis.io.atom_data.macro_atom_data import MacroAtomData
from tardis.io.atom_data.nlte_data import NLTEData
from tardis.io.atom_data.util import (
    resolve_atom_data_fname,
    set_atom_data_attributes,
)
from tardis.plasma.properties.continuum_processes import (
    get_ground_state_multi_index,
)


class AtomDataNotPreparedError(Exception):
    pass


class AtomDataMissingError(Exception):
    pass


logger = logging.getLogger(__name__)


class EnergyData:
    def active_data(self):
        pass


class AtomData:
    """
    Class for storing atomic data

    Parameters
    ----------
    atom_data : pandas.DataFrame
    A DataFrame containing the *basic atomic data* with:
        index : atomic_number
        columns : symbol, name, mass[u].

    ionization_data : pandas.DataFrame
    A DataFrame containing the *ionization data* with:
        index : atomic_number, ion_number
        columns : ionization_energy[eV].
    It is important to note here is that `ion_number` describes the *final ion state*
    e.g. H I - H II is described with ion=1

    levels : pandas.DataFrame
    A DataFrame containing the *levels data* with:
        index : numerical index
        columns : atomic_number, ion_number, level_number, energy[eV], g[1], metastable.

    lines : pandas.DataFrame
    A DataFrame containing the *lines data* with:
        index : numerical index
        columns : line_id, atomic_number, ion_number, level_number_lower, level_number_upper,
        wavelength[angstrom], nu[Hz], f_lu[1], f_ul[1], B_ul[?], B_ul[?], A_ul[1/s].

    macro_atom_data :
    A DataFrame containing the *macro atom data* with:
        index : numerical index
        columns : atomic_number, ion_number, source_level_number, destination_level_number,
        transition_line_id, transition_type, transition_probability;

    macro_atom_references :
    A DataFrame containing  the *macro atom references* with:
        index : numerical index
        columns : atomic_number, ion_number, source_level_number, count_down, count_up, count_total.
    Refer to the docs: http://tardis.readthedocs.io/en/latest/physics/plasma/macroatom.html

    collision_data : (pandas.DataFrame, np.array)
    A DataFrame containing the *electron collisions data* with:
        index : atomic_number, ion_number, level_number_lower, level_number_upper
        columns : e_col_id, delta_e, g_ratio, c_ul;

    collision_data_temperatures : np.array
        An array with the collision temperatures.

    zeta_data :
    A DataFrame containing the *zeta data* for the
    nebular ionization calculation
    (i.e., the fraction of recombinations that go directly to the
    ground state) with:
        index : atomic_number, ion_charge
        columns : temperatures[K]

    synpp_refs : ?

    photoionization_data : pandas.DataFrame
    A DataFrame containing the *photoionization data* with:
        index : numerical index
        columns : atomic_number, ion_number, level_number, nu[Hz], x_sect[cm^2]

    two_photon_data : pandas.DataFrame
    A DataFrame containing the *two photon decay data* with:
        index: atomic_number, ion_number, level_number_lower, level_number_upper
        columns: A_ul[1/s], nu0[Hz], alpha, beta, gamma

    decay_radiation_data : pandas.DataFrame
    A dataframe containing the *decay radiation data* with:
        index: Isotope names
        columns: atomic_number, element, Rad energy, Rad intensity decay mode.
        Curated from nndc

    Attributes
    ----------
    prepared : bool
    atom_data : pandas.DataFrame
    ionization_data : pandas.DataFrame
    macro_atom_data_all : pandas.DataFrame
    macro_atom_references_all : pandas.DataFrame
    collision_data : pandas.DataFrame
    collision_data_temperatures : numpy.array
    zeta_data : pandas.DataFrame
    synpp_refs : pandas.DataFrame
    photoionization_data : pandas.DataFrame
    two_photon_data : pandas.DataFrame
    decay_radiation_data : pandas.DataFrame

    Methods
    -------
    from_hdf:
        Function to read the atom data from a TARDIS atom HDF Store
    prepare_atom_data:
        Prepares the atom data to set the lines, levels and if requested macro
        atom data.  This function mainly cuts the `levels` and `lines` by
        discarding any data that is not needed (any data for atoms that are not
        needed

    Notes
    -----
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
        "yg_data",
        "two_photon_data",
        "linelist",
        "decay_radiation_data",
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
        fname : Path, optional
            Path to the HDFStore file or name of known atom data file
            (default: None)
        """
        dataframes = {}
        nonavailable = []

        fname = resolve_atom_data_fname(fname)

        with pd.HDFStore(fname, "r") as store:
            for name in cls.hdf_names:
                try:
                    dataframes[name] = store.select(name)
                except KeyError:
                    logger.debug(f"Dataframe does not contain {name} column")
                    nonavailable.append(name)

            if "metadata" in store:
                carsus_version_str = (
                    store["metadata"].loc[("format", "version")].value
                )
                carsus_version = tuple(map(int, carsus_version_str.split(".")))
                if carsus_version == (1, 0):
                    # Checks for various collisional data from Carsus files
                    if "collisions_data" in store:
                        try:
                            dataframes["collision_data_temperatures"] = store[
                                "collisions_metadata"
                            ].temperatures
                            if "cmfgen" in store["collisions_metadata"].dataset:
                                dataframes["yg_data"] = store["collisions_data"]
                                dataframes["collision_data"] = "dummy value"
                            elif (
                                "chianti"
                                in store["collisions_metadata"].dataset
                            ):
                                dataframes["collision_data"] = store[
                                    "collisions_data"
                                ]
                            else:
                                raise KeyError(
                                    "Atomic Data Collisions Not a Valid Chanti or CMFGEN Carsus Data File"
                                )
                        except KeyError as e:
                            logger.warning(
                                "Atomic Data is not a Valid Carsus Atomic Data File"
                            )
                            raise
                    dataframes["levels"] = store["levels_data"]
                    dataframes["lines"] = store["lines_data"]
                else:
                    raise ValueError(
                        f"Current carsus version, {carsus_version}, is not supported."
                    )
            if "linelist" in store:
                dataframes["linelist"] = store["linelist"]

            atom_data = cls(**dataframes)

            set_atom_data_attributes(atom_data, store, "uuid1")
            set_atom_data_attributes(atom_data, store, "md5")
            set_atom_data_attributes(atom_data, store, "database_version")

            # TODO: strore data sources as attributes in carsus

            logger.info(
                f"Reading Atom Data with: UUID = {atom_data.uuid1} MD5  = {atom_data.md5} "
            )
            if nonavailable:
                logger.info(
                    "Non provided Atomic Data: {0}".format(
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
        yg_data=None,
        two_photon_data=None,
        linelist=None,
        decay_radiation_data=None,
    ):
        self.prepared = False

        # CONVERT VALUES TO CGS UNITS

        # Convert atomic masses to CGS
        # We have to use constants.u because astropy uses
        # different values for the unit u and the constant.
        # This is changed in later versions of astropy (
        # the value of constants.u is used in all cases)
        atom_data.loc[:, "mass"] = atom_data["mass"].values * const.u.cgs.value

        # SET ATTRIBUTES

        self.atom_data = atom_data
        self.ionization_data = IonizationData(ionization_data).data
        self.levels = LevelsData(levels).data
        self.lines = LinesData(lines).data
        self.macro_atom_data_class = MacroAtomData(
            macro_atom_data, macro_atom_references
        )
        # replace this along with _check_related refactor
        self.macro_atom_data_all = (
            self.macro_atom_data_class.transition_probability_data
        )
        self.macro_atom_references_all = (
            self.macro_atom_data_class.block_reference_data
        )

        self.zeta_data = ZetaData(zeta_data).data

        collision_data_class = CollisionData(
            collision_data, collision_data_temperatures, yg_data
        )
        self.collision_data = collision_data_class.data
        self.collision_data_temperatures = collision_data_class.temperatures

        self.synpp_refs = synpp_refs

        self.photoionization_data = PhotoIonizationData(
            photoionization_data
        ).data

        self.yg_data = collision_data_class.yg

        self.two_photon_data = two_photon_data

        if linelist is not None:
            self.linelist = LineList(linelist).data

        if decay_radiation_data is not None:
            self.decay_radiation_data = DecayRadiationData(
                decay_radiation_data
            ).data
        self._check_related()

    def _check_related(self):
        """
        Check that either all or none of the related dataframes are given.
        """
        for group in self.related_groups:
            check_list = [name for name in group if getattr(self, name) is None]

            if len(check_list) != 0 and len(check_list) != len(group):
                raise AtomDataMissingError(
                    f'The following dataframes from the related group [{", ".join(group)}]'
                    f'were not given: {", ".join(check_list)}'
                )

    def prepare_atom_data(
        self,
        selected_atomic_numbers,
        line_interaction_type,
        nlte_species,
        continuum_interaction_species,
    ):
        """
        Prepares the atom data to set the lines, levels and if requested macro
        atom data.  This function mainly cuts the `levels` and `lines` by
        discarding any data that is not needed (any data for atoms that are not
        needed

        Parameters
        ----------
        selected_atoms : set
            set of selected atom numbers, e.g. set([14, 26])
        line_interaction_type : str
            can be 'scatter', 'downbranch' or 'macroatom'
        """
        if not self.prepared:
            self.prepared = True
        else:
            raise AtomDataNotPreparedError("AtomData was already prepared")
        self.selected_atomic_numbers = selected_atomic_numbers

        self._check_selected_atomic_numbers()

        self.nlte_species = nlte_species

        self.levels_index = pd.Series(
            np.arange(len(self.levels), dtype=int), index=self.levels.index
        )

        # cutting levels_lines
        self.lines = self.lines[
            self.lines.index.isin(
                self.selected_atomic_numbers, level="atomic_number"
            )
        ]

        self.lines = self.lines.sort_values(by="wavelength")

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

        self.prepare_macro_atom_data(
            line_interaction_type,
            tmp_lines_lower2level_idx,
            tmp_lines_upper2level_idx,
        )
        if len(continuum_interaction_species) > 0:
            self.prepare_continuum_interaction_data(
                continuum_interaction_species
            )

        self.nlte_data = NLTEData(self, nlte_species)

    def prepare_continuum_interaction_data(self, continuum_interaction_species):
        """
        Prepares the atom data for the continuum interaction

        Parameters
        ----------
        continuum_interaction : ContinuumInteraction
            The continuum interaction object
        """
        # photoionization_data = atomic_data.photoionization_data.set_index(
        #    ["atomic_number", "ion_number", "level_number"]
        # )
        mask_selected_species = self.photoionization_data.index.droplevel(
            "level_number"
        ).isin(continuum_interaction_species)
        self.photoionization_data = self.photoionization_data[
            mask_selected_species
        ]
        self.photo_ion_block_references = np.pad(
            self.photoionization_data.nu.groupby(level=[0, 1, 2])
            .count()
            .values.cumsum(),
            [1, 0],
        )
        self.photo_ion_unique_index = self.photoionization_data.index.unique()
        self.nu_ion_threshold = (
            self.photoionization_data.groupby(level=[0, 1, 2]).first().nu
        )
        self.energy_photo_ion_lower_level = self.levels.loc[
            self.photo_ion_unique_index
        ].energy

        source_idx = self.macro_atom_references.loc[
            self.photo_ion_unique_index
        ].references_idx
        destination_idx = self.macro_atom_references.loc[
            get_ground_state_multi_index(self.photo_ion_unique_index)
        ].references_idx
        self.photo_ion_levels_idx = pd.DataFrame(
            {
                "source_level_idx": source_idx.values,
                "destination_level_idx": destination_idx.values,
            },
            index=self.photo_ion_unique_index,
        )

        self.level2continuum_edge_idx = pd.Series(
            np.arange(len(self.nu_ion_threshold)),
            self.nu_ion_threshold.sort_values(ascending=False).index,
            name="continuum_idx",
        )

        level_idxs2continuum_idx = self.photo_ion_levels_idx.copy()
        level_idxs2continuum_idx[
            "continuum_idx"
        ] = self.level2continuum_edge_idx
        self.level_idxs2continuum_idx = level_idxs2continuum_idx.set_index(
            ["source_level_idx", "destination_level_idx"]
        )

    def prepare_macro_atom_data(
        self,
        line_interaction_type,
        tmp_lines_lower2level_idx,
        tmp_lines_upper2level_idx,
    ):
        if (
            self.macro_atom_data_class.transition_probability_data is not None
            and not line_interaction_type == "scatter"
        ):
            new_data = self.macro_atom_data_class.active_data(
                self.selected_atomic_numbers,
                line_interaction_type,
                self.lines_index,
                tmp_lines_lower2level_idx,
                tmp_lines_upper2level_idx,
            )

            self.macro_atom_data = new_data.transition_probability_data
            self.macro_atom_references = new_data.block_reference_data
            self.lines_lower2macro_reference_idx = (
                new_data.lines_lower2macro_reference_idx
            )
            self.lines_upper2macro_reference_idx = (
                new_data.lines_upper2macro_reference_idx
            )

            # TODO: where does this go?!
            if self.yg_data is not None:
                self.yg_data = self.yg_data.reindex(
                    self.selected_atomic_numbers, level=0
                )

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

            msg = f"For atomic numbers {missing_numbers_str} there is no atomic data."
            raise AtomDataMissingError(msg)

    def __repr__(self):
        return f"<Atomic Data UUID={self.uuid1} MD5={self.md5} Lines={self.lines.line_id.count():d} Levels={self.levels.energy.count():d}>"
