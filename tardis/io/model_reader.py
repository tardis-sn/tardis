# reading different model files

from tardis.io.config_reader import ConfigurationNameSpace
from tardis.montecarlo.base import MontecarloRunner
from tardis.util.base import parse_quantity, is_valid_nuclide_or_elem

import warnings
import numpy as np
from numpy import recfromtxt, genfromtxt
import pandas as pd
from astropy import units as u
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z
import h5py

import logging

# Adding logging support
logger = logging.getLogger(__name__)


class ConfigurationError(Exception):
    pass


def read_density_file(filename, filetype):
    """
    read different density file formats

    Parameters
    ----------
    filename : str
        filename or path of the density file
    filetype : str
        type of the density file

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    velocity : np.ndarray
        the array containing the velocities
    unscaled_mean_densities : np.ndarray
        the array containing the densities
    """
    file_parsers = {
        "artis": read_artis_density,
        "simple_ascii": read_simple_ascii_density,
        "cmfgen_model": read_cmfgen_density,
    }

    electron_densities = None
    temperature = None
    if filetype == "cmfgen_model":
        (
            time_of_model,
            velocity,
            unscaled_mean_densities,
            electron_densities,
            temperature,
        ) = read_cmfgen_density(filename)
    else:
        (time_of_model, velocity, unscaled_mean_densities) = file_parsers[
            filetype
        ](filename)

    v_inner = velocity[:-1]
    v_outer = velocity[1:]

    invalid_volume_mask = (v_outer - v_inner) <= 0
    if invalid_volume_mask.sum() > 0:
        message = "\n".join(
            [
                f"cell {i:d}: v_inner {v_inner_i:s}, v_outer " f"{v_outer_i:s}"
                for i, v_inner_i, v_outer_i in zip(
                    np.arange(len(v_outer))[invalid_volume_mask],
                    v_inner[invalid_volume_mask],
                    v_outer[invalid_volume_mask],
                )
            ]
        )
        raise ConfigurationError(
            "Invalid volume of following cell(s):\n" f"{message:s}"
        )

    return (
        time_of_model,
        velocity,
        unscaled_mean_densities,
        electron_densities,
        temperature,
    )


def read_abundances_file(
    abundance_filename,
    abundance_filetype,
    inner_boundary_index=None,
    outer_boundary_index=None,
):
    """
    read different density file formats

    Parameters
    ----------
    abundance_filename : str
        filename or path of the density file
    abundance_filetype : str
        type of the density file
    inner_boundary_index : int
        index of the inner shell, default None
    outer_boundary_index : int
        index of the outer shell, default None
    """

    file_parsers = {
        "simple_ascii": read_simple_ascii_abundances,
        "artis": read_simple_ascii_abundances,
        "cmfgen_model": read_cmfgen_composition,
        "custom_composition": read_csv_composition,
    }

    isotope_abundance = pd.DataFrame()
    if abundance_filetype in ["cmfgen_model", "custom_composition"]:
        index, abundances, isotope_abundance = file_parsers[abundance_filetype](
            abundance_filename
        )
    else:
        index, abundances = file_parsers[abundance_filetype](abundance_filename)

    if outer_boundary_index is not None:
        outer_boundary_index_m1 = outer_boundary_index - 1
    else:
        outer_boundary_index_m1 = None
    index = index[inner_boundary_index:outer_boundary_index]
    abundances = abundances.loc[
        :, slice(inner_boundary_index, outer_boundary_index_m1)
    ]
    abundances.columns = np.arange(len(abundances.columns))
    return index, abundances, isotope_abundance


def read_uniform_abundances(abundances_section, no_of_shells):
    """
    Parameters
    ----------
    abundances_section : config.model.abundances
    no_of_shells : int

    Returns
    -------
    abundance : pandas.DataFrame
    isotope_abundance : pandas.DataFrame
    """
    abundance = pd.DataFrame(
        columns=np.arange(no_of_shells),
        index=pd.Index(np.arange(1, 120), name="atomic_number"),
        dtype=np.float64,
    )

    isotope_index = pd.MultiIndex(
        [[]] * 2, [[]] * 2, names=["atomic_number", "mass_number"]
    )
    isotope_abundance = pd.DataFrame(
        columns=np.arange(no_of_shells), index=isotope_index, dtype=np.float64
    )

    for element_symbol_string in abundances_section:
        if element_symbol_string == "type":
            continue
        try:
            if element_symbol_string in Z_DICT.values():
                z = elem_to_Z(element_symbol_string)
                abundance.loc[z] = float(
                    abundances_section[element_symbol_string]
                )
            else:
                nuc = Nuclide(element_symbol_string)
                mass_no = nuc.A
                z = nuc.Z
                isotope_abundance.loc[(z, mass_no), :] = float(
                    abundances_section[element_symbol_string]
                )

        except RuntimeError as err:
            raise RuntimeError(
                f"Abundances are not defined properly in config file : {err.args}"
            )

    return abundance, isotope_abundance


def read_simple_ascii_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    5 s
    #index velocity [km/s] density [g/cm^3]
    0 1.1e4 1.6e8
    1 1.2e4 1.7e8

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    data : pandas.DataFrame
        data frame containing index, velocity (in km/s) and density
    """

    with open(fname) as fh:
        time_of_model_string = fh.readline().strip()
        time_of_model = parse_quantity(time_of_model_string)

    data = recfromtxt(
        fname,
        skip_header=1,
        names=("index", "velocity", "density"),
        dtype=(int, float, float),
    )
    velocity = (data["velocity"] * u.km / u.s).to("cm/s")
    mean_density = (data["density"] * u.Unit("g/cm^3"))[1:]

    return time_of_model, velocity, mean_density


def read_artis_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    5
    #index velocity [km/s] log10(density) [log10(g/cm^3)]
    0 1.1e4 1.6e8
    1 1.2e4 1.7e8

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    data : pandas.DataFrame
        data frame containing index, velocity (in km/s) and density
    """

    with open(fname) as fh:
        for i, line in enumerate(open(fname)):
            if i == 0:
                no_of_shells = np.int64(line.strip())
            elif i == 1:
                time_of_model = u.Quantity(float(line.strip()), "day").to("s")
            elif i == 2:
                break

    artis_model_columns = [
        "index",
        "velocities",
        "mean_densities_0",
        "ni56_fraction",
        "co56_fraction",
        "fe52_fraction",
        "cr48_fraction",
    ]
    artis_model = recfromtxt(
        fname,
        skip_header=2,
        usecols=(0, 1, 2, 4, 5, 6, 7),
        unpack=True,
        dtype=[(item, np.float64) for item in artis_model_columns],
    )

    velocity = u.Quantity(artis_model["velocities"], "km/s").to("cm/s")
    mean_density = u.Quantity(10 ** artis_model["mean_densities_0"], "g/cm^3")[
        1:
    ]

    return time_of_model, velocity, mean_density


def read_cmfgen_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    The first density describes the mean density in the center of the model and is not used.
    The file consists of a header row and next row contains unit of the respective attributes
    Note that the first column has to contain a running index

    Example:

    index velocity densities electron_densities temperature
    - km/s g/cm^3 /cm^3 K
    0 871.66905 4.2537191e-09 2.5953807e+14 7.6395577
    1 877.44269 4.2537191e-09 2.5953807e+14 7.6395577

    Rest columns contain abundances of elements and isotopes

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    time_of_model : astropy.units.Quantity
        time at which the model is valid
    velocity : np.ndarray
    mean_density : np.ndarray
    electron_densities : np.ndarray
    temperature : np.ndarray
    """
    warnings.warn(
        "The current CMFGEN model parser is deprecated", DeprecationWarning
    )

    df = pd.read_csv(fname, comment="#", delimiter=r"\s+", skiprows=[0, 2])

    with open(fname) as fh:
        for row_index, line in enumerate(fh):
            if row_index == 0:
                time_of_model_string = line.strip().replace("t0:", "")
                time_of_model = parse_quantity(time_of_model_string)
            elif row_index == 2:
                quantities = line.split()

    velocity = u.Quantity(df["velocity"].values, quantities[1]).to("cm/s")
    temperature = u.Quantity(df["temperature"].values, quantities[2])[1:]
    mean_density = u.Quantity(df["densities"].values, quantities[3])[1:]
    electron_densities = u.Quantity(
        df["electron_densities"].values, quantities[4]
    )[1:]

    return (
        time_of_model,
        velocity,
        mean_density,
        electron_densities,
        temperature,
    )


def read_simple_ascii_abundances(fname):
    """
    Reading an abundance file of the following structure (example; lines starting with hash will be ignored):
    The first line of abundances describe the abundances in the center of the model and are not used.
    #index element1, element2, ..., element30
    0 0.4 0.3, .. 0.2

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    index : np.ndarray
        containing the indices
    abundances : pandas.DataFrame
        data frame containing index, element1 - element30 and columns according to the shells
    """
    data = np.loadtxt(fname)

    index = data[1:, 0].astype(int)
    abundances = pd.DataFrame(
        data[1:, 1:].transpose(), index=np.arange(1, data.shape[1])
    )

    return index, abundances


def read_cmfgen_composition(fname, delimiter=r"\s+"):
    """Read composition from a CMFGEN model file

    The CMFGEN file format contains information about the ejecta state in the
    first four columns and the following ones contain elemental and isotopic
    abundances.

    WARNING : deprecated

    fname : str
        filename of the csv file
    """

    warnings.warn(
        "The current CMFGEN model parser is deprecated", DeprecationWarning
    )

    return read_csv_isotope_abundances(
        fname, delimiter=delimiter, skip_columns=4, skip_rows=[0, 2, 3]
    )


def read_csv_composition(fname, delimiter=r"\s+"):
    """Read composition from a simple CSV file

    The CSV file can contain specific isotopes or elemental abundances in the
    different columns. The first row must contain the header in which the
    contents of each column is specified by the elemental symbol (for elemental
    abundances) or by the symbol plus mass number (for isotopic abundances).

    Example: C O Fe Ni56 Co

    The i-th row specifies the composition in the i-th shell

    fname : str
        filename of the csv file
    """

    return read_csv_isotope_abundances(
        fname, delimiter=delimiter, skip_columns=0, skip_rows=[1]
    )


def read_csv_isotope_abundances(
    fname, delimiter=r"\s+", skip_columns=0, skip_rows=[1]
):
    """
    A generic parser for a TARDIS composition stored as a CSV file

    The parser can read in both elemental and isotopic abundances. The first
    column is always expected to contain a running index, labelling the grid
    cells. The parser also allows for additional information to be stored in
    the first skip_columns columns. These will be ignored if skip_columns > 0.
    Note that the first column, containing the cell index is not taken into
    account here.

    Specific header lines can be skipped by the skip_rows keyword argument

    It is expected that the first row of the date block (after skipping the
    rows specified in skip_rows) specifies the different elements and isotopes.
    Each row after contains the composition in the corresponding grid shell.
    The first composition row describes the composition of the photosphere and
    is essentially ignored (for the default value of skip_rows).

    Example:

    Index C   O   Ni56
    0     1   1   1
    1     0.4 0.3 0.2

    Parameters
    ----------
    fname : str
        filename or path with filename

    Returns
    -------
    index : np.ndarray
    abundances : pandas.DataFrame
    isotope_abundance : pandas.MultiIndex
    """

    df = pd.read_csv(
        fname, comment="#", sep=delimiter, skiprows=skip_rows, index_col=0
    )
    df = df.transpose()

    abundance = pd.DataFrame(
        columns=np.arange(df.shape[1]),
        index=pd.Index([], name="atomic_number"),
        dtype=np.float64,
    )

    isotope_index = pd.MultiIndex(
        [[]] * 2, [[]] * 2, names=["atomic_number", "mass_number"]
    )
    isotope_abundance = pd.DataFrame(
        columns=np.arange(df.shape[1]), index=isotope_index, dtype=np.float64
    )

    for element_symbol_string in df.index[skip_columns:]:
        if element_symbol_string in Z_DICT.values():
            z = elem_to_Z(element_symbol_string)
            abundance.loc[z, :] = df.loc[element_symbol_string].tolist()
        else:
            nuc = Nuclide(element_symbol_string)
            z = nuc.Z
            mass_no = nuc.A
            isotope_abundance.loc[(z, mass_no), :] = df.loc[
                element_symbol_string
            ].tolist()

    return abundance.index, abundance, isotope_abundance


def parse_csv_abundances(csvy_data):
    """
    A parser for the csv data part of a csvy model file. This function filters out columns that are not abundances.

    Parameters
    ----------
    csvy_data : pandas.DataFrame

    Returns
    -------
    index : np.ndarray
    abundances : pandas.DataFrame
    isotope_abundance : pandas.MultiIndex
    """

    abundance_col_names = [
        name for name in csvy_data.columns if is_valid_nuclide_or_elem(name)
    ]
    df = csvy_data.loc[:, abundance_col_names]

    df = df.transpose()

    abundance = pd.DataFrame(
        columns=np.arange(df.shape[1]),
        index=pd.Index([], name="atomic_number"),
        dtype=np.float64,
    )

    isotope_index = pd.MultiIndex(
        [[]] * 2, [[]] * 2, names=["atomic_number", "mass_number"]
    )
    isotope_abundance = pd.DataFrame(
        columns=np.arange(df.shape[1]), index=isotope_index, dtype=np.float64
    )

    for element_symbol_string in df.index[0:]:
        if element_symbol_string in Z_DICT.values():
            z = elem_to_Z(element_symbol_string)
            abundance.loc[z, :] = df.loc[element_symbol_string].tolist()
        else:
            nuc = Nuclide(element_symbol_string)
            z = nuc.Z
            mass_no = nuc.A
            isotope_abundance.loc[(z, mass_no), :] = df.loc[
                element_symbol_string
            ].tolist()

    return abundance.index, abundance, isotope_abundance


def runner_to_dict(runner):
    """
    Retrieves all the data from a runner object and returns a dictionary.

    Parameters
    ----------
    runner : tardis.montecarlo.MontecarloRunner

    Returns
    -------
    runner_dict : dict
    integrator_settings : dict
    v_packet_settings : dict
    virtual_spectrum_spawn_range : dict
    """

    runner_dict = {
        "Edotlu_estimator": runner.Edotlu_estimator,
        "bf_heating_estimator": runner.bf_heating_estimator,
        "disable_electron_scattering": runner.disable_electron_scattering,
        "enable_full_relativity": runner.enable_full_relativity,
        "enable_reflective_inner_boundary": runner.enable_reflective_inner_boundary,
        "inner_boundary_albedo": runner.inner_boundary_albedo,
        "input_energy": runner.input_energy,
        "input_mu": runner.input_mu,
        "input_nu": runner.input_nu,
        "input_r_cgs": runner.input_r,
        "j_blue_estimator": runner.j_blue_estimator,
        "j_estimator": runner.j_estimator,
        "last_interaction_in_nu": runner.last_interaction_in_nu,
        "last_interaction_type": runner.last_interaction_type,
        "last_line_interaction_in_id": runner.last_line_interaction_in_id,
        "last_line_interaction_out_id": runner.last_line_interaction_out_id,
        "last_line_interaction_shell_id": runner.last_line_interaction_shell_id,
        "line_interaction_type": runner.line_interaction_type,
        "nu_bar_estimator": runner.nu_bar_estimator,
        "photo_ion_estimator": runner.photo_ion_estimator,
        "photo_ion_estimator_statistics": runner.photo_ion_estimator_statistics,
        "r_inner": runner.r_inner_cgs,
        "r_outer": runner.r_outer_cgs,
        "seed": runner.seed,
        "spectrum_frequency_cgs": runner.spectrum_frequency,
        "spectrum_method": runner.spectrum_method,
        "stim_recomb_cooling_estimator": runner.stim_recomb_cooling_estimator,
        "stim_recomb_estimator": runner.stim_recomb_estimator,
        "time_of_simulation_cgs": runner.time_of_simulation,
        "use_gpu": runner.use_gpu,
        "v_inner": runner.v_inner_cgs,
        "v_outer": runner.v_outer_cgs,
        "virt_logging": runner.virt_logging,
        "virt_packet_energies": runner.virt_packet_energies,
        "virt_packet_initial_mus": runner.virt_packet_initial_mus,
        "virt_packet_initial_rs": runner.virt_packet_initial_rs,
        "virt_packet_last_interaction_in_nu": runner.virt_packet_last_interaction_in_nu,
        "virt_packet_last_interaction_type": runner.virt_packet_last_interaction_type,
        "virt_packet_last_line_interaction_in_id": runner.virt_packet_last_line_interaction_in_id,
        "virt_packet_last_line_interaction_out_id": runner.virt_packet_last_line_interaction_out_id,
        "virt_packet_nus": runner.virt_packet_nus,
        "volume_cgs": runner.volume,
    }

    try:
        runner_dict["single_packet_seed"] = runner.single_packet_seed
    except AttributeError:
        runner_dict["single_packet_seed"] = None

    for key, value in runner_dict.items():
        if key.endswith("_cgs"):
            runner_dict[key] = [value.cgs.value, value.unit.to_string()]

    integrator_settings = runner.integrator_settings
    v_packet_settings = runner.v_packet_settings
    virtual_spectrum_spawn_range = runner.virtual_spectrum_spawn_range

    return (
        runner_dict,
        integrator_settings,
        v_packet_settings,
        virtual_spectrum_spawn_range,
    )


def store_runner_to_hdf(runner, fname):
    """
    Stores data from MontecarloRunner object into a hdf file.

    Parameters
    ----------
    runner : tardis.montecarlo.MontecarloRunner
    filename : str
    """

    with h5py.File(fname, "a") as f:
        runner_group = f.require_group("runner")
        runner_group.clear()

        (
            runner_data,
            integrator_settings,
            v_packet_settings,
            virtual_spectrum_spawn_range,
        ) = runner_to_dict(runner)

        for key, value in runner_data.items():
            if key.endswith("_cgs"):
                runner_group.create_dataset(key, data=value[0])
                runner_group.create_dataset(key + "_unit", data=value[1])
            else:
                if value is not None:
                    runner_group.create_dataset(key, data=value)

        integrator_settings_group = runner_group.create_group(
            "integrator_settings"
        )
        for key, value in integrator_settings.items():
            integrator_settings_group.create_dataset(key, data=value)

        v_packet_settings_group = runner_group.create_group("v_packet_settings")
        for key, value in v_packet_settings.items():
            v_packet_settings_group.create_dataset(key, data=value)

        virtual_spectrum_spawn_range_group = runner_group.create_group(
            "virtual_spectrum_spawn_range"
        )
        for key, value in virtual_spectrum_spawn_range.items():
            virtual_spectrum_spawn_range_group.create_dataset(key, data=value)


def runner_from_hdf(fname):
    """
    Creates a MontecarloRunner object using data stored in a hdf file.

    Parameters
    ----------
    fname : str

    Returns
    -------
    new_runner : tardis.montecarlo.MontecarloRunner
    """

    d = {"single_packet_seed": None}

    # Loading data from hdf file
    with h5py.File(fname, "r") as f:
        runner_group = f["runner"]
        for key, value in runner_group.items():
            if not key.endswith("_unit"):
                if isinstance(value, h5py.Dataset):
                    if isinstance(value[()], bytes):
                        d[key] = value[()].decode("utf-8")
                    else:
                        d[key] = value[()]
                else:
                    data_inner = {}
                    for key_inner, value_inner in value.items():
                        if isinstance(value_inner[()], bytes):
                            data_inner[key] = value_inner[()].decode("utf-8")
                        else:
                            data_inner[key] = value_inner[()]
                        data_inner[key_inner] = value_inner[()]
                    d[key] = data_inner

        for key, value in runner_group.items():
            if key.endswith("_unit"):
                d[key[:-5]] = [d[key[:-5]], value[()]]

    # Converting cgs data to astropy quantities
    for key, value in d.items():
        if key.endswith("_cgs"):
            d[key] = u.Quantity(value[0], unit=u.Unit(value[1].decode("utf-8")))

    # Converting virtual spectrum spawn range values to astropy quantities
    vssr = d["virtual_spectrum_spawn_range"]
    d["virtual_spectrum_spawn_range"] = {
        "start": u.Quantity(vssr["start"], unit=u.Unit("Angstrom")),
        "end": u.Quantity(vssr["end"], unit=u.Unit("Angstrom")),
    }

    # Converting dictionaries to ConfigurationNameSpace
    d["integrator_settings"] = ConfigurationNameSpace(d["integrator_settings"])
    d["v_packet_settings"] = ConfigurationNameSpace(d["v_packet_settings"])
    d["virtual_spectrum_spawn_range"] = ConfigurationNameSpace(
        d["virtual_spectrum_spawn_range"]
    )

    # Creating a runner object and storing data
    new_runner = MontecarloRunner(
        seed=d["seed"],
        spectrum_frequency=d["spectrum_frequency_cgs"],
        virtual_spectrum_spawn_range=d["virtual_spectrum_spawn_range"],
        disable_electron_scattering=d["disable_electron_scattering"],
        enable_reflective_inner_boundary=d["enable_reflective_inner_boundary"],
        enable_full_relativity=d["enable_full_relativity"],
        inner_boundary_albedo=d["inner_boundary_albedo"],
        line_interaction_type=d["line_interaction_type"],
        integrator_settings=d["integrator_settings"],
        v_packet_settings=d["v_packet_settings"],
        spectrum_method=d["spectrum_method"],
        virtual_packet_logging=d["virt_logging"],
        single_packet_seed=d["single_packet_seed"],
        use_gpu=d["use_gpu"],
    )

    new_runner.Edotlu_estimator = d["Edotlu_estimator"]
    new_runner.bf_heating_estimator = d["bf_heating_estimator"]
    new_runner.input_energy = d["input_energy"]
    new_runner.input_mu = d["input_mu"]
    new_runner.input_nu = d["input_nu"]
    new_runner.input_r = d["input_r_cgs"]
    new_runner.j_blue_estimator = d["j_blue_estimator"]
    new_runner.j_estimator = d["j_estimator"]
    new_runner.last_interaction_in_nu = d["last_interaction_in_nu"]
    new_runner.last_interaction_type = d["last_interaction_type"]
    new_runner.last_line_interaction_in_id = d["last_line_interaction_in_id"]
    new_runner.last_line_interaction_out_id = d["last_line_interaction_out_id"]
    new_runner.last_line_interaction_shell_id = d[
        "last_line_interaction_shell_id"
    ]
    new_runner.nu_bar_estimator = d["nu_bar_estimator"]
    new_runner.photo_ion_estimator = d["photo_ion_estimator"]
    new_runner.photo_ion_estimator_statistics = d[
        "photo_ion_estimator_statistics"
    ]
    new_runner.r_inner_cgs = d["r_inner"]
    new_runner.r_outer_cgs = d["r_outer"]
    new_runner.stim_recomb_cooling_estimator = d[
        "stim_recomb_cooling_estimator"
    ]
    new_runner.stim_recomb_estimator = d["stim_recomb_estimator"]
    new_runner.time_of_simulation = d["time_of_simulation_cgs"]
    new_runner.v_inner_cgs = d["v_inner"]
    new_runner.v_outer_cgs = d["v_outer"]
    new_runner.virt_packet_energies = d["virt_packet_energies"]
    new_runner.virt_packet_initial_mus = d["virt_packet_initial_mus"]
    new_runner.virt_packet_initial_rs = d["virt_packet_initial_rs"]
    new_runner.virt_packet_last_interaction_in_nu = d[
        "virt_packet_last_interaction_in_nu"
    ]
    new_runner.virt_packet_last_interaction_type = d[
        "virt_packet_last_interaction_type"
    ]
    new_runner.virt_packet_last_line_interaction_in_id = d[
        "virt_packet_last_line_interaction_in_id"
    ]
    new_runner.virt_packet_last_line_interaction_out_id = d[
        "virt_packet_last_line_interaction_out_id"
    ]
    new_runner.virt_packet_nus = d["virt_packet_nus"]
    new_runner.volume = d["volume_cgs"]

    return new_runner


def model_to_dict(model):
    """
    Retrieves all the data from a Radial1DModel object and returns a dictionary.

    Parameters
    ----------
    runner : tardis.model.Radial1DModel

    Returns
    -------
    model_dict : dict
    homologous_density : dict
    isotope_abundance : dict
    """
    model_dict = {
        "velocity_cgs": model.velocity,
        "abundance": model.abundance,
        "time_explosion_cgs": model.time_explosion,
        "t_inner_cgs": model.t_inner,
        "t_radiative_cgs": model.t_radiative,
        "dilution_factor": model.dilution_factor,
        "v_boundary_inner_cgs": model.v_boundary_inner,
        "v_boundary_outer_cgs": model.v_boundary_outer,
        "w": model.w,
        "t_rad_cgs": model.t_rad,
        "r_inner_cgs": model.r_inner,
        "density_cgs": model.density,
    }

    for key, value in model_dict.items():
        if hasattr(value, "unit"):
            model_dict[key] = [value.cgs.value, value.unit.to_string()]

    homologous_density = model.homologous_density.__dict__
    isotope_abundance = model.raw_isotope_abundance.__dict__

    return model_dict, homologous_density, isotope_abundance


def store_model_to_hdf(model, fname):
    """
    Stores data from Radial1DModel object into a hdf file.

    Parameters
    ----------
    model : tardis.model.Radial1DModel
    filename : str
    """
    with h5py.File(fname, "a") as f:
        model_group = f.require_group("model")
        model_group.clear()

        model_dict, homologous_density, isotope_abundance = model_to_dict(model)

        for key, value in model_dict.items():
            if key.endswith("_cgs"):
                model_group.create_dataset(key, data=value[0])
                model_group.create_dataset(key + "_unit", data=value[1])
            else:
                model_group.create_dataset(key, data=value)

        homologous_density_group = model_group.create_group(
            "homologous_density"
        )
        for key, value in homologous_density.items():
            homologous_density_group.create_dataset(key, data=value)


def model_from_hdf(fname):
    """
    Creates a Radial1DModel object using data stored in a hdf file.

    Parameters
    ----------
    fname : str

    Returns
    -------
    new_model : tardis.model.Radial1DModel
    """

    from tardis.model import Radial1DModel, HomologousDensity

    d = {}

    # Loading data from hdf file
    with h5py.File(fname, "r") as f:
        model_group = f["model"]
        for key, value in model_group.items():
            if not key.endswith("_unit"):
                if type(value) == h5py._hl.dataset.Dataset:
                    if isinstance(value[()], bytes):
                        d[key] = value[()].decode("utf-8")
                    else:
                        d[key] = value[()]
                else:
                    data_inner = {}
                    for key_inner, value_inner in value.items():
                        if isinstance(value_inner[()], bytes):
                            data_inner[key] = value_inner[()].decode("utf-8")
                        else:
                            data_inner[key] = value_inner[()]
                        data_inner[key_inner] = value_inner[()]
                    d[key] = data_inner

        for key, value in model_group.items():
            if key.endswith("_unit"):
                d[key[:-5]] = [d[key[:-5]], value[()]]

    # Converting cgs data to astropy quantities
    for key, value in d.items():
        if key.endswith("_cgs"):
            d[key] = u.Quantity(value[0], unit=u.Unit(value[1].decode("utf-8")))

    homologous_density = HomologousDensity(
        d["homologous_density"]["density_0"], d["homologous_density"]["time_0"]
    )

    new_model = Radial1DModel(
        velocity=d["velocity_cgs"],
        homologous_density=homologous_density,
        abundance=d["abundance"],
        isotope_abundance=None,
        time_explosion=d["time_explosion_cgs"],
        t_inner=d["t_inner_cgs"],
        t_radiative=d["t_radiative_cgs"],
        dilution_factor=d["dilution_factor"],
        v_boundary_inner=d["v_boundary_inner_cgs"],
        v_boundary_outer=d["v_boundary_outer_cgs"],
    )

    new_model.t_rad = d["t_rad_cgs"]
    new_model.w = d["w"]

    return new_model
