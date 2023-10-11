import logging
from astropy import units as u
import numpy as np
from radioactivedecay import Nuclide
from radioactivedecay.utils import Z_DICT, elem_to_Z
import yaml
import pandas as pd
from tardis.io.util import YAMLLoader
from tardis.model.geometry.radial1d import HomologousRadial1DGeometry
from tardis.util.base import is_valid_nuclide_or_elem, quantity_linspace

YAML_DELIMITER = "---"

logger = logging.getLogger(__name__)


def load_csvy(fname):
    """
    Parameters
    ----------
    fname : string
            Path to csvy file

    Returns
    -------
    yaml_dict : dictionary
                YAML part of the csvy file
    data : pandas.dataframe
            csv data from csvy file
    """
    with open(fname) as fh:
        yaml_lines = []
        yaml_end_ind = -1
        for i, line in enumerate(fh):
            if i == 0:
                assert (
                    line.strip() == YAML_DELIMITER
                ), "First line of csvy file is not '---'"
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_end_ind = i
                break
        else:
            raise ValueError(f"End {YAML_DELIMITER} not found")
        yaml_dict = yaml.load("".join(yaml_lines[1:-1]), YAMLLoader)
        try:
            data = pd.read_csv(fname, skiprows=yaml_end_ind + 1)
        except pd.errors.EmptyDataError as e:
            logger.debug(f"Could not Read CSV. Setting Dataframe to None")
            data = None

    return yaml_dict, data


def load_yaml_from_csvy(fname):
    """
    Parameters
    ----------
    fname : string
            Path to csvy file

    Returns
    -------
    yaml_dict : dictionary
                YAML part of the csvy file
    """
    with open(fname) as fh:
        yaml_lines = []
        yaml_end_ind = -1
        for i, line in enumerate(fh):
            if i == 0:
                assert (
                    line.strip() == YAML_DELIMITER
                ), "First line of csvy file is not '---'"
            yaml_lines.append(line)
            if i > 0 and line.strip() == YAML_DELIMITER:
                yaml_end_ind = i
                break
        else:
            raise ValueError(f"End {YAML_DELIMITER} not found")
        yaml_dict = yaml.load("".join(yaml_lines[1:-1]), YAMLLoader)
    return yaml_dict


def load_csv_from_csvy(fname):
    """
    Parameters
    ----------
    fname : string
            Path to csvy file

    Returns
    -------
    data : pandas.dataframe
           csv data from csvy file
    """
    yaml_dict, data = load_csvy(fname)
    return data


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


def parse_csvy_geometry(
    config, csvy_model_config, csvy_model_data, time_explosion
):
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

    geometry = HomologousRadial1DGeometry(
        velocity[:-1],  # r_inner
        velocity[1:],  # r_outer
        v_inner_boundary=v_boundary_inner,
        v_outer_boundary=v_boundary_outer,
        time_explosion=time_explosion,
    )
    return geometry
