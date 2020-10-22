import yaml
import pandas as pd
from tardis.io.util import YAMLLoader

YAML_DELIMITER = "---"


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
            raise ValueError("End %s not found" % (YAML_DELIMITER))
        yaml_dict = yaml.load("".join(yaml_lines[1:-1]), YAMLLoader)
        try:
            data = pd.read_csv(fname, skiprows=yaml_end_ind + 1)
        except pd.errors.EmptyDataError as e:
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
            raise ValueError("End %s not found" % (YAML_DELIMITER))
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
