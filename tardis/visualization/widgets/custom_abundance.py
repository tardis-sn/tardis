"""Class to create and display Custom Abundance Widget."""
import os
import logging
import numpy as np
import pandas as pd
import ipywidgets as ipw

from tardis.io.config_reader import Configuration
from tardis.io.config_validator import validate_dict
from tardis.io.parsers.csvy import load_csvy
from tardis.io.model_reader import (
    read_density_file,
    read_uniform_abundances,
    parse_csv_abundances,
)

logger = logging.getLogger(__name__)


class CustomAbundanceWidget:
    def __init__(self, config, abundance):
        self.config = config  # Save the config for output.

        self.abundance = abundance
        self.no_of_shell = abundance.shape[1]

    @classmethod
    def from_csvy(cls, fname):
        csvy_model_config, csvy_model_data = load_csvy(fname)
        base_dir = os.path.abspath(os.path.dirname(__file__))
        schema_dir = os.path.join(base_dir, "../..", "io", "schemas")
        csvy_schema_file = os.path.join(schema_dir, "csvy_model.yml")
        csvy_model_config = Configuration(
            validate_dict(csvy_model_config, schemapath=csvy_schema_file)
        )

        if hasattr(csvy_model_config, "velocity"):
            no_of_shells = csvy_model_config.velocity.num

        if hasattr(csvy_model_config, "abundance"):
            abundances_section = csvy_model_config.abundance
            abundance, isotope_abundance = read_uniform_abundances(
                abundances_section, no_of_shells
            )
        else:
            _, abundance, isotope_abundance = parse_csv_abundances(
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

        # Combine element and isotope to one DataFrame
        abundance["mass_number"] = 0.0
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        norm_factor = abundance.sum(axis=0)

        if np.any(np.abs(norm_factor - 1) > 1e-12):
            logger.warning(
                "Abundances have not been normalized to 1." " - normalizing"
            )
            abundance /= norm_factor

        return cls(csvy_model_config, abundance)

    @classmethod
    def from_yml(cls, fname):
        config = Configuration.from_yaml(fname)

        if hasattr(config, "csvy_model"):
            raise ValueError(
                "The configuration file contains 'csvy_model' module. 'from_yml()' only supports uniform abundance, please call 'from_csvy()' instead."
            )
        else:
            abundances_section = config.model.abundances

            if abundances_section.type == "uniform":
                structure = config.model.structure

                if structure.type == "specific":
                    no_of_shells = structure.velocity.num

                elif structure.type == "file":
                    if os.path.isabs(structure.filename):
                        structure_fname = structure.filename
                    else:
                        structure_fname = os.path.join(
                            config.config_dirname, structure.filename
                        )

                    _, velocity, _, _, _ = read_density_file(
                        structure_fname, structure.filetype
                    )
                    # Note: This is the number of shells *without* taking in mind the
                    #       v boundaries.
                    no_of_shells = len(velocity) - 1

                else:
                    raise NotImplementedError

                abundance, isotope_abundance = read_uniform_abundances(
                    abundances_section, no_of_shells
                )

                abundance = abundance.replace(np.nan, 0.0)
                abundance = abundance[abundance.sum(axis=1) > 0]
                isotope_abundance = isotope_abundance.replace(np.nan, 0.0)
                isotope_abundance = isotope_abundance[
                    isotope_abundance.sum(axis=1) > 0
                ]

                # integrate element and isotope to one DataFrame
                abundance["mass_number"] = 0.0
                abundance.set_index("mass_number", append=True, inplace=True)
                abundance = pd.concat([abundance, isotope_abundance])
                abundance.sort_index(inplace=True)

                norm_factor = abundance.sum(axis=0)

                if np.any(np.abs(norm_factor - 1) > 1e-12):
                    logger.warning(
                        "Abundances have not been normalized to 1."
                        " - normalizing"
                    )
                    abundance /= norm_factor

            else:
                raise ValueError(
                    "'from_yml()' only supports uniform abundance currently."
                )

        return cls(config, abundance)

    @classmethod
    def from_hdf(cls, fname):
        pass

    @classmethod
    def from_simulation(cls, sim):
        abundance = sim.model.raw_abundance.copy()
        isotope_abundance = sim.model.raw_isotope_abundance.copy()

        # integrate element and isotope to one DataFrame
        abundance["mass_number"] = 0.0
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        return cls(config=None, abundance=abundance)
