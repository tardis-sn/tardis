"""Class to create and display Custom Abundance Widget."""
import os
import logging
import numpy as np
import pandas as pd
from astropy import units as u

from tardis.util.base import quantity_linspace
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
    def __init__(self, config, abundance, velocity):
        self.config = config  # Save the config for output.

        self.abundance = abundance
        self.velocity = np.array(velocity) # unit: cm/s
        self.no_of_shell = abundance.shape[1]

    @classmethod
    def from_csvy(cls, fpath):
        csvy_model_config, csvy_model_data = load_csvy(fpath)
        base_dir = os.path.abspath(os.path.dirname(__file__))
        schema_dir = os.path.join(base_dir, "../..", "io", "schemas")
        csvy_schema_file = os.path.join(schema_dir, "csvy_model.yml")
        csvy_model_config = Configuration(
            validate_dict(csvy_model_config, schemapath=csvy_schema_file)
        )

        if hasattr(csvy_model_config, "velocity"):
            velocity = quantity_linspace(csvy_model_config.velocity.start,
                                         csvy_model_config.velocity.stop,
                                         csvy_model_config.velocity.num + 1).cgs
        else:
            velocity_field_index = [field['name'] for field in csvy_model_config.datatype.fields].index('velocity')
            velocity_unit = u.Unit(csvy_model_config.datatype.fields[velocity_field_index]['unit'])
            velocity = csvy_model_data['velocity'].values * velocity_unit
            velocity = velocity.to('cm/s')

        no_of_shells = len(velocity) - 1

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

        # Combine elements and isotopes to one DataFrame
        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        norm_factor = abundance.sum(axis=0)

        if np.any(np.abs(norm_factor - 1) > 1e-12):
            logger.warning(
                "Abundances have not been normalized to 1." " - normalizing"
            )
            abundance /= norm_factor

        return cls(config=csvy_model_config, abundance=abundance, velocity=velocity)

    @classmethod
    def from_yml(cls, fpath):
        config = Configuration.from_yaml(fpath)

        if hasattr(config, "csvy_model"):
            raise ValueError(
                "The configuration file contains 'csvy_model' module. 'from_yml()' only supports uniform abundance, please call 'from_csvy()' instead."
            )
        else:
            abundances_section = config.model.abundances

            if abundances_section.type == "uniform":
                structure = config.model.structure

                if structure.type == "specific":
                    velocity = quantity_linspace(structure.velocity.start,
                                         structure.velocity.stop,
                                         structure.velocity.num + 1).cgs

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

                else:
                    raise NotImplementedError
                # Note: This is the number of shells *without* taking in mind the
                #       v boundaries.
                no_of_shells = len(velocity) - 1

                abundance, isotope_abundance = read_uniform_abundances(
                    abundances_section, no_of_shells
                )

                abundance = abundance.replace(np.nan, 0.0)
                abundance = abundance[abundance.sum(axis=1) > 0]
                isotope_abundance = isotope_abundance.replace(np.nan, 0.0)
                isotope_abundance = isotope_abundance[
                    isotope_abundance.sum(axis=1) > 0
                ]

                # Combine elements and isotopes to one DataFrame
                abundance["mass_number"] = ""
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

        return cls(config=config, abundance=abundance, velocity=velocity)

    @classmethod
    def from_hdf(cls, fpath):
        with pd.HDFStore(fpath, "r") as hdf:
            abundance = (hdf['/simulation/plasma/abundance'])
            v_inner = hdf['/simulation/model/v_inner'].to_numpy()
            v_outer = hdf['/simulation/model/v_outer'].to_numpy()
            velocity = np.append(v_inner, v_outer[-1])

        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        
        return cls(config=None, abundance=abundance, velocity=velocity)


    @classmethod
    def from_simulation(cls, sim):
        abundance = sim.model.raw_abundance.copy()
        isotope_abundance = sim.model.raw_isotope_abundance.copy()

        # integrate element and isotope to one DataFrame
        abundance["mass_number"] = ""
        abundance.set_index("mass_number", append=True, inplace=True)
        abundance = pd.concat([abundance, isotope_abundance])
        abundance.sort_index(inplace=True)

        velocity = sim.model.velocity

        return cls(config=None, abundance=abundance, velocity=velocity)
