from tardis.io.model.readers.snec.xg_files import read_xg_file, XGData
import numpy as np
import yaml
import pandas as pd
import numpy.testing as npt
from dataclasses import dataclass
import xarray as xr

with open("parser_config/snec_xg_output_quantities.yml", "r") as fh:
    SNEC_XG_OUTPUT_METADATA = yaml.safe_load(fh)

with open("parser_config/snec_initial_composition.yml") as fh:
    SNEC_INITIAL_COMPOSITION_METADATA = yaml.safe_load(fh)

with open("parser_config/snec_initial_quantities.yml") as fh:
    SNEC_INITIAL_QUANTITIES_METADATA = yaml.safe_load(fh)

with open("parser_config/snec_em_output_metadata.yml") as fh:
    SNEC_EM_OUTPUT_METADATA = yaml.safe_load(fh)


@dataclass
class SNECOutput:
    xg_data: XGData
    initial_composition: pd.DataFrame
    initial_quantities: pd.DataFrame
    em_output: pd.DataFrame

    def to_xr_dataset(self) -> xr.Dataset:
        """Convert SNECOutput to xarray.Dataset combining xg_data, em_output, initial_composition, and initial_quantities."""
        # Base dataset from xg_data
        ds = self.xg_data.to_xr_dataset()

        # Merge emission output over time, interpolating to xg_data time using nearest
        em_ds = self.em_output.set_index("time").to_xarray()
        em_ds = em_ds.reindex(time=ds.time, method="nearest")
        ds = ds.merge(em_ds)

        # Merge initial composition as coordinates/data_vars over cell_id
        comp_ds = self.initial_composition.set_index("cell_id").to_xarray()
        ds = ds.merge(comp_ds)

        # Merge initial quantities as coordinates/data_vars over cell_id
        quant_ds = self.initial_quantities.set_index("cell_id").to_xarray()
        ds = ds.merge(quant_ds)

        return ds


def read_snec_output_xg(snec_output_dir, show_progress=False):
    all_xg_data = {}
    mass_xg_data = read_xg_file(
        snec_output_dir / "output" / "mass.xg",
        column_names=["radius", "enclosed_mass"],
    )

    for output_quantity, metadata in SNEC_XG_OUTPUT_METADATA.items():
        xg_file = snec_output_dir / "output" / f"{output_quantity}.xg"
        if not xg_file.exists():
            raise FileNotFoundError(f"File {xg_file} does not exist.")

        xg_data = read_xg_file(
            xg_file,
            column_names=["enclosed_mass", output_quantity],
            show_progress=show_progress,
        )

        # Ensure timestamps match
        assert np.all(
            np.isclose(mass_xg_data.timestamps, xg_data.timestamps)
        ), f"Time stamps do not match for {output_quantity} and mass."

        # Add metadata to XGData
        xg_data.metadata = metadata
        all_xg_data[output_quantity] = xg_data

    merged_data_blocks = []
    for i, time_stamp in enumerate(mass_xg_data.timestamps):
        merged_df = mass_xg_data.data_blocks[i][["radius"]].copy()
        merged_df["enclosed_mass"] = mass_xg_data.data_blocks[i]["enclosed_mass"]

        for quantity_name, xg_data in all_xg_data.items():
            merged_df[quantity_name] = xg_data.data_blocks[i][quantity_name]

        merged_data_blocks.append(merged_df)

    # Create a single XGData object with all merged data
    merged_xg_data = XGData(
        timestamps=mass_xg_data.timestamps,
        data_blocks=merged_data_blocks,
        metadata=SNEC_XG_OUTPUT_METADATA,  # Use the full SNEC_XG_OUTPUT_QUANTITIES as metadata
    )
    return merged_xg_data


def read_snec_dat_output(snec_output_dir, snec_dat_names, first_column_name):
    """
    Load the initial composition data from SNEC output.

    Parameters
    ----------
    snec_output_dir : Path
        Path to the directory containing the SNEC output files.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the initial composition for each cell, indexed by cell ID.
    """
    first_composition_name = snec_dat_names[0]
    snec_initial_composition_df = pd.read_csv(
        snec_output_dir / "output" / f"{first_composition_name}.dat",
        names=[first_column_name, first_composition_name],
        sep=r"\s+",
    )

    for snec_initial_composition_name in snec_dat_names[1:]:
        current_snec_initial_composition_df = pd.read_csv(
            snec_output_dir / "output" / f"{snec_initial_composition_name}.dat",
            names=[first_column_name, snec_initial_composition_name],
            sep=r"\s+",
        )
        npt.assert_array_almost_equal(
            snec_initial_composition_df.iloc[:, 0],
            current_snec_initial_composition_df.iloc[:, 0],
        )

        snec_initial_composition_df[snec_initial_composition_name] = (
            current_snec_initial_composition_df[snec_initial_composition_name]
        )

    return snec_initial_composition_df


def read_snec_initial_composition(snec_output_dir):
    """
    Load the initial composition data from SNEC output.

    Parameters
    ----------
    snec_output_dir : Path
        Path to the directory containing the SNEC output files.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the initial composition for each cell, indexed by cell ID.
    """
    return read_snec_dat_output(
        snec_output_dir,
        list(SNEC_INITIAL_COMPOSITION_METADATA.keys()),
        first_column_name="cell_id",
    )


def read_snec_initial_quantities(snec_output_dir):
    """
    Load the initial quantities data from SNEC output.

    Parameters
    ----------
    snec_output_dir : Path
        Path to the directory containing the SNEC output files.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the initial quantities for each cell, indexed by cell ID.
    """
    return read_snec_dat_output(
        snec_output_dir,
        list(SNEC_INITIAL_QUANTITIES_METADATA.keys()),
        first_column_name="cell_id",
    )


def read_snec_em_output(snec_output_dir):
    """
    Load the output quantities data from SNEC output.

    Parameters
    ----------
    snec_output_dir : Path
        Path to the directory containing the SNEC output files.

    Returns
    -------
    pd.DataFrame
        DataFrame containing the output quantities for each cell, indexed by cell ID.
    """

    em_output_metadata = [
        item for item in SNEC_EM_OUTPUT_METADATA.keys() if not item.startswith("index_")
    ]
    em_index_output_metadata = [
        item for item in SNEC_EM_OUTPUT_METADATA.keys() if item.startswith("index_")
    ]

    em_output_df = read_snec_dat_output(
        snec_output_dir,
        em_output_metadata,
        first_column_name="time",
    )

    em_index_output_df = read_snec_dat_output(
        snec_output_dir,
        em_index_output_metadata,
        first_column_name="time",
    )

    npt.assert_allclose(em_output_df["time"], em_index_output_df["time"], rtol=1e-9)

    return em_output_df.join(em_index_output_df.iloc[:, 1:])


def read_snec_output(snec_output_dir, show_progress=False):
    """
    Read SNEC output files and return a SNECOutput dataclass instance.

    Parameters
    ----------
    snec_output_dir : Path
        Path to the directory containing the SNEC output files.
    show_progress : bool, optional
        If True, show progress for reading xg files. Default is False.

    Returns
    -------
    SNECOutput
        Dataclass instance containing the SNEC output data.
    """
    return SNECOutput(
        xg_data=read_snec_output_xg(snec_output_dir, show_progress),
        initial_composition=read_snec_initial_composition(snec_output_dir),
        initial_quantities=read_snec_initial_quantities(snec_output_dir),
        em_output=read_snec_em_output(snec_output_dir),
    )
