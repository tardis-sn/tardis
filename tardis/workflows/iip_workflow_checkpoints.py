import json
import logging
import re
from collections.abc import Callable
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import __version__
from tardis.plasma.radiation_field import DilutePlanckianRadiationField
from tardis.workflows.type_iip_workflow import TypeIIPWorkflow

logger = logging.getLogger(__name__)

# Package this into a reader class that handles each type in the schema with a specific reader method
# That way we don't have the huge if else blocks to check every possible thing
SCALAR_READERS: dict[str, Callable[[str], object]] = {
    "none": lambda value: None,
    "bool": lambda value: value == "True",
    "int": int,
    "float": float,
    "complex": complex,
    "str": lambda value: value,
}

INITIALIZATION_DERIVED_PLASMA_OUTPUTS = frozenset(
    {
        "abundance",
        "atomic_data",
        "atomic_mass",
        "continuum_data",
        "density",
        "excitation_energy",
        "f_lu",
        "g",
        "ionization_data",
        "levels",
        "lines",
        "lines_lower_level_index",
        "lines_multi_index",
        "lines_upper_level_index",
        "link_t_rad_t_electron",
        "metastability",
        "nlte_data",
        "nu",
        "number_density",
        "photo_ion_cross_sections",
        "photo_ion_index_sorted",
        "selected_atoms",
        "T_Yg",
        "time_explosion",
        "wavelength_cm",
        "Yg_index",
        "Yg_interp",
        "yg_allowed_index",
        "yg_data",
        "yg_forbidden_index",
        "zeta_data",
    }
)

# There's surely a better way to do this with Path natively, this is just redundant and clumsy
def iteration_checkpoint_path(
    checkpoint_path: str | Path, completed_iterations: int
) -> Path:
    checkpoint_base_path = Path(checkpoint_path)
    return checkpoint_base_path.with_name(
        f"{checkpoint_base_path.name}_{completed_iterations}.h5"
    )


def base_checkpoint_path(checkpoint_path: str | Path) -> Path:
    checkpoint_path = Path(checkpoint_path)
    return checkpoint_path.with_name(re.sub(r"_\d+$", "", checkpoint_path.stem))


# HDF serialization


def _write_hdf_value(
    store: pd.HDFStore, path: str, value: object
) -> dict[str, str]:
    # Rewrite into schema reader writer class
    schema = {
        "path": path,
        "kind": "",
        "value": "",
        "unit": "",
        "index_names": "",
    }

    if value is None:
        schema["kind"] = "none"
        return schema

    if isinstance(value, np.generic):
        value = value.item()

    if isinstance(value, bool | int | float | complex | str):
        schema["kind"] = type(value).__name__
        schema["value"] = value if isinstance(value, str) else repr(value)
        return schema

    if isinstance(value, u.Quantity):
        schema["kind"] = f"quantity_{value.ndim}d"
        schema["unit"] = str(value.unit)
        value = value.value

    if isinstance(value, np.ndarray):
        if not schema["kind"]:
            schema["kind"] = f"array_{value.ndim}d"
        pandas_value = (
            pd.Series(value, copy=False)
            if value.ndim == 1
            else pd.DataFrame(value, copy=False)
        )
    elif isinstance(value, pd.DataFrame):
        schema["kind"] = "dataframe"
        pandas_value = value
    elif isinstance(value, pd.Series):
        schema["kind"] = "series"
        pandas_value = value
    elif isinstance(value, pd.MultiIndex):
        schema["kind"] = "multiindex"
        schema["index_names"] = json.dumps(list(value.names))
        pandas_value = value.to_frame(index=False)
    elif isinstance(value, list):
        schema["kind"] = "array_list"
        schema["value"] = str(len(value))
        if not value:
            return schema
        pandas_value = pd.DataFrame(np.stack(value), copy=False)
    else:
        raise TypeError(f"Unsupported checkpoint value: {type(value)!r}")

    store.put(path, pandas_value, format="fixed")
    return schema


def _read_hdf_value(store: pd.HDFStore, schema: object) -> object:
    # Should be rewritten to a single dict lookup. Collect all the reader functions
    # into a dict or something similar. (use a class)
    kind = schema.kind
    if kind in SCALAR_READERS:
        return SCALAR_READERS[kind](schema.value)

    if kind == "array_list" and schema.value == "0":
        return []

    value = store[schema.path]
    if kind in {"dataframe", "series"}:
        return value
    if kind == "multiindex":
        return pd.MultiIndex.from_frame(
            value,
            names=json.loads(schema.index_names),
        )

    array = value.to_numpy()
    if kind == "array_list":
        result = [row.copy() for row in array]
    elif kind.startswith("quantity_"):
        result = array * u.Unit(schema.unit)
    else:
        result = array
    return result


def _write_checkpoint_file(
    checkpoint_path: Path, checkpoint: dict[str, object]
) -> None:
    workflow_state = checkpoint["workflow_state"]
    plasma_state = checkpoint["plasma"]
    metadata = pd.DataFrame(
        [
            {
                "tardis_version": __version__,
                "completed_iterations": workflow_state["completed_iterations"],
                "consecutive_converges_count": workflow_state[
                    "consecutive_converges_count"
                ],
                "converged": workflow_state["converged"],
                "niter": plasma_state["niter"],
                "niter_ly": plasma_state["niter_ly"],
                "plasma_converged": plasma_state["plasma_converged"],
            }
        ]
    )
    # This is a horrible list comprehension that needs to be rewritten
    # NO **{}
    values = {
        "workflow/t_radiative": workflow_state["t_radiative"],
        "workflow/dilution_factor": workflow_state["dilution_factor"],
        **{
            f"workflow/estimators/{name}": value
            for name, value in workflow_state[
                "normalized_continuum_estimators"
            ].items()
        },
        **{
            f"plasma/outputs/{name}": value
            for name, value in plasma_state["outputs"].items()
        },
    }

    with pd.HDFStore(
        checkpoint_path,
        mode="w",
        complevel=9,
        complib="blosc",
    ) as store:
        schema = [
            _write_hdf_value(store, path, value)
            for path, value in values.items()
        ]
        store.put("metadata", metadata, format="table")
        store.put("schema", pd.DataFrame(schema), format="table")


def save_checkpoint(
    workflow: TypeIIPWorkflow,
    normalized_continuum_estimators: dict[str, object],
) -> None:
    checkpoint_base_path = workflow.configuration.checkpoints.path
    if checkpoint_base_path is None:
        return

    checkpoint_path = iteration_checkpoint_path(
        checkpoint_base_path,
        workflow.completed_iterations,
    )
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    plasma = workflow.plasma_solver
    checkpoint = {
        "workflow_state": {
            "completed_iterations": workflow.completed_iterations,
            "consecutive_converges_count": (
                workflow.consecutive_converges_count
            ),
            "converged": workflow.converged,
            "t_radiative": workflow.simulation_state.t_radiative,
            "dilution_factor": workflow.simulation_state.dilution_factor,
            "normalized_continuum_estimators": (
                normalized_continuum_estimators
            ),
        },
        "plasma": {
            "niter": plasma.niter,
            "niter_ly": plasma.niter_ly,
            "plasma_converged": plasma.plasma_converged,
            "outputs": {
                output: getattr(plasma_property, output)
                for output, plasma_property in plasma.outputs_dict.items()
                if output not in INITIALIZATION_DERIVED_PLASMA_OUTPUTS
            },
        },
    }
    _write_checkpoint_file(checkpoint_path, checkpoint)
    logger.info("Saved IIP checkpoint to %s", checkpoint_path)


def load_checkpoint(checkpoint_path: str | Path) -> dict[str, object]:
    with pd.HDFStore(checkpoint_path, mode="r") as store:
        metadata = store["metadata"].iloc[0]
        values = {
            schema.path: _read_hdf_value(store, schema)
            for schema in store["schema"].itertuples(index=False)
        }

    estimator_prefix = "workflow/estimators/"
    output_prefix = "plasma/outputs/"
    return {
        "workflow_state": {
            "completed_iterations": int(metadata.completed_iterations),
            "consecutive_converges_count": int(
                metadata.consecutive_converges_count
            ),
            "converged": bool(metadata.converged),
            "t_radiative": values["workflow/t_radiative"],
            "dilution_factor": values["workflow/dilution_factor"],
            "normalized_continuum_estimators": {
                path.removeprefix(estimator_prefix): value
                for path, value in values.items()
                if path.startswith(estimator_prefix)
            },
        },
        "plasma": {
            "niter": int(metadata.niter),
            "niter_ly": int(metadata.niter_ly),
            "plasma_converged": bool(metadata.plasma_converged),
            "outputs": {
                path.removeprefix(output_prefix): value
                for path, value in values.items()
                if path.startswith(output_prefix)
            },
        },
    }


def resume_from_checkpoint(
    workflow: TypeIIPWorkflow, checkpoint_path: str | Path
) -> None:
    checkpoint = load_checkpoint(checkpoint_path)
    workflow_state = checkpoint["workflow_state"]
    plasma_state = checkpoint["plasma"]
    plasma = workflow.plasma_solver

    plasma.niter = plasma_state["niter"]
    plasma.niter_ly = plasma_state["niter_ly"]
    plasma.plasma_converged = plasma_state["plasma_converged"]
    for output, value in plasma_state["outputs"].items():
        setattr(plasma.outputs_dict[output], output, value)

    workflow.completed_iterations = workflow_state["completed_iterations"]
    workflow.consecutive_converges_count = workflow_state[
        "consecutive_converges_count"
    ]
    workflow.converged = workflow_state["converged"]
    workflow.simulation_state.radiation_field_state = (
        DilutePlanckianRadiationField(
            workflow_state["t_radiative"],
            workflow_state["dilution_factor"],
        )
    )
    workflow.solve_continuum_state(
        workflow_state["normalized_continuum_estimators"]
    )
