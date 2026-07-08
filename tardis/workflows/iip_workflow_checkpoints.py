import logging
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any

import h5py
import numpy as np
import pandas as pd
from astropy import units as u

from tardis.iip_plasma.standard_plasmas import LegacyPlasmaArray
from tardis.plasma.radiation_field import DilutePlanckianRadiationField

logger = logging.getLogger(__name__)

STRING_DTYPE = h5py.string_dtype(encoding="utf-8")

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


def save_checkpoint(
    workflow: Any,
    normalized_continuum_estimators: dict[str, object] | None = None,
    estimated_values: dict[str, object] | None = None,
) -> None:
    """
    Write the current solved IIP checkpoint to HDF5.

    Parameters
    ----------
    workflow : Any
        IIP workflow with a ``plasma_solver``, ``completed_iterations``, and
        optional base ``checkpoint_path``.
    normalized_continuum_estimators : dict[str, object], optional
        Continuum estimators needed to rebuild the continuum state on resume.
    estimated_values : dict[str, object], optional
        Convergence estimates from the Monte Carlo iteration.
    """
    checkpoint_base_path = getattr(workflow, "checkpoint_path", None)
    if checkpoint_base_path is None:
        return

    checkpoint_path = iteration_checkpoint_path(
        checkpoint_base_path,
        workflow.completed_iterations,
    )
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(checkpoint_path, "w") as hdf_file:
        hdf_file.attrs["serialization"] = "hdf5"

        write_workflow_state(
            hdf_file.create_group("workflow_state"),
            workflow,
            normalized_continuum_estimators,
            estimated_values,
        )
        write_plasma_checkpoint(
            hdf_file.create_group("plasma"),
            workflow.plasma_solver,
        )

    logger.info("Saved IIP checkpoint to %s", checkpoint_path)


def iteration_checkpoint_path(
    checkpoint_path: str | Path, completed_iterations: int
) -> Path:
    """Return the HDF5 checkpoint path for a completed iteration count."""
    checkpoint_base_path = Path(checkpoint_path)
    return checkpoint_base_path.with_name(
        f"{checkpoint_base_path.name}_{completed_iterations}.h5"
    )


def base_checkpoint_path(checkpoint_path: str | Path) -> Path:
    """Return the extensionless checkpoint series base path."""
    checkpoint_path = Path(checkpoint_path)
    base_name = re.sub(r"_\d+$", "", checkpoint_path.stem)
    return checkpoint_path.with_name(base_name)


def write_workflow_state(
    group: h5py.Group,
    workflow: Any,
    normalized_continuum_estimators: dict[str, object] | None,
    estimated_values: dict[str, object] | None,
) -> None:
    """Write workflow-level resume state into an HDF5 group."""
    write_mapping(
        group,
        {
            "completed_iterations": workflow.completed_iterations,
            "consecutive_converges_count": workflow.consecutive_converges_count,
            "converged": workflow.converged,
            "simulation_state_t_radiative": (
                workflow.simulation_state.t_radiative
            ),
            "simulation_state_dilution_factor": (
                workflow.simulation_state.dilution_factor
            ),
            "normalized_continuum_estimators": normalized_continuum_estimators,
            "estimated_values": estimated_values,
        },
    )


def write_plasma_checkpoint(
    group: h5py.Group, plasma: LegacyPlasmaArray
) -> None:
    """Write plasma property state and output values to HDF5."""
    write_mapping(group.create_group("state"), plasma_state(plasma))
    write_property_sequence(
        group.create_group("properties"),
        get_unique_plasma_properties(plasma),
    )
    write_property_sequence(
        group.create_group("previous_iteration_properties"),
        plasma.previous_iteration_properties,
    )
    write_mapping(group.create_group("outputs"), plasma_outputs(plasma))


def plasma_state(plasma: LegacyPlasmaArray) -> dict[str, object]:
    """Return scalar plasma attributes that are independent of properties."""
    skipped_keys = {
        "graph",
        "input_properties",
        "outputs_dict",
        "plasma_properties",
        "previous_iteration_properties",
        "_plasma_properties_dict",
        "_topological_sort_indices",
        "_topological_sort_order",
        "_update_list_cache",
    }
    return filter_checkpoint_mapping(plasma.__dict__, skipped_keys)


def plasma_outputs(plasma: LegacyPlasmaArray) -> dict[str, object]:
    """Return dynamic plasma outputs needed to resume from a checkpoint."""
    return {
        output: value
        for output, plasma_property in plasma.outputs_dict.items()
        if output not in INITIALIZATION_DERIVED_PLASMA_OUTPUTS
        and is_checkpoint_value_supported(
            value := getattr(plasma_property, output)
        )
    }


def write_property_sequence(
    group: h5py.Group, plasma_properties: list[object]
) -> None:
    """Write ordered plasma property objects to an HDF5 group."""
    group.attrs["count"] = len(plasma_properties)
    for index, plasma_property in enumerate(plasma_properties):
        property_group = group.create_group(f"property_{index:04d}")
        property_group.attrs["name"] = plasma_property.name
        property_group.attrs["class_path"] = (
            f"{plasma_property.__class__.__module__}."
            f"{plasma_property.__class__.__qualname__}"
        )
        write_value(
            property_group.create_group("outputs"),
            list(plasma_property.outputs),
        )
        write_mapping(
            property_group.create_group("state"),
            plasma_property_state(plasma_property),
        )


def plasma_property_state(plasma_property: object) -> dict[str, object]:
    """Return checkpointable non-output state for a plasma property."""
    if set(plasma_property.outputs).issubset(
        INITIALIZATION_DERIVED_PLASMA_OUTPUTS
    ):
        return {}

    return filter_checkpoint_mapping(
        plasma_property.__dict__,
        skipped_keys={
            "plasma_parent",
            *plasma_property.outputs,
        },
    )


def filter_checkpoint_mapping(
    mapping: dict[str, object],
    skipped_keys: set[str] | None = None,
) -> dict[str, object]:
    """Return checkpointable values from a mapping."""
    if skipped_keys is None:
        skipped_keys = set()

    return {
        key: value
        for key, value in mapping.items()
        if key not in skipped_keys
        and not callable(value)
        and is_checkpoint_value_supported(value)
    }


def is_checkpoint_value_supported(value: object) -> bool:
    """Return whether a value can be stored by the HDF5 checkpoint writer."""
    if value is None or isinstance(value, str | bool | int | float | complex):
        is_supported = True
    elif isinstance(value, np.generic):
        is_supported = True
    elif isinstance(value, u.Quantity):
        is_supported = is_checkpoint_value_supported(value.value)
    elif isinstance(value, pd.DataFrame):
        is_supported = (
            is_checkpoint_value_supported(value.index)
            and is_checkpoint_value_supported(value.columns)
            and is_checkpoint_value_supported(value.to_numpy())
        )
    elif isinstance(value, pd.Series):
        is_supported = (
            is_checkpoint_value_supported(value.index)
            and is_checkpoint_value_supported(value.name)
            and is_checkpoint_value_supported(value.to_numpy())
        )
    elif isinstance(value, pd.MultiIndex):
        is_supported = is_checkpoint_value_supported(
            value.names
        ) and is_checkpoint_value_supported(
            [tuple(item) for item in value.tolist()]
        )
    elif isinstance(value, pd.Index):
        is_supported = is_checkpoint_value_supported(
            value.name
        ) and is_checkpoint_value_supported(value.tolist())
    elif isinstance(value, np.ndarray):
        if value.dtype.kind in "biufc":
            is_supported = True
        else:
            is_supported = is_checkpoint_value_supported(value.tolist())
    elif isinstance(value, tuple | list):
        is_supported = all(
            is_checkpoint_value_supported(item) for item in value
        )
    elif isinstance(value, dict):
        is_supported = all(
            is_checkpoint_value_supported(key)
            and is_checkpoint_value_supported(dict_value)
            for key, dict_value in value.items()
        )
    else:
        is_supported = False
    return is_supported


def write_mapping(group: h5py.Group, mapping: dict[str, object]) -> None:
    """Write a dictionary as ordered key-value items."""
    group.attrs["kind"] = "dict"
    group.attrs["count"] = len(mapping)
    for index, (key, value) in enumerate(mapping.items()):
        item_group = group.create_group(f"item_{index:04d}")
        write_value(item_group.create_group("key"), key)
        write_value(item_group.create_group("value"), value)


def write_value(group: h5py.Group, value: object) -> None:
    """Write a checkpoint-supported Python object into an HDF5 group."""
    if value is None:
        group.attrs["kind"] = "none"
    elif isinstance(value, str):
        group.attrs["kind"] = "string"
        group.attrs["value"] = value
    elif isinstance(value, bool):
        group.attrs["kind"] = "bool"
        group.attrs["value"] = value
    elif isinstance(value, int) and not isinstance(value, bool):
        group.attrs["kind"] = "int"
        group.attrs["value"] = value
    elif isinstance(value, float):
        group.attrs["kind"] = "float"
        group.attrs["value"] = value
    elif isinstance(value, complex):
        group.attrs["kind"] = "complex"
        group.create_dataset("value", data=value)
    elif isinstance(value, np.generic):
        write_value(group, value.item())
    elif isinstance(value, u.Quantity):
        group.attrs["kind"] = "quantity"
        group.attrs["unit"] = str(value.unit)
        write_value(group.create_group("value"), value.value)
    elif isinstance(value, pd.DataFrame):
        write_dataframe(group, value)
    elif isinstance(value, pd.Series):
        write_series(group, value)
    elif isinstance(value, pd.MultiIndex):
        write_multi_index(group, value)
    elif isinstance(value, pd.Index):
        write_index(group, value)
    elif isinstance(value, np.ndarray):
        write_ndarray(group, value)
    elif isinstance(value, tuple):
        write_sequence(group, value, "tuple")
    elif isinstance(value, list):
        write_sequence(group, value, "list")
    elif isinstance(value, dict):
        write_mapping(group, value)
    else:
        raise TypeError(f"Unsupported checkpoint value type: {type(value)!r}")


def write_dataframe(group: h5py.Group, dataframe: pd.DataFrame) -> None:
    """Write a pandas DataFrame into an HDF5 group."""
    group.attrs["kind"] = "dataframe"
    write_value(group.create_group("index"), dataframe.index)
    write_value(group.create_group("columns"), dataframe.columns)
    write_value(group.create_group("data"), dataframe.to_numpy())
    group.create_dataset(
        "dtypes",
        data=np.asarray(
            [str(dtype) for dtype in dataframe.dtypes], dtype=STRING_DTYPE
        ),
    )


def write_series(group: h5py.Group, series: pd.Series) -> None:
    """Write a pandas Series into an HDF5 group."""
    group.attrs["kind"] = "series"
    group.attrs["dtype"] = str(series.dtype)
    write_value(group.create_group("index"), series.index)
    write_value(group.create_group("name"), series.name)
    write_value(group.create_group("data"), series.to_numpy())


def write_index(group: h5py.Group, index: pd.Index) -> None:
    """Write a pandas Index into an HDF5 group."""
    group.attrs["kind"] = "index"
    write_value(group.create_group("name"), index.name)
    write_value(group.create_group("values"), index.tolist())


def write_multi_index(group: h5py.Group, index: pd.MultiIndex) -> None:
    """Write a pandas MultiIndex into an HDF5 group."""
    group.attrs["kind"] = "multi_index"
    write_value(group.create_group("names"), index.names)
    write_value(
        group.create_group("values"),
        [tuple(item) for item in index.tolist()],
    )


def write_ndarray(group: h5py.Group, array: np.ndarray) -> None:
    """Write a numpy array into an HDF5 group."""
    if array.dtype.kind in "biufc":
        group.attrs["kind"] = "ndarray"
        group.attrs["dtype"] = str(array.dtype)
        group.create_dataset("data", data=array)
    else:
        group.attrs["kind"] = "ndarray_object"
        group.attrs["dtype"] = str(array.dtype)
        group.attrs["shape"] = array.shape
        write_value(group.create_group("data"), array.tolist())


def write_sequence(
    group: h5py.Group, sequence: tuple[object, ...] | list[object], kind: str
) -> None:
    """Write an ordered Python sequence into an HDF5 group."""
    group.attrs["kind"] = kind
    group.attrs["count"] = len(sequence)
    for index, item in enumerate(sequence):
        write_value(group.create_group(f"item_{index:04d}"), item)


def load_checkpoint(checkpoint_path: str | Path) -> dict[str, object]:
    """Load an IIP workflow checkpoint from HDF5."""
    with h5py.File(checkpoint_path, "r") as hdf_file:
        return {
            "serialization": read_text(hdf_file.attrs["serialization"]),
            "workflow_state": read_mapping(hdf_file["workflow_state"]),
            "plasma": read_plasma_checkpoint(hdf_file["plasma"]),
        }


def read_plasma_checkpoint(group: h5py.Group) -> dict[str, object]:
    """Read plasma checkpoint state from an HDF5 group."""
    return {
        "state": read_mapping(group["state"]),
        "properties": read_property_sequence(group["properties"]),
        "previous_iteration_properties": read_property_sequence(
            group["previous_iteration_properties"]
        ),
        "outputs": read_mapping(group["outputs"]),
    }


def read_property_sequence(group: h5py.Group) -> list[dict[str, object]]:
    """Read ordered plasma property state records from HDF5."""
    return [
        {
            "name": read_text(group[f"property_{index:04d}"].attrs["name"]),
            "class_path": read_text(
                group[f"property_{index:04d}"].attrs["class_path"]
            ),
            "outputs": read_value(group[f"property_{index:04d}"]["outputs"]),
            "state": read_mapping(group[f"property_{index:04d}"]["state"]),
        }
        for index in range(int(group.attrs["count"]))
    ]


def read_mapping(group: h5py.Group) -> dict[object, object]:
    """Read a dictionary written by ``write_mapping``."""
    return {
        read_value(group[f"item_{index:04d}"]["key"]): read_value(
            group[f"item_{index:04d}"]["value"]
        )
        for index in range(int(group.attrs["count"]))
    }


def read_value(group: h5py.Group) -> object:
    """Read a checkpoint-supported object from an HDF5 group."""
    value_kind = read_text(group.attrs["kind"])
    if value_kind == "none":
        value = None
    elif value_kind == "string":
        value = read_text(group.attrs["value"])
    elif value_kind == "bool":
        value = bool(group.attrs["value"])
    elif value_kind == "int":
        value = int(group.attrs["value"])
    elif value_kind == "float":
        value = float(group.attrs["value"])
    elif value_kind == "complex":
        value = complex(group["value"][()])
    elif value_kind == "quantity":
        value = read_value(group["value"]) * u.Unit(
            read_text(group.attrs["unit"])
        )
    elif value_kind == "dataframe":
        value = read_dataframe(group)
    elif value_kind == "series":
        value = read_series(group)
    elif value_kind == "index":
        value = pd.Index(
            read_value(group["values"]),
            name=read_value(group["name"]),
        )
    elif value_kind == "multi_index":
        value = pd.MultiIndex.from_tuples(
            [tuple(item) for item in read_value(group["values"])],
            names=read_value(group["names"]),
        )
    elif value_kind == "ndarray":
        value = np.asarray(group["data"], dtype=read_text(group.attrs["dtype"]))
    elif value_kind == "ndarray_object":
        value = np.asarray(read_value(group["data"]), dtype=object).reshape(
            tuple(group.attrs["shape"])
        )
    elif value_kind == "tuple":
        value = tuple(read_sequence(group))
    elif value_kind == "list":
        value = read_sequence(group)
    elif value_kind == "dict":
        value = read_mapping(group)
    else:
        raise ValueError(f"Unknown checkpoint value type: {value_kind}")
    return value


def read_dataframe(group: h5py.Group) -> pd.DataFrame:
    """Read a pandas DataFrame from an HDF5 group."""
    dataframe = pd.DataFrame(
        read_value(group["data"]),
        index=read_value(group["index"]),
        columns=read_value(group["columns"]),
    )
    for column, dtype in zip(
        dataframe.columns, read_string_array(group["dtypes"])
    ):
        try:
            dataframe[column] = dataframe[column].astype(dtype)
        except (TypeError, ValueError) as error:
            del error
            continue
    return dataframe


def read_series(group: h5py.Group) -> pd.Series:
    """Read a pandas Series from an HDF5 group."""
    series = pd.Series(
        read_value(group["data"]),
        index=read_value(group["index"]),
        name=read_value(group["name"]),
    )
    try:
        return series.astype(read_text(group.attrs["dtype"]))
    except (TypeError, ValueError) as error:
        del error
        return series


def read_sequence(group: h5py.Group) -> list[object]:
    """Read a sequence written by ``write_sequence``."""
    return [
        read_value(group[f"item_{index:04d}"])
        for index in range(int(group.attrs["count"]))
    ]


def read_string_array(dataset: h5py.Dataset) -> list[str]:
    """Read a UTF-8 HDF5 string dataset."""
    return [read_text(value) for value in dataset[()]]


def read_text(value: object) -> str:
    """Return a Python string from h5py bytes or string scalars."""
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def get_unique_plasma_properties(plasma: LegacyPlasmaArray) -> list[object]:
    """Return plasma property objects without multi-output duplicates."""
    unique_properties = []
    seen_properties = set()
    for plasma_property in plasma.plasma_properties:
        property_id = id(plasma_property)
        if property_id in seen_properties:
            continue
        seen_properties.add(property_id)
        unique_properties.append(plasma_property)
    return unique_properties


def restore_plasma_checkpoint(
    plasma: LegacyPlasmaArray, checkpoint: dict[str, object]
) -> None:
    """Restore a compatible IIP plasma from an HDF5 checkpoint payload."""
    plasma_checkpoint = checkpoint["plasma"]
    plasma._update_list_cache = {}

    for key, value in plasma_checkpoint["state"].items():
        setattr(plasma, key, value)

    for plasma_property, property_state in zip(
        get_unique_plasma_properties(plasma),
        plasma_checkpoint["properties"],
    ):
        restore_property_state(plasma_property, property_state["state"])
        if hasattr(plasma_property, "plasma_parent"):
            plasma_property.plasma_parent = plasma

    for previous_property, property_state in zip(
        plasma.previous_iteration_properties,
        plasma_checkpoint["previous_iteration_properties"],
    ):
        restore_property_state(previous_property, property_state["state"])

    for output, value in plasma_checkpoint["outputs"].items():
        setattr(plasma.outputs_dict[output], output, value)


def restore_property_state(
    plasma_property: object, state: dict[str, object]
) -> None:
    """Restore checkpointed attributes on a plasma property object."""
    for key, value in state.items():
        setattr(plasma_property, key, value)


def resume_from_checkpoint(
    workflow: Any, checkpoint_path: str | Path
) -> dict[str, object]:
    """Resume an IIP workflow from an HDF5 checkpoint file."""
    checkpoint = load_checkpoint(checkpoint_path)
    restore_workflow_checkpoint(workflow, checkpoint)
    workflow._resumed_from_checkpoint = True
    return checkpoint


def restore_workflow_checkpoint(
    workflow: Any, checkpoint: dict[str, object]
) -> None:
    """Restore plasma and workflow state from a checkpoint payload."""
    restore_plasma_checkpoint(workflow.plasma_solver, checkpoint)

    workflow_state = checkpoint["workflow_state"]
    workflow.completed_iterations = workflow_state["completed_iterations"]
    workflow.consecutive_converges_count = workflow_state[
        "consecutive_converges_count"
    ]
    workflow.converged = workflow_state["converged"]

    t_radiative = workflow_state["simulation_state_t_radiative"]
    dilution_factor = workflow_state["simulation_state_dilution_factor"]
    if t_radiative is not None:
        workflow.simulation_state.t_radiative = t_radiative
    if dilution_factor is not None:
        workflow.simulation_state.dilution_factor = dilution_factor
    if (
        workflow.simulation_state.t_radiative is not None
        and workflow.simulation_state.dilution_factor is not None
    ):
        workflow.simulation_state.radiation_field_state = (
            DilutePlanckianRadiationField(
                workflow.simulation_state.t_radiative,
                workflow.simulation_state.dilution_factor,
            )
        )

    workflow.normalized_continuum_estimators = workflow_state[
        "normalized_continuum_estimators"
    ]
    if workflow.normalized_continuum_estimators is not None:
        workflow.solve_continuum_state(workflow.normalized_continuum_estimators)

    estimated_values = workflow_state["estimated_values"]
    if estimated_values is not None:
        workflow.converged = workflow.check_convergence(estimated_values)
        workflow.completed_iterations += 1
