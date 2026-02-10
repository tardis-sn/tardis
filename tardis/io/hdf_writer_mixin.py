from __future__ import annotations

import re
from importlib.metadata import version as get_version
from pathlib import Path
from typing import Any, Union, Optional, Dict, List

import numpy as np
import pandas as pd

# Use importlib.metadata to get version - works even with namespace collision (Issue #3021)
__version__ = get_version("tardis")
from tardis.io.util import logger


class HDFWriterMixin:
    """
    Mixin class that provides HDF writing functionality to other classes.
    """

    def __new__(cls, *args: Any, **kwargs: Any):
        """
        Create a new instance and initialize HDF properties.

        Parameters
        ----------
        *args : Any
            Positional arguments passed to the parent class.
        **kwargs : Any
            Keyword arguments passed to the parent class.

        Returns
        -------
        HDFWriterMixin
            The new instance with initialized HDF properties.
        """
        instance = super().__new__(cls)
        instance.optional_hdf_properties = []
        instance.__init__(*args, **kwargs)
        return instance

    @staticmethod
    def to_hdf_util(
        path_or_buf: Union[str, pd.HDFStore],
        path: str,
        elements: Dict[str, Any],
        overwrite: bool,
        complevel: int = 9,
        complib: str = "blosc",
        format=None,
    ) -> None:
        """
        A function to uniformly store TARDIS data to an HDF file.

        Scalars will be stored in a Series under path/scalars.
        1D arrays will be stored under path/property_name as distinct Series.
        2D arrays will be stored under path/property_name as distinct DataFrames.

        Units will be stored as their CGS value.

        Parameters
        ----------
        path_or_buf : str or pandas.HDFStore
            Path or buffer to the HDF file.
        path : str
            Path inside the HDF file to store the `elements`.
        elements : dict
            A dict of property names and their values to be stored.
        overwrite : bool
            If the HDF file path already exists, whether to overwrite it or not.
        complevel : int, optional
            Compression level for the HDF file, by default 9.
        complib : str, optional
            Compression library to use, by default "blosc".

        Notes
        -----
        `overwrite` option doesn't have any effect when `path_or_buf` is an
        HDFStore because the user decides on the mode in which they have
        opened the HDFStore ('r', 'w' or 'a').

        Raises
        ------
        FileExistsError
            If the HDF file already exists and overwrite is False.
        """
        if (
            isinstance(path_or_buf, str)
            and Path(path_or_buf).exists()
            and not overwrite
        ):
            raise FileExistsError(
                "The specified HDF file already exists. If you still want "
                "to overwrite it, set function parameter overwrite=True"
            )

        try:  # when path_or_buf is a str, the HDFStore should get created
            buf = pd.HDFStore(
                path_or_buf, complevel=complevel, complib=complib  # type: ignore[arg-type]
            )
        except TypeError as e:
            if str(e) == "Expected bytes, got HDFStore":
                # when path_or_buf is an HDFStore buffer instead
                logger.debug(
                    "Expected bytes, got HDFStore. Changing path to HDF buffer"
                )
                buf = path_or_buf
            else:
                raise e

        if not buf.is_open:
            buf.open()

        scalars = {}
        for key, value in elements.items():
            if value is None:
                value = "none"
            if hasattr(value, "cgs"):
                value = value.cgs.value
            if np.isscalar(value):
                scalars[key] = value
            elif hasattr(value, "shape"):
                if value.ndim == 1:
                    # This try,except block is only for model.plasma.levels
                    try:
                        pd.Series(value).to_hdf(
                            buf, key=str(Path(path) / key), format=format
                        )
                    except NotImplementedError:
                        logger.debug(
                            "Could not convert SERIES to HDF. Converting DATAFRAME to HDF"
                        )
                        pd.DataFrame(value).to_hdf(
                            buf, key=str(Path(path) / key), format=format
                        )
                else:
                    pd.DataFrame(value).to_hdf(
                        buf, key=str(Path(path) / key), format=format
                    )
            else:  # value is a TARDIS object like model, transport or plasma
                try:
                    value.to_hdf(
                        buf, path, name=key, overwrite=overwrite, format=format
                    )
                except AttributeError:
                    logger.debug(
                        "Could not convert VALUE to HDF. Converting DATA (Dataframe) to HDF"
                    )
                    data = pd.DataFrame([value])
                    data.to_hdf(buf, key=str(Path(path) / key), format=format)

        metadata = {"tardis_version": __version__}
        pd.Series(metadata).to_hdf(
            buf,
            key=str(Path(path) / "metadata"),
        )

        if scalars:
            pd.Series(scalars).to_hdf(buf, key=str(Path(path) / "scalars"))

        if buf.is_open:
            buf.close()

    def get_properties(self) -> Dict[str, Any]:
        """
        Get properties for HDF storage.

        Returns
        -------
        Dict[str, Any]
            Dictionary of property names and their values.
        """
        data = {name: getattr(self, name) for name in self.full_hdf_properties}
        return data

    @property
    def full_hdf_properties(self) -> List[str]:
        """
        Get the full list of HDF properties.

        Returns
        -------
        List[str]
            List of all HDF property names.
        """
        if hasattr(self, "virt_logging") and self.virt_logging:
            self.hdf_properties.extend(self.vpacket_hdf_properties)

        return self.optional_hdf_properties + self.hdf_properties

    @staticmethod
    def convert_to_snake_case(s: str) -> str:
        """
        Convert a camelCase or PascalCase string to snake_case.

        Parameters
        ----------
        s : str
            The string to convert.

        Returns
        -------
        str
            The string converted to snake_case.
        """
        s1 = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", s)
        return re.sub("([a-z0-9])([A-Z])", r"\1_\2", s1).lower()

    def to_hdf(
        self,
        file_path_or_buf: Union[str, pd.HDFStore],
        path: str = "",
        name: Optional[str] = None,
        overwrite: bool = False,
        format: Optional[str] = None,
    ) -> None:
        """
        Save the object to an HDF file.

        Parameters
        ----------
        file_path_or_buf : str or pandas.HDFStore
            Path or buffer to the HDF file.
        path : str, optional
            Path inside the HDF file to store the `elements`, by default "".
        name : str, optional
            Group inside the HDF file to which the `elements` need to be saved.
            If None, will use the class name converted to snake_case.
        overwrite : bool, optional
            If the HDF file path already exists, whether to overwrite it or not,
            by default False.
        """
        if name is None:
            try:
                name = self.hdf_name
            except AttributeError:
                name = self.convert_to_snake_case(self.__class__.__name__)
                logger.debug(
                    "self.hdf_name not present, setting name to %s for HDF", name
                )

        data = self.get_properties()
        buff_path = str(Path(path) / (name or ""))
        self.to_hdf_util(
            file_path_or_buf, buff_path, data, overwrite=overwrite, format=format
        )


class PlasmaWriterMixin(HDFWriterMixin):
    """
    Mixin class for writing plasma data to HDF files.
    """

    def get_properties(self) -> Dict[str, Any]:
        """
        Get plasma properties for HDF storage.

        Returns
        -------
        Dict[str, Any]
            Dictionary of plasma property names and their values.
        """
        data = {}
        if self.collection:
            properties = [
                name
                for name in self.plasma_properties
                if isinstance(name, tuple(self.collection))
            ]
        else:
            properties = self.plasma_properties
        for prop in properties:
            for output in prop.outputs:
                data[output] = getattr(prop, output)
        data["atom_data_uuid"] = self.atomic_data.uuid1
        if "atomic_data" in data:
            data.pop("atomic_data")
        if "nlte_data" in data:
            logger.warning("nlte_data can't be saved")
            data.pop("nlte_data")
        return data

    def to_hdf(
        self,
        file_path_or_buf: Union[str, pd.HDFStore],
        path: str = "",
        name: Optional[str] = None,
        collection: Any = None,
        overwrite: bool = False,
        format: Optional[str] = None,
    ) -> None:
        """
        Save the plasma object to an HDF file.

        Parameters
        ----------
        file_path_or_buf : str or pandas.HDFStore
            Path or buffer to the HDF file.
        path : str, optional
            Path inside the HDF file to store the `elements`, by default "".
        name : str, optional
            Group inside the HDF file to which the `elements` need to be saved.
        collection : Any, optional
            `None` or a `PlasmaPropertyCollection` of which members are
            the property types which will be stored. If `None` then
            all types of properties will be stored. This acts like a filter,
            for example if a value of `property_collections.basic_inputs` is
            given, only those input parameters will be stored to the HDF file.
        overwrite : bool, optional
            If the HDF file path already exists, whether to overwrite it or not,
            by default False.
        """
        self.collection = collection
        super().to_hdf(file_path_or_buf, path, name, overwrite, format=format)
