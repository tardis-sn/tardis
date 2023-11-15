import re
from typing import Any, List, Tuple

import numpy as np
import pandas as pd
from syrupy.data import SnapshotCollection
from syrupy.extensions.single_file import SingleFileSnapshotExtension, WriteMode
from syrupy.location import PyTestLocation
from syrupy.types import SerializableData, SnapshotIndex


class SingleFileSanitizedNames(SingleFileSnapshotExtension):
    # changing write mode to text helps avoid an error message
    # that comes when files are serialised in syrupy in bytes
    # either way we won't be serialising files in most cases in bytes
    _write_mode = WriteMode.TEXT
    _file_extension = "txt"

    # would change names of all snapshots generated
    # that use this class- making filenames compliant with python standards.
    @classmethod
    def get_snapshot_name(
        cls, *, test_location: "PyTestLocation", index: "SnapshotIndex"
    ) -> str:
        original_name = SingleFileSnapshotExtension.get_snapshot_name(
            test_location=test_location, index=index
        )
        double_under = r"[:\[\]{}]"
        no_space = r'[,"\']'  # quotes and commas

        name = re.sub(double_under, "__", original_name)
        name = re.sub(no_space, "", name)

        return f"{name}"


class NumpySnapshotExtenstion(SingleFileSanitizedNames):
    _file_extension = "npy"

    def matches(self, *, serialized_data, snapshot_data):
        try:
            if (
                np.testing.assert_allclose(
                    np.array(snapshot_data), np.array(serialized_data)
                )
                is not None
            ):
                return False
            else:
                return True

        except:
            return False

    def _read_snapshot_data_from_location(
        self, *, snapshot_location: str, snapshot_name: str, session_id: str
    ):
        # see https://github.com/tophat/syrupy/blob/f4bc8453466af2cfa75cdda1d50d67bc8c4396c3/src/syrupy/extensions/base.py#L139
        try:
            return np.load(snapshot_location)

        except OSError:
            return None

    @classmethod
    def _write_snapshot_collection(
        cls, *, snapshot_collection: SnapshotCollection
    ) -> None:
        # see https://github.com/tophat/syrupy/blob/f4bc8453466af2cfa75cdda1d50d67bc8c4396c3/src/syrupy/extensions/base.py#L161

        filepath, data = (
            snapshot_collection.location,
            next(iter(snapshot_collection)).data,
        )

        np.save(filepath, data)

    def serialize(self, data: SerializableData, **kwargs: Any) -> str:
        return data


class PandasSnapshotExtenstion(SingleFileSanitizedNames):
    _file_extension = "h5"

    def matches(self, *, serialized_data, snapshot_data):
        try:
            comparer = {
                pd.Series: pd.testing.assert_series_equal,
                pd.DataFrame: pd.testing.assert_frame_equal,
            }
            try:
                comp_func = comparer[type(serialized_data)]
            except KeyError:
                raise ValueError(
                    "Can only compare Series and Dataframes with PandasSnapshotExtenstion."
                )

            if comp_func(serialized_data, snapshot_data) is not None:
                return False
            else:
                return True

        except:
            return False

    def _read_snapshot_data_from_location(
        self, *, snapshot_location: str, snapshot_name: str, session_id: str
    ):
        # see https://github.com/tophat/syrupy/blob/f4bc8453466af2cfa75cdda1d50d67bc8c4396c3/src/syrupy/extensions/base.py#L139
        try:
            data = pd.read_hdf(snapshot_location)
            return data

        except OSError:
            return None

    @classmethod
    def _write_snapshot_collection(
        cls, *, snapshot_collection: SnapshotCollection
    ) -> None:
        # see https://github.com/tophat/syrupy/blob/f4bc8453466af2cfa75cdda1d50d67bc8c4396c3/src/syrupy/extensions/base.py#L161
        filepath, data = (
            snapshot_collection.location,
            next(iter(snapshot_collection)).data,
        )
        data.to_hdf(filepath, "/data")

    def serialize(self, data: SerializableData, **kwargs: Any) -> str:
        return data
