from pathlib import Path

import pytest
from syrupy.extensions.amber import AmberSnapshotExtension
from syrupy.location import PyTestLocation

SNAPSHOT_LOCATION = "plasma"


@pytest.fixture
def snapshot_pd(snapshot, tardis_snapshot_path, pandas_snapshot_extention):
    refpath = tardis_snapshot_path.joinpath(SNAPSHOT_LOCATION)

    class PandasSnapshotExtenstionRefdata(pandas_snapshot_extention):
        @classmethod
        def dirname(cls, *, test_location: "PyTestLocation") -> str:
            return str(Path(test_location.filepath).parent.joinpath(refpath))

    return snapshot.use_extension(PandasSnapshotExtenstionRefdata)


@pytest.fixture
def snapshot_np(snapshot, tardis_snapshot_path, numpy_snapshot_extension):
    refpath = tardis_snapshot_path.joinpath(SNAPSHOT_LOCATION)

    class NumpySnapshotExtenstionRefdata(numpy_snapshot_extension):
        @classmethod
        def dirname(cls, *, test_location: "PyTestLocation") -> str:
            return str(Path(test_location.filepath).parent.joinpath(refpath))

    return snapshot.use_extension(NumpySnapshotExtenstionRefdata)


@pytest.fixture
def snapshot(snapshot, tardis_snapshot_path):
    refpath = tardis_snapshot_path.joinpath(SNAPSHOT_LOCATION)

    # TODO: try this with SingleFileSnapshotExtension
    # test files from single test function would be grouped together
    class AmberSnapshotExtensionRefdata(AmberSnapshotExtension):
        @classmethod
        def dirname(cls, *, test_location: "PyTestLocation") -> str:
            return str(Path(test_location.filepath).parent.joinpath(refpath))

    return snapshot.use_extension(AmberSnapshotExtensionRefdata)
