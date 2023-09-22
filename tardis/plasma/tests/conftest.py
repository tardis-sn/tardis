from pathlib import Path
import pytest
from syrupy.location import PyTestLocation
from syrupy.extensions.amber import AmberSnapshotExtension

SNAPSHOT_LOCATION = "snapshot_data/plasma"

@pytest.fixture
def snapshot_pd(snapshot, tardis_ref_path, pandas_snapshot_extention):
    refpath = tardis_ref_path.joinpath(SNAPSHOT_LOCATION)

    class PandasSnapshotExtenstionRefdata(pandas_snapshot_extention):
        @classmethod
        def dirname(cls, *, test_location: "PyTestLocation") -> str:
            return str(Path(test_location.filepath).parent.joinpath(refpath))
    return snapshot.use_extension(PandasSnapshotExtenstionRefdata)


@pytest.fixture
def snapshot_np(snapshot, tardis_ref_path, numpy_snapshot_extension):
    refpath = tardis_ref_path.joinpath(SNAPSHOT_LOCATION)
    class NumpySnapshotExtenstionRefdata(numpy_snapshot_extension):
            @classmethod
            def dirname(cls, *, test_location: "PyTestLocation") -> str:
                return str(Path(test_location.filepath).parent.joinpath(refpath))
    return snapshot.use_extension(NumpySnapshotExtenstionRefdata)

@pytest.fixture
def snapshot(snapshot, tardis_ref_path):
    refpath = tardis_ref_path.joinpath(SNAPSHOT_LOCATION)
    # TODO: try this with SingleFileSnapshotExtension
    class AmberSnapshotExtensionRefdata(AmberSnapshotExtension):
            @classmethod
            def dirname(cls, *, test_location: "PyTestLocation") -> str:
                return str(Path(test_location.filepath).parent.joinpath(refpath))
    return snapshot.use_extension(AmberSnapshotExtensionRefdata)

