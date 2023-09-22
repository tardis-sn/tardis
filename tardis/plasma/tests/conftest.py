from pathlib import Path
import pytest
from syrupy.location import PyTestLocation

MODULE_NAME = "plasma"

@pytest.fixture
def snapshot_pd(snapshot, tardis_ref_path, pandas_snapshot_extention):
    refpath = tardis_ref_path.joinpath(MODULE_NAME)

    class PandasSnapshotExtenstionRefdata(pandas_snapshot_extention):
        @classmethod
        def dirname(cls, *, test_location: "PyTestLocation") -> str:
            return str(Path(test_location.filepath).parent.joinpath(refpath))
    return snapshot.use_extension(PandasSnapshotExtenstionRefdata)


@pytest.fixture
def snapshot_np(snapshot, tardis_ref_path, numpy_snapshot_extension):
    refpath = tardis_ref_path.joinpath(MODULE_NAME)
    class NumpySnapshotExtenstionRefdata(numpy_snapshot_extension):
            @classmethod
            def dirname(cls, *, test_location: "PyTestLocation") -> str:
                return str(Path(test_location.filepath).parent.joinpath(refpath))
    return snapshot.use_extension(NumpySnapshotExtenstionRefdata)

