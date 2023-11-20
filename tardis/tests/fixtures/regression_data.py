from pathlib import Path
import pytest
import re
import os


class RegressionData:
    def __init__(self, request) -> None:
        self.request = request
        tardis_ref_path = request.config.getoption("--tardis-refdata")
        self.tardis_ref_path = Path(
            os.path.expandvars(os.path.expanduser(tardis_ref_path))
        )
        self.enable_generate_reference = request.config.getoption(
            "--generate-reference"
        )

    @property
    def module_name(self):
        return self.request.node.module.__name__

    @property
    def test_name(self):
        return self.request.node.name

    @property
    def regression_data_fname_prefix(self):
        double_under = re.compile(r"[:\[\]{}]")
        no_space = re.compile(r'[,"\']')  # quotes and commas

        name = double_under.sub("__", self.test_name)
        name = no_space.sub("", name)
        return name

    @property
    def relative_regression_data_dir(self):
        return Path(self.module_name.replace(".", "/"))

    def check_data(self, data):
        full_fname_prefix = (
            self.relative_regression_data_dir
            / self.regression_data_fname_prefix
        )
        if self.enable_generate_reference:
            if hasattr(data, "to_hdf"):
                data.to_hdf(
                    full_fname_prefix.with_suffix(".h5"),
                )
            pytest.skip("Skipping test to generate reference data")

        else:
            if hasattr(data, "to_hdf"):
                ref_data = pd.read_hdf(
                    self.tardis_ref_path / f"{full_fname_prefix}.h5"
                )
            return ref_data


@pytest.fixture(scope="function")
def regression_data(request):
    return RegressionData(request)
