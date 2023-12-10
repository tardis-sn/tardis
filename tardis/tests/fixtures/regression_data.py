import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tardis.io.util import HDFWriterMixin


class RegressionData:
    def __init__(self, request) -> None:
        self.request = request
        regression_data_path = request.config.getoption(
            "--tardis-regression-data"
        )
        if regression_data_path is None:
            pytest.skip("--tardis-regression-data was not specified")
        self.regression_data_path = Path(
            os.path.expandvars(os.path.expanduser(regression_data_path))
        )
        self.enable_generate_reference = request.config.getoption(
            "--generate-reference"
        )
        self.fname = f"{self.fname_prefix}.UNKNOWN_FORMAT"

    @property
    def module_name(self):
        return self.request.node.module.__name__

    @property
    def test_name(self):
        return self.request.node.name

    @property
    def fname_prefix(self):
        double_under = re.compile(r"[:\[\]{}]")
        no_space = re.compile(r'[,"\']')  # quotes and commas

        name = double_under.sub("__", self.test_name)
        name = no_space.sub("", name)
        return name

    @property
    def relative_regression_data_dir(self):
        relative_data_dir = Path(self.module_name.replace(".", "/"))
        if self.request.cls is not None:
            relative_data_dir /= HDFWriterMixin.convert_to_snake_case(
                self.request.cls.__name__
            )
        return relative_data_dir

    @property
    def absolute_regression_data_dir(self):
        return self.regression_data_path / self.relative_regression_data_dir

    @property
    def fpath(self):
        return self.absolute_regression_data_dir / self.fname

    def sync_dataframe(self, data, key="data"):
        self.fname = f"{self.fname_prefix}.h5"
        fpath = self.absolute_regression_data_dir / self.fname
        if self.enable_generate_reference:
            fpath.parent.mkdir(parents=True, exist_ok=True)
            data.to_hdf(
                fpath,
                key=key,
            )
            pytest.skip("Skipping test to generate reference data")
        else:
            return pd.read_hdf(fpath, key=key)

    def sync_ndarray(self, data):
        fpath = self.absolute_regression_data_dir / f"{self.fname_prefix}.npy"
        if self.enable_generate_reference:
            fpath.parent.mkdir(parents=True, exist_ok=True)
            np.save(fpath, data)
            pytest.skip("Skipping test to generate reference data")
        else:
            return np.load(fpath)

    def sync_str(self, data):
        fpath = self.absolute_regression_data_dir / f"{self.fname_prefix}.txt"
        if self.enable_generate_reference:
            fpath.parent.mkdir(parents=True, exist_ok=True)
            with fpath.open("w") as fh:
                fh.write(data)
            pytest.skip(
                f"Skipping test to generate regression_data {fpath} data"
            )
        else:
            with fpath.open("r") as fh:
                return fh.read()

    def sync_hdf_store(self, tardis_module):
        self.fname = f"{self.fname_prefix}.h5"
        if self.enable_generate_reference:
            self.fpath.parent.mkdir(parents=True, exist_ok=True)
            with pd.HDFStore(self.fpath, mode="w") as store:
                tardis_module.to_hdf(store, overwrite=True)
            pytest.skip(
                f"Skipping test to generate regression_data {self.fpath} data"
            )
        else:
            return pd.HDFStore(self.fpath, mode="r")


@pytest.fixture(scope="function")
def regression_data(request):
    return RegressionData(request)
