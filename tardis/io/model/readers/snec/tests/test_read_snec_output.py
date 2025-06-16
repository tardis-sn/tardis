import astropy.units as u
import pytest

from tardis.io.model.readers.snec.snec_output import read_snec_output_xg
from tardis.io.model.readers.snec.xg_files import XGData
from tardis.tests.fixtures.regression_data import RegressionData

@pytest.fixture
def regression_test_snec_dir(regression_data: RegressionData):
    return regression_data.regression_data_path / "testdata" / "MESA_STIR_MESA_SNEC"

def test_read_snec_output_xg(regression_test_snec_dir):
    xg_data = read_snec_output_xg(regression_test_snec_dir, show_progress=False)
    assert isinstance(xg_data, XGData)
    assert len(xg_data.timestamps) == 1001
    assert len(xg_data.data_blocks) == 1001
    assert xg_data.timestamps.unit == u.s
    expected_columns = {'radius', 'enclosed_mass', 'mass', 'vel', 'rho', 'temp', 'logT',
                       'tau', 'lum', 'p_rad', 'press', 'E_shell', 'Ni_deposit_function',
                       'ye', 'free_electron_frac', 'photosphere_tracer', 'time_diff',
                       'delta_time', 'time_exp', 'Q', 'kappa', 'kappa_table', 'eps',
                       'logR_op', 'cs2', 'H_1', 'H_2', 'He_1', 'He_2', 'He_3'}
    actual_columns = set(xg_data.data_blocks[0].columns)
    assert expected_columns.symmetric_difference(actual_columns) == set()
