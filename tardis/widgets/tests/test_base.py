import pytest
from tardis.widgets.base import BaseShellInfo, SimulationShellInfo, HDFShellInfo


@pytest.fixture(scope='class')
def base_shell_info(simulation_verysimple):
    return BaseShellInfo(simulation_verysimple.model.t_radiative,
                         simulation_verysimple.model.w,
                         simulation_verysimple.plasma.abundance,
                         simulation_verysimple.plasma.number_density,
                         simulation_verysimple.plasma.ion_number_density,
                         simulation_verysimple.plasma.level_number_density)


@pytest.fixture(scope='class')
def simulation_shell_info(simulation_verysimple):
    return SimulationShellInfo(simulation_verysimple)


@pytest.fixture(scope="class")
def hdf_shell_info(hdf_file_path, simulation_verysimple):
    simulation_verysimple.to_hdf(hdf_file_path)  # save sim at hdf_file_path
    return HDFShellInfo(hdf_file_path)


class TestBaseShellInfo:
    def test_shells_data(self, base_shell_info, simulation_verysimple):
        shells_data = base_shell_info.shells_data()
        assert shells_data.shape[0] == len(
            simulation_verysimple.model.t_radiative)

    def test_element_count_data(self, base_shell_info, simulation_verysimple):
        element_count_data = base_shell_info.element_count(1)
        assert element_count_data.shape[0] == simulation_verysimple.plasma.abundance.shape[0]


class TestSimulationShellInfo(TestBaseShellInfo):
    # Override the base_shell_info fixture to use value of simulation_shell_info fixture
    @pytest.fixture
    def base_shell_info(self, simulation_shell_info):
        return simulation_shell_info


class TestHDFShellInfo(TestBaseShellInfo):
    # Override the base_shell_info fixture to use value of hdf_shell_info fixture
    @pytest.fixture
    def base_shell_info(self, hdf_shell_info):
        return hdf_shell_info
