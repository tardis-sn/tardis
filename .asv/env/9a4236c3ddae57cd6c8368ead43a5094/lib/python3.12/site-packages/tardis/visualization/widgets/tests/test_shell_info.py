import pytest
import numpy as np
import pandas.testing as pdt

from tardis.tests.test_util import monkeysession
from tardis.visualization.widgets.shell_info import (
    BaseShellInfo,
    SimulationShellInfo,
    HDFShellInfo,
    ShellInfoWidget,
)


@pytest.fixture(scope="class")
def base_shell_info(simulation_verysimple):
    return BaseShellInfo(
        simulation_verysimple.simulation_state.t_radiative,
        simulation_verysimple.simulation_state.dilution_factor,
        simulation_verysimple.plasma.abundance,
        simulation_verysimple.plasma.number_density,
        simulation_verysimple.plasma.ion_number_density,
        simulation_verysimple.plasma.level_number_density,
    )


@pytest.fixture(scope="class")
def simulation_shell_info(simulation_verysimple):
    return SimulationShellInfo(simulation_verysimple)


@pytest.fixture(scope="class")
def hdf_shell_info(hdf_file_path, simulation_verysimple):
    simulation_verysimple.to_hdf(
        hdf_file_path, overwrite=True
    )  # save sim at hdf_file_path
    return HDFShellInfo(hdf_file_path)


class TestBaseShellInfo:
    def test_shells_data(self, base_shell_info, simulation_verysimple):
        shells_data = base_shell_info.shells_data()
        assert shells_data.shape == (
            len(simulation_verysimple.simulation_state.t_radiative),
            2,
        )
        assert np.allclose(
            shells_data.iloc[:, 0].map(np.float64),
            simulation_verysimple.simulation_state.t_radiative.value,
        )
        assert np.allclose(
            shells_data.iloc[:, 1].map(np.float64),
            simulation_verysimple.simulation_state.dilution_factor,
        )

    @pytest.mark.parametrize("shell_num", [1, 20])
    def test_element_count_data(
        self, base_shell_info, simulation_verysimple, shell_num
    ):
        element_count_data = base_shell_info.element_count(1)
        assert element_count_data.shape == (
            len(simulation_verysimple.plasma.abundance[shell_num - 1]),
            2,
        )
        assert np.allclose(
            element_count_data.iloc[:, -1].map(np.float64),
            simulation_verysimple.plasma.abundance[shell_num - 1],
        )

    @pytest.mark.parametrize(("atomic_num", "shell_num"), [(12, 1), (20, 20)])
    def test_ion_count_data(
        self, base_shell_info, simulation_verysimple, atomic_num, shell_num
    ):
        ion_count_data = base_shell_info.ion_count(atomic_num, shell_num)
        sim_ion_number_density = (
            simulation_verysimple.plasma.ion_number_density[shell_num - 1].loc[
                atomic_num
            ]
        )
        sim_element_number_density = (
            simulation_verysimple.plasma.number_density.loc[
                atomic_num, shell_num - 1
            ]
        )
        assert ion_count_data.shape == (len(sim_ion_number_density), 2)
        assert np.allclose(
            ion_count_data.iloc[:, -1].map(np.float64),
            sim_ion_number_density / sim_element_number_density,
        )

    @pytest.mark.parametrize(
        ("ion_num", "atomic_num", "shell_num"), [(2, 12, 1), (3, 20, 20)]
    )
    def test_level_count_data(
        self,
        base_shell_info,
        simulation_verysimple,
        ion_num,
        atomic_num,
        shell_num,
    ):
        level_count_data = base_shell_info.level_count(
            ion_num, atomic_num, shell_num
        )
        sim_level_number_density = (
            simulation_verysimple.plasma.level_number_density[
                shell_num - 1
            ].loc[atomic_num, ion_num]
        )
        sim_ion_number_density = (
            simulation_verysimple.plasma.ion_number_density[shell_num - 1].loc[
                atomic_num, ion_num
            ]
        )
        assert level_count_data.shape == (len(sim_level_number_density), 1)
        assert np.allclose(
            level_count_data.iloc[:, 0].map(np.float64),
            sim_level_number_density / sim_ion_number_density,
        )


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


class TestShellInfoWidget:
    # Indices of each table to select for testing
    select_shell_num = 4
    select_atomic_num = 12
    select_ion_num = 3

    @pytest.fixture(scope="class")
    def shell_info_widget(self, base_shell_info, monkeysession):
        shell_info_widget = ShellInfoWidget(base_shell_info)
        monkeysession.setattr(
            "tardis.visualization.widgets.shell_info.is_notebook", lambda: True
        )
        # To attach event listeners to table widgets of shell_info_widget
        _ = shell_info_widget.display()
        return shell_info_widget

    def test_selection_on_shells_table(
        self, base_shell_info, shell_info_widget
    ):
        shell_info_widget.shells_table.change_selection([self.select_shell_num])

        expected_element_count = base_shell_info.element_count(
            self.select_shell_num
        )
        pdt.assert_frame_equal(
            expected_element_count, shell_info_widget.element_count_table.df
        )

        expected_ion_count = base_shell_info.ion_count(
            expected_element_count.index[0], self.select_shell_num
        )
        pdt.assert_frame_equal(
            expected_ion_count, shell_info_widget.ion_count_table.df
        )

        expected_level_count = base_shell_info.level_count(
            expected_ion_count.index[0],
            expected_element_count.index[0],
            self.select_shell_num,
        )
        pdt.assert_frame_equal(
            expected_level_count, shell_info_widget.level_count_table.df
        )

    def test_selection_on_element_count_table(
        self, base_shell_info, shell_info_widget
    ):
        shell_info_widget.element_count_table.change_selection(
            [self.select_atomic_num]
        )

        expected_ion_count = base_shell_info.ion_count(
            self.select_atomic_num, self.select_shell_num
        )
        pdt.assert_frame_equal(
            expected_ion_count, shell_info_widget.ion_count_table.df
        )

        expected_level_count = base_shell_info.level_count(
            expected_ion_count.index[0],
            self.select_atomic_num,
            self.select_shell_num,
        )
        pdt.assert_frame_equal(
            expected_level_count, shell_info_widget.level_count_table.df
        )

    def test_selection_on_ion_count_table(
        self, base_shell_info, shell_info_widget
    ):
        shell_info_widget.ion_count_table.change_selection(
            [self.select_ion_num]
        )

        expected_level_count = base_shell_info.level_count(
            self.select_ion_num, self.select_atomic_num, self.select_shell_num
        )
        pdt.assert_frame_equal(
            expected_level_count, shell_info_widget.level_count_table.df
        )
