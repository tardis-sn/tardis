import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.transport.montecarlo.packets.radiative_packet import InteractionType
from tardis.visualization.plot_util import (
    axis_label_in_latex,
    create_wavelength_mask,
    expand_species_list,
    extract_and_process_packet_data,
    get_mid_point_idx,
    get_spectrum_data,
    parse_species_list_util,
    to_rgb255_string,
)
from tardisbase.testing.regression_data.regression_data import PlotDataHDF


class TestPlotUtil:
    """Test utility functions used in plotting."""

    species_list = [["Si II", "Ca II", "C", "Fe I-V"]]
    packets_mode = ["real"]

    @pytest.mark.parametrize(
        "label_text,unit,only_text,expected",
        [
            (
                "Luminosity",
                u.erg / u.s,
                True,
                r"$\text{Luminosity}\,[\mathrm{erg\,s^{-1}}]$",
            ),
            (
                "Luminosity",
                u.erg / u.s,
                False,
                r"$Luminosity\,[\mathrm{erg\,s^{-1}}]$",
            ),
        ],
    )
    def test_axis_label_in_latex(self, label_text, unit, only_text, expected):
        result = axis_label_in_latex(label_text, unit, only_text)
        assert result == expected

    @pytest.mark.parametrize(
        "arr,expected_idx",
        [
            (np.array([1, 3, 7, 10]), 2),
            (np.array([10, 7, 3, 1]), 1),
            (np.array([1]), 0),
        ],
    )
    def test_get_mid_point_idx(self, arr, expected_idx):
        idx = get_mid_point_idx(arr)
        assert idx == expected_idx

    @pytest.mark.parametrize(
        "rgba,expected",
        [
            ((0.2, 0.4, 0.6, 1.0), "rgb(51, 102, 153)"),
            ((0.0, 0.0, 0.0, 1.0), "rgb(0, 0, 0)"),
        ],
    )
    def test_to_rgb255_string(self, rgba, expected):
        result = to_rgb255_string(rgba)
        assert result == expected

    @pytest.fixture(scope="module")
    def packet_data(self, simulation_simple_tracked):
        data = extract_and_process_packet_data(simulation_simple_tracked, "real")
        data["packets_df"]["last_interaction_type"] = data["packets_df"][
            "last_interaction_type"
        ].astype(str)
        data["packets_df_line_interaction"]["last_interaction_type"] = data[
            "packets_df_line_interaction"
        ]["last_interaction_type"].astype(str)
        return data

    def test_extract_and_process_packet_data_packets_df(
        self, packet_data, regression_data
    ):
        expected = regression_data.sync_dataframe(
            packet_data["packets_df"], key="packets_df"
        )
        pd.testing.assert_frame_equal(packet_data["packets_df"], expected)

    def test_extract_and_process_packet_data_line_interaction(
        self, packet_data, regression_data
    ):
        expected = regression_data.sync_dataframe(
            packet_data["packets_df_line_interaction"],
            key="packets_df_line_interaction",
        )
        pd.testing.assert_frame_equal(
            packet_data["packets_df_line_interaction"], expected
        )

    @pytest.mark.parametrize(
        "input_species, expected_output",
        [
            (
                ["Fe II", "Ca", "Si I-V"],
                ["Fe II", "Ca", "Si I", "Si II", "Si III", "Si IV", "Si V"],
            ),
            (
                ["Fe II", "Ca", "Si III-V"],
                ["Fe II", "Ca", "Si III", "Si IV", "Si V"],
            ),
            (["Fe 2"], None),
            (["Si 1-5"], None),
        ],
    )
    def test_expand_species_list(self, input_species, expected_output):
        if expected_output is None:
            with pytest.raises(ValueError):
                expand_species_list(input_species)
        else:
            assert expand_species_list(input_species) == expected_output

    @pytest.mark.parametrize(
        "input_species, expected_species_mapped, expected_species_list, expected_keep_colour, expected_full_species_list",
        [
            (
                ["Fe II", "Ca", "Si I-V"],
                {
                    (26, 1): [(26, 1)],
                    (20, 0): [(20, i) for i in range(20)],
                    (14, 0): [(14, 0)],
                    (14, 1): [(14, 1)],
                    (14, 2): [(14, 2)],
                    (14, 3): [(14, 3)],
                    (14, 4): [(14, 4)],
                },
                [(26, 1)]
                + [(20, i) for i in range(20)]
                + [(14, i) for i in range(5)],
                [20],
                ["Fe II", "Ca", "Si I", "Si II", "Si III", "Si IV", "Si V"],
            )
        ],
    )
    def test_parse_species_list_util(
        self,
        input_species,
        expected_species_mapped,
        expected_species_list,
        expected_keep_colour,
        expected_full_species_list,
    ):
        (
            species_mapped_result,
            species_list_result,
            keep_colour_result,
            full_species_list_result,
        ) = parse_species_list_util(input_species)
        assert set(species_mapped_result.keys()) == set(
            expected_species_mapped.keys()
        )
        for key in expected_species_mapped:
            assert key in species_mapped_result
            np.testing.assert_array_equal(
                species_mapped_result[key], expected_species_mapped[key]
            )
        np.testing.assert_array_equal(
            species_list_result, expected_species_list
        )
        np.testing.assert_array_equal(keep_colour_result, expected_keep_colour)
        assert full_species_list_result == expected_full_species_list

    def test_get_spectrum_data(self, simulation_simple_tracked):
        actual_data = get_spectrum_data("real", simulation_simple_tracked)
        packets_type = "spectrum_real_packets"

        expected_data = {
            "spectrum_delta_frequency": getattr(
                simulation_simple_tracked.spectrum_solver, packets_type
            ).delta_frequency,
            "spectrum_frequency_bins": getattr(
                simulation_simple_tracked.spectrum_solver, packets_type
            )._frequency,
            "spectrum_luminosity_density_lambda": getattr(
                simulation_simple_tracked.spectrum_solver, packets_type
            ).luminosity_density_lambda,
            "spectrum_wavelength": getattr(
                simulation_simple_tracked.spectrum_solver, packets_type
            ).wavelength,
        }

        for key, expected_value in expected_data.items():
            np.testing.assert_allclose(
                actual_data[key].value,
                expected_value.value,
            )

    @pytest.fixture(scope="module")
    def masked_packet_data(self, packet_data):
        mask = create_wavelength_mask(
            {"real": packet_data},
            "real",
            [3000, 9000] * u.AA,
            "packets_df",
            "nus",
        )
        return pd.DataFrame({"mask": mask})

    def test_create_wavelength_mask(
        self, masked_packet_data, regression_data
    ):
        expected = regression_data.sync_dataframe(masked_packet_data, key="mask")
        pd.testing.assert_frame_equal(masked_packet_data, expected)
