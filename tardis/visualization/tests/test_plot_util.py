import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.tests.fixtures.regression_data import PlotDataHDF
from tardis.visualization.plot_util import (
    axis_label_in_latex,
    extract_and_process_packet_data,
    get_mid_point_idx,
    get_spectrum_data,
    parse_species_list_util,
    to_rgb255_string,
)


class TestPlotUtil:
    """Test utility functions used in plotting."""

    species_list = [["Si II", "Ca II", "C", "Fe I-V"]]
    packets_mode = ["real", "virtual"]

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

    @pytest.mark.parametrize("packets_mode", ["real", "virtual"])
    def test_extract_and_process_packet_data(
        self, simulation_simple, packets_mode
    ):
        actual_data = extract_and_process_packet_data(
            simulation_simple, packets_mode
        )

        transport_state = simulation_simple.transport.transport_state
        lines_df = (
            simulation_simple.plasma.atomic_data.lines.reset_index().set_index(
                "line_id"
            )
        )

        if packets_mode == "virtual":
            vpt = transport_state.vpacket_tracker
            expected_data = {
                "last_interaction_type": vpt.last_interaction_type,
                "last_line_interaction_in_id": vpt.last_interaction_in_id,
                "last_line_interaction_out_id": vpt.last_interaction_out_id,
                "last_line_interaction_in_nu": vpt.last_interaction_in_nu,
                "last_interaction_in_r": vpt.last_interaction_in_r,
                "nus": u.Quantity(vpt.nus, "Hz"),
                "energies": u.Quantity(vpt.energies, "erg"),
                "lambdas": u.Quantity(vpt.nus, "Hz").to(
                    "angstrom", u.spectral()
                ),
            }
        else:
            mask = transport_state.emitted_packet_mask
            nus = u.Quantity(
                transport_state.packet_collection.output_nus[mask], u.Hz
            )
            expected_data = {
                "last_interaction_type": transport_state.last_interaction_type[
                    mask
                ],
                "last_line_interaction_in_id": transport_state.last_line_interaction_in_id[
                    mask
                ],
                "last_line_interaction_out_id": transport_state.last_line_interaction_out_id[
                    mask
                ],
                "last_line_interaction_in_nu": transport_state.last_interaction_in_nu[
                    mask
                ],
                "last_interaction_in_r": transport_state.last_interaction_in_r[
                    mask
                ],
                "nus": nus,
                "energies": transport_state.packet_collection.output_energies[
                    mask
                ],
                "lambdas": nus.to("angstrom", u.spectral()),
            }

        expected_df = pd.DataFrame(expected_data)

        line_mask = (expected_df["last_interaction_type"] > -1) & (
            expected_df["last_line_interaction_in_id"] > -1
        )
        expected_df_line_interaction = expected_df.loc[line_mask].copy()
        expected_df_line_interaction["last_line_interaction_atom"] = (
            lines_df["atomic_number"]
            .iloc[expected_df_line_interaction["last_line_interaction_out_id"]]
            .to_numpy()
        )
        expected_df_line_interaction["last_line_interaction_species"] = (
            lines_df["atomic_number"]
            .iloc[expected_df_line_interaction["last_line_interaction_out_id"]]
            .to_numpy()
            * 100
            + lines_df["ion_number"]
            .iloc[expected_df_line_interaction["last_line_interaction_out_id"]]
            .to_numpy()
        )

        pd.testing.assert_frame_equal(
            actual_data["packets_df"].reset_index(drop=True),
            expected_df.reset_index(drop=True),
        )

        pd.testing.assert_frame_equal(
            actual_data["packets_df_line_interaction"].reset_index(drop=True),
            expected_df_line_interaction.reset_index(drop=True),
        )

    @pytest.fixture(scope="function")
    def generate_parse_species_hdf(self):
        (
            species_mapped_result,
            species_list_result,
            keep_colour_result,
            full_species_list,
        ) = parse_species_list_util(self.species_list[0])

        plot_object = np.array(
            [
                list(item)
                for sublist in species_mapped_result.values()
                for item in sublist
            ]
        )

        property_group = {
            "species_mapped": plot_object,
            "species_list": np.array(species_list_result),
            "keep_colour": np.array(keep_colour_result),
            "full_species_list": np.array(full_species_list),
        }

        return PlotDataHDF(**property_group)

    def test_parse_species_list_util(
        self, regression_data, generate_parse_species_hdf
    ):
        expected = regression_data.sync_hdf_store(generate_parse_species_hdf)

        for key in [
            "species_mapped",
            "species_list",
            "keep_colour",
            "full_species_list",
        ]:
            expected_values = expected.get("plot_data_hdf/" + key)
            actual_values = getattr(generate_parse_species_hdf, key)

            if key in ["species_list", "full_species_list"]:
                np.testing.assert_array_equal(expected_values, actual_values)
            else:
                np.testing.assert_allclose(expected_values, actual_values)
        expected.close()

    @pytest.mark.parametrize("packets_mode", ["real", "virtual"])
    def test_get_spectrum_data(self, simulation_simple, packets_mode):
        actual_data = get_spectrum_data(packets_mode, simulation_simple)
        packets_type = f"spectrum_{packets_mode}_packets"

        expected_data = {
            "spectrum_delta_frequency": getattr(
                simulation_simple.spectrum_solver, packets_type
            ).delta_frequency,
            "spectrum_frequency_bins": getattr(
                simulation_simple.spectrum_solver, packets_type
            )._frequency,
            "spectrum_luminosity_density_lambda": getattr(
                simulation_simple.spectrum_solver, packets_type
            ).luminosity_density_lambda,
            "spectrum_wavelength": getattr(
                simulation_simple.spectrum_solver, packets_type
            ).wavelength,
        }

        for key, expected_value in expected_data.items():
            np.testing.assert_allclose(
                actual_data[key].value,
                expected_value.value,
            )
