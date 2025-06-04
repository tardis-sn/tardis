import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from tardis.tests.fixtures.regression_data import PlotDataHDF
from tardis.visualization.plot_util import process_line_interactions
from tardis.visualization.plot_util import process_line_interactions
from tardis.visualization.plot_util import (
    axis_label_in_latex,
    create_wavelength_mask,
    expand_species_list,
    extract_and_process_packet_data,
    get_mid_point_idx,
    get_spectrum_data_from_spectrum_solver,
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
            simulation_simple.transport.transport_state,
            simulation_simple.plasma,
            packets_mode=packets_mode,
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
    def test_expand_species_list(
        self, input_species, expected_output
    ):
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

    @pytest.mark.parametrize("packets_mode", ["real", "virtual"])
    def test_get_spectrum_data(self, simulation_simple, packets_mode):
        actual_data = get_spectrum_data_from_spectrum_solver(
            simulation_simple.spectrum_solver, packets_mode
        )
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

    @pytest.fixture(scope="module")
    def generate_masked_dataframe_hdf(self, simulation_simple):
        packet_data = {
            "real": extract_and_process_packet_data(
                transport_state=simulation_simple.transport.transport_state,
                plasma=simulation_simple.plasma,
                packets_mode="real",
            ),
            "virtual": extract_and_process_packet_data(
                transport_state=simulation_simple.transport.transport_state,
                plasma=simulation_simple.plasma,
                packets_mode="virtual",
            ),
        }
        masked_data = {
            mode: PlotDataHDF(
                masked_df=pd.DataFrame(create_wavelength_mask(
                    packet_data, mode, [3000, 9000] * u.AA, df_key="packets_df", column_name="nus"
                ))
            )
            for mode in ["real", "virtual"]
        }
        return masked_data

    @pytest.mark.parametrize("mode", ["real", "virtual"])
    def test_create_wavelength_mask(self, generate_masked_dataframe_hdf, regression_data, mode):
        expected = regression_data.sync_dataframe(generate_masked_dataframe_hdf[mode].masked_df, key=mode)
        actual = generate_masked_dataframe_hdf[mode].masked_df
        pd.testing.assert_frame_equal(actual, expected)
        def test_packet_energies_and_nus_copy():
            # Create a packets_df with packet_energies and packet_nus only

            original = pd.DataFrame(
                {
                    "packet_energies": [10.0, 20.0],
                    "packet_nus": [1e3, 2e3],
                    "last_interaction_type": [-1, -1],
                    "last_line_interaction_in_id": [-1, -1],
                    "last_line_interaction_out_id": [-1, -1],
                }
            )
            packet_data = {"packets_df": original.copy()}
            # lines_df can be empty since no line interactions are expected
            lines_df = pd.DataFrame(columns=["atomic_number", "ion_number"])

            # Should not raise or warn, just copy columns
            process_line_interactions(packet_data, lines_df)

            df = packet_data["packets_df"]
            # energies and nus columns should be created and match originals
            assert "energies" in df.columns
            assert "nus" in df.columns
            assert df["energies"].tolist() == original["packet_energies"].tolist()
            assert df["nus"].tolist() == original["packet_nus"].tolist()

            # No interactions → empty line-interaction DataFrame
            assert "packets_df_line_interaction" in packet_data
            assert packet_data["packets_df_line_interaction"].empty

        def test_process_line_interactions_warns_and_maps_atom_and_species():

            # lines_df mapping line_id 10 → atomic_number=1, ion_number=3
            lines_df = pd.DataFrame(
                {"atomic_number": [1], "ion_number": [3]}, index=[10]
            )

            # Create a packets_df with two rows:
            # row 0 should be selected (both last_interaction_type and in_id > -1)
            # row 1 filtered out
            packets_df = pd.DataFrame(
                {
                    "energies": [5.0, 6.0],
                    "nus": [5e3, 6e3],
                    "last_interaction_type": [1, -1],
                    "last_line_interaction_in_id": [7, -1],
                    "last_line_interaction_out_id": [10, 10],
                }
            )
            packet_data = {"packets_df": packets_df.copy()}

            # Expect DeprecationWarning because 'energies' is present
            with pytest.warns(DeprecationWarning):
                process_line_interactions(packet_data, lines_df)

            df_line = packet_data["packets_df_line_interaction"]
            # Only the first row should remain
            assert len(df_line) == 1
            row = df_line.reset_index(drop=True).iloc[0]

            # Verify atom and species columns
            assert row["last_line_interaction_atom"] == 1
            # species = atomic_number * 100 + ion_number = 1*100 + 3
            assert row["last_line_interaction_species"] == 103
