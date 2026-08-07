import numpy as np
import numpy.testing as npt
import pytest

from tardis.simulation import Simulation


@pytest.mark.parametrize(
    "input_params,sliced",
    [
        ("scatter", False),
        ("macroatom", False),
        ("macroatom", True),
        ("downbranch", False),
        ("downbranch", True),
    ],
)
def test_opacity_state_to_numba(
    nb_simulation_verysimple: Simulation,
    input_params: str,
    sliced: bool,
) -> None:
    line_interaction_type = input_params
    plasma = nb_simulation_verysimple.plasma
    macro_atom_state = (
        None
        if line_interaction_type == "scatter"
        else nb_simulation_verysimple.macro_atom_state
    )
    actual = nb_simulation_verysimple.opacity_state.to_numba(
        macro_atom_state,
        line_interaction_type,
    )

    if sliced:
        index = slice(2, 5)
        actual = actual[index]
    else:
        index = ...

    npt.assert_allclose(
        actual.electron_density, plasma.electron_densities.values[index]
    )
    npt.assert_allclose(actual.line_list_nu, plasma.atomic_data.lines.nu.values)
    npt.assert_allclose(
        actual.tau_sobolev, plasma.tau_sobolevs.values[:, index]
    )
    if line_interaction_type == "scatter":
        empty = np.zeros(1, dtype=np.int64)
        npt.assert_allclose(
            actual.transition_probabilities, np.zeros((1, 1), dtype=np.float64)
        )
        npt.assert_allclose(actual.linf2macro_level_upper, empty)
        npt.assert_allclose(actual.macro_block_edge_index, empty)
        npt.assert_allclose(actual.transition_type, empty)
        npt.assert_allclose(actual.destination_level_id, empty)
        npt.assert_allclose(actual.transition_line_id, empty)
    else:
        npt.assert_allclose(
            actual.transition_probabilities,
            macro_atom_state.transition_probabilities.values[:, index],
        )
        npt.assert_allclose(
            actual.line2macro_level_upper,
            macro_atom_state.line2macro_level_upper.values,
        )
        npt.assert_allclose(
            actual.macro_block_edge_index,
            macro_atom_state.macro_block_edge_index,
        )
        npt.assert_allclose(
            actual.transition_type,
            macro_atom_state.transition_metadata.transition_type.values,
        )
        npt.assert_allclose(
            actual.destination_level_id,
            macro_atom_state.transition_metadata.destination_level_idx.values,
        )
        npt.assert_allclose(
            actual.transition_line_id,
            macro_atom_state.transition_metadata.transition_line_idx.values,
        )


def test_VPacketCollection_add_packet(verysimple_3vpacket_collection):
    assert verysimple_3vpacket_collection.length == 0

    nus = [3.0e15, 0.0, 1e15, 1e5]
    energies = [0.4, 0.1, 0.6, 1e10]
    initial_mus = [0.1, 0, 1, 0.9]
    initial_rs = [3e42, 4.5e45, 0, 9.0e40]
    last_interaction_in_nus = np.array(
        [3.0e15, 0.0, 1e15, 1e5], dtype=np.float64
    )
    last_interaction_in_rs = np.array(
        [3e42, 4.5e45, 0, 9.0e40], dtype=np.float64
    )
    last_interaction_types = np.array([1, 1, 3, 2], dtype=np.int64)
    last_interaction_in_ids = np.array([100, 0, 1, 1000], dtype=np.int64)
    last_interaction_out_ids = np.array([1201, 123, 545, 1232], dtype=np.int64)
    last_interaction_shell_ids = np.array([2, -1, 6, 0], dtype=np.int64)

    for (
        nu,
        energy,
        initial_mu,
        initial_r,
        last_interaction_in_nu,
        last_interaction_in_r,
        last_interaction_type,
        last_interaction_in_id,
        last_interaction_out_id,
        last_interaction_shell_id,
    ) in zip(
        nus,
        energies,
        initial_mus,
        initial_rs,
        last_interaction_in_nus,
        last_interaction_in_rs,
        last_interaction_types,
        last_interaction_in_ids,
        last_interaction_out_ids,
        last_interaction_shell_ids,
    ):
        verysimple_3vpacket_collection.add_packet(
            nu,
            energy,
            initial_mu,
            initial_r,
            last_interaction_in_nu,
            last_interaction_in_r,
            last_interaction_type,
            last_interaction_in_id,
            last_interaction_out_id,
            last_interaction_shell_id,
        )

    npt.assert_array_equal(
        verysimple_3vpacket_collection.nus[
            : verysimple_3vpacket_collection.idx
        ],
        nus,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.energies[
            : verysimple_3vpacket_collection.idx
        ],
        energies,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.initial_mus[
            : verysimple_3vpacket_collection.idx
        ],
        initial_mus,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.initial_rs[
            : verysimple_3vpacket_collection.idx
        ],
        initial_rs,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.last_interaction_in_nu[
            : verysimple_3vpacket_collection.idx
        ],
        last_interaction_in_nus,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.last_interaction_in_r[
            : verysimple_3vpacket_collection.idx
        ],
        last_interaction_in_rs,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.last_interaction_type[
            : verysimple_3vpacket_collection.idx
        ],
        last_interaction_types,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.last_interaction_in_id[
            : verysimple_3vpacket_collection.idx
        ],
        last_interaction_in_ids,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.last_interaction_out_id[
            : verysimple_3vpacket_collection.idx
        ],
        last_interaction_out_ids,
    )
    npt.assert_array_equal(
        verysimple_3vpacket_collection.last_interaction_shell_id[
            : verysimple_3vpacket_collection.idx
        ],
        last_interaction_shell_ids,
    )
    assert verysimple_3vpacket_collection.length == 9
