import pytest
import logging
import os
import pandas as pd
import numpy as np

from tardis.io.config_reader import Configuration
from tardis.simulation import Simulation
from tardis.io.logger.logger import LOGGING_LEVELS
from tardis import run_tardis


def test_logging_simulation(atomic_data_fname, caplog):
    """
    Testing the logs for simulations runs
    """
    config = Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple.yml"
    )
    config["atom_data"] = atomic_data_fname

    simulation = Simulation.from_config(config)

    simulation.run()

    for record in caplog.records:
        assert record.levelno >= logging.INFO


# Testing Configuration of Logger via run_tardis() Function
@pytest.mark.parametrize(
    ["log_level", "specific_log_level"],
    [
        ("Info", False),
        ("INFO", False),
        ("INFO", True),
        ("DEBUG", False),
        ("DEBUG", True),
        ("WARNING", True),
        ("ERROR", True),
        ("CRITICAL", True),
        ("NOTSET", False),
    ],
)
class TestSimulationLogging:
    """
    Class implemented for testing the logging configuration available via run_tardis()
    Tests Functional Arguments : log_level & specific
    Tests YAML Parameters : logging_level & specific_logging
    """

    def test_logging_config(
        self, atomic_data_fname, caplog, log_level, specific_log_level
    ):
        config = Configuration.from_yaml(
            "tardis/io/tests/data/tardis_configv1_verysimple_logger.yml"
        )
        config["atom_data"] = atomic_data_fname

        caplog.clear()
        run_tardis(
            config=config,
            log_level=log_level,
            specific_log_level=specific_log_level,
        )
        for record in caplog.records:
            if specific_log_level == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]

    def test_logging_config_yaml(
        self, atomic_data_fname, caplog, log_level, specific_log_level
    ):
        config = Configuration.from_yaml(
            "tardis/io/tests/data/tardis_configv1_verysimple_logger.yml"
        )
        config["atom_data"] = atomic_data_fname
        config["debug"]["log_level"] = log_level
        config["debug"]["specific_log_level"] = specific_log_level

        caplog.clear()
        run_tardis(config=config)
        for record in caplog.records:
            if specific_log_level == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]

    def test_logging_both_specified(
        self, atomic_data_fname, caplog, log_level, specific_log_level
    ):
        config = Configuration.from_yaml(
            "tardis/io/tests/data/tardis_configv1_verysimple_logger.yml"
        )
        config["atom_data"] = atomic_data_fname
        config["debug"]["log_level"] = log_level
        config["debug"]["specific_log_level"] = specific_log_level

        caplog.clear()
        run_tardis(
            config=config,
            log_level=log_level,
            specific_log_level=specific_log_level,
        )
        for record in caplog.records:
            if specific_log_level == True:
                assert record.levelno == LOGGING_LEVELS[log_level.upper()]
            else:
                assert record.levelno >= LOGGING_LEVELS[log_level.upper()]


@pytest.fixture
def config():
    return Configuration.from_yaml(
        "tardis/io/tests/data/tardis_configv1_verysimple_tracking.yml"
    )


@pytest.fixture
def tracker_ref_path(tardis_ref_path):
    return os.path.abspath(os.path.join(tardis_ref_path, "rpacket_tracking.h5"))


@pytest.fixture
def tracking_refdata(
    config, atomic_data_fname, tracker_ref_path, generate_reference
):
    config["atom_data"] = atomic_data_fname

    simulation = Simulation.from_config(config)
    simulation.run()

    track_df = simulation.runner.rpacket_tracker
    key = "tracking"

    if not generate_reference:
        return simulation
    else:
        track_df.to_hdf(tracker_ref_path, key=key, mode="w")
        pytest.skip("Reference data was generated during this run.")


@pytest.fixture
def read_comparison_refdata(tracker_ref_path):
    return pd.read_hdf(tracker_ref_path)


# @pytest.mark.parametrize(
#     ["no_of_packets", "initial_seed", "last_seed", "iterations"],
#     [(1200, 2850180890, 2683780343, 3)],
# )
def test_tracking_dataframe(
    config,
    tracking_refdata,
    # no_of_packets,
    # initial_seed,
    # last_seed,
    # iterations,
):
    sim = tracking_refdata

    # Initial Test to check if the data frame is generated or not
    assert config["montecarlo"]["tracking"]["track_rpacket"] == True
    assert isinstance(sim.runner.rpacket_tracker, pd.DataFrame)
    # # assert (
    # #     len(sim.runner.rpacket_tracker["Packet Seed"].unique()) == no_of_packets
    # # )
    # assert sim.runner.rpacket_tracker["Packet Seed"].iloc[0] == initial_seed
    # assert sim.runner.rpacket_tracker["Packet Seed"].iloc[-1] == last_seed
    # assert len(sim.runner.rpacket_tracker["Iteration"].unique()) == iterations


def test_compare_dataframe(
    tracking_refdata,
    read_comparison_refdata,
):
    sim = tracking_refdata
    comparison_df = read_comparison_refdata

    pd.testing.assert_frame_equal(
        sim.runner.rpacket_tracker,
        comparison_df,
        check_dtype=True,
        check_column_type=True,
        check_exact=True,
    )

    assert isinstance(comparison_df, pd.DataFrame)
    assert len(comparison_df["Packet Seed"].unique()) == len(
        sim.runner.rpacket_tracker["Packet Seed"].unique()
    )
    assert len(comparison_df["Iteration"].unique()) == len(
        sim.runner.rpacket_tracker["Iteration"].unique()
    )


def test_parallel_dataframe(
    config,
    atomic_data_fname,
    read_comparison_refdata,
):
    comparison_df = read_comparison_refdata

    config["atom_data"] = atomic_data_fname
    config["montecarlo"]["nthreads"] = 3

    sim = Simulation.from_config(config)
    sim.run()

    assert isinstance(sim.runner.rpacket_tracker, pd.DataFrame)
    assert len(comparison_df["Packet Seed"].unique()) == len(
        sim.runner.rpacket_tracker["Packet Seed"].unique()
    )
    assert len(comparison_df["Iteration"].unique()) == len(
        sim.runner.rpacket_tracker["Iteration"].unique()
    )

    config["montecarlo"]["nthreads"] = 1
