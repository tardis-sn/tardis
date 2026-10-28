"""Check configuration of continuum macro-atom worker threads."""

from typing import Any

import pytest
from jsonschema.exceptions import ValidationError

from tardis.io.configuration.config_reader import Configuration


@pytest.mark.parametrize("montecarlo_nthreads", [1, 4])
def test_continuum_threads_default_independently(
    tardis_config_verysimple: dict[str, Any], montecarlo_nthreads: int
) -> None:
    """Default continuum workers to one regardless of transport threads."""
    tardis_config_verysimple["montecarlo"]["nthreads"] = montecarlo_nthreads

    config = Configuration.from_config_dict(tardis_config_verysimple)

    assert config.plasma.continuum_interaction.get("nthreads") == 1
    assert config.montecarlo.nthreads == montecarlo_nthreads


@pytest.mark.parametrize("continuum_nthreads", [1, 3])
def test_continuum_threads_can_be_configured(
    tardis_config_verysimple: dict[str, Any], continuum_nthreads: int
) -> None:
    """Accept continuum worker counts independently of transport threads."""
    tardis_config_verysimple["montecarlo"]["nthreads"] = 4
    tardis_config_verysimple["plasma"]["continuum_interaction"] = {
        "nthreads": continuum_nthreads
    }

    config = Configuration.from_config_dict(tardis_config_verysimple)

    assert config.plasma.continuum_interaction.nthreads == continuum_nthreads
    assert config.montecarlo.nthreads == 4


@pytest.mark.parametrize("continuum_nthreads", [0, -1, 1.5, "2", True])
def test_continuum_threads_require_positive_integer(
    tardis_config_verysimple: dict[str, Any],
    continuum_nthreads: float | str | bool,
) -> None:
    """Reject worker counts that cannot specify a positive thread limit."""
    tardis_config_verysimple["plasma"]["continuum_interaction"] = {
        "nthreads": continuum_nthreads
    }

    with pytest.raises(ValidationError) as error:
        Configuration.from_config_dict(tardis_config_verysimple)

    assert list(error.value.absolute_path) == [
        "plasma",
        "continuum_interaction",
        "nthreads",
    ]
