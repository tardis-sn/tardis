import pytest
from traitlets.config.loader import Config


@pytest.fixture()
def jp_server_config():
    return Config({"ServerApp": {"terminals_enabled": False}})


async def test_not_enabled(jp_configurable_serverapp):
    assert jp_configurable_serverapp().terminals_enabled is False
    assert jp_configurable_serverapp().web_app.settings["terminals_available"] is False
    assert "terminal_manager" not in jp_configurable_serverapp().web_app.settings
