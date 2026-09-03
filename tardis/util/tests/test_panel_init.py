from unittest.mock import Mock

import pytest

import tardis.util.panel_init as panel_init


def test_notebook_uses_default_panel_comms(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    panel_extension = Mock()
    monkeypatch.setattr(panel_init.pn, "extension", panel_extension)

    panel_init.notebook()

    panel_extension.assert_called_once_with()
