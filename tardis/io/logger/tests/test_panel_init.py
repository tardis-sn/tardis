"""Tests for Panel communication initialization used by logger widgets."""

from unittest.mock import Mock

import pytest

from tardis.util import panel_init


def test_notebook_uses_default_panel_comms(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep standard notebooks on Panel's default communication backend."""
    extension = Mock()
    monkeypatch.setattr(panel_init.pn, "extension", extension)

    panel_init.notebook()

    extension.assert_called_once_with()


def test_vscode_uses_vscode_panel_comms(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Keep VS Code widgets on Panel's VS Code communication backend."""
    extension = Mock()
    monkeypatch.setattr(panel_init.pn, "extension", extension)

    panel_init.vscode()

    extension.assert_called_once_with(comms="vscode")
