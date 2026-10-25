"""Local server configuration needed by an actual browser session."""

import sys

import panel as pn
import pytest

from tardis.visualization.widgets.resonance_explorer.__main__ import main


def test_local_server_accepts_its_own_websocket_origin(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Claim: Firefox at 127.0.0.1 can open the Bokeh document socket."""
    options: dict[str, object] = {}

    def capture_server(*_args: object, **kwargs: object) -> None:
        options.update(kwargs)

    monkeypatch.setattr(pn, "serve", capture_server)
    monkeypatch.setattr(
        sys, "argv", ["resonance_explorer", "--no-browser", "--port", "5118"]
    )
    main()
    assert options["address"] == "127.0.0.1"
    assert options["websocket_origin"] == ["127.0.0.1:5118"]
