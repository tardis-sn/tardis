"""Launch the local single-shell resonance explorer."""

from __future__ import annotations

import argparse

import panel as pn

from .app import build_app


def main() -> None:
    """Serve one explorer per browser session on localhost only."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, default=5006, help="Local TCP port")
    parser.add_argument(
        "--no-browser",
        action="store_true",
        help="Start the server without opening a browser",
    )
    arguments = parser.parse_args()
    pn.serve(
        build_app,
        address="127.0.0.1",
        port=arguments.port,
        websocket_origin=[f"127.0.0.1:{arguments.port}"],
        show=not arguments.no_browser,
        resources="inline",
        title="Single-shell resonance explorer",
    )


if __name__ == "__main__":
    main()
