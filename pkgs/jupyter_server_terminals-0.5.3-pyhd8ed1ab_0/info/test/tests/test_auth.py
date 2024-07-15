"""Tests for authorization"""
import asyncio
from typing import Dict

import pytest
from jupyter_server.auth.authorizer import Authorizer
from jupyter_server.auth.utils import HTTP_METHOD_TO_AUTH_ACTION, match_url_to_resource
from tornado.httpclient import HTTPClientError
from tornado.websocket import WebSocketHandler
from traitlets.config.loader import Config


class AuthorizerforTesting(Authorizer):
    # Set these class attributes from within a test
    # to verify that they match the arguments passed
    # by the REST API.
    permissions: Dict[str, str] = {}  # noqa: RUF012

    def normalize_url(self, path):
        """Drop the base URL and make sure path leads with a /"""
        base_url = self.parent.base_url
        # Remove base_url
        if path.startswith(base_url):
            path = path[len(base_url) :]
        # Make sure path starts with /
        if not path.startswith("/"):
            path = "/" + path
        return path

    def is_authorized(self, handler, user, action, resource):
        # Parse Request
        method = "WEBSOCKET" if isinstance(handler, WebSocketHandler) else handler.request.method
        url = self.normalize_url(handler.request.path)

        # Map request parts to expected action and resource.
        expected_action = HTTP_METHOD_TO_AUTH_ACTION[method]
        expected_resource = match_url_to_resource(url)

        # Assert that authorization layer returns the
        # correct action + resource.
        assert action == expected_action
        assert resource == expected_resource

        # Now, actually apply the authorization layer.
        return all(
            [
                action in self.permissions.get("actions", []),
                resource in self.permissions.get("resources", []),
            ]
        )


@pytest.fixture()
def jp_server_config():
    return Config(
        {
            "ServerApp": {
                "jpserver_extensions": {"jupyter_server_terminals": True},
                "authorizer_class": AuthorizerforTesting,
            }
        }
    )


@pytest.fixture()
def send_request(jp_fetch, jp_ws_fetch):
    """Send to Jupyter Server and return response code."""

    async def _(url, **fetch_kwargs):
        fetch = jp_ws_fetch if url.endswith("channels") or "/websocket/" in url else jp_fetch

        try:
            r = await fetch(url, **fetch_kwargs, allow_nonstandard_methods=True)
            code = r.code
        except HTTPClientError as err:
            code = err.code
        else:
            if fetch is jp_ws_fetch:
                r.close()

        print(code, url, fetch_kwargs)
        return code

    return _


HTTP_REQUESTS = [
    {
        "method": "POST",
        "url": "/api/terminals",
        "body": "",
    },
    {
        "method": "GET",
        "url": "/api/terminals",
    },
    {
        "method": "GET",
        "url": "/terminals/websocket/{term_name}",
    },
    {
        "method": "DELETE",
        "url": "/api/terminals/{term_name}",
    },
]

HTTP_REQUESTS_PARAMETRIZED = [(req["method"], req["url"], req.get("body")) for req in HTTP_REQUESTS]

# -------- Test scenarios -----------


@pytest.mark.parametrize("method, url, body", HTTP_REQUESTS_PARAMETRIZED)  # noqa: PT006
@pytest.mark.parametrize("allowed", (True, False))  # noqa: PT007
async def test_authorized_requests(
    request,
    io_loop,
    send_request,
    jp_serverapp,
    method,
    url,
    body,
    allowed,
):
    term_manager = jp_serverapp.web_app.settings["terminal_manager"]
    request.addfinalizer(lambda: io_loop.run_sync(term_manager.terminate_all))
    term_model = term_manager.create()
    term_name = term_model["name"]

    url = url.format(term_name=term_name)
    if allowed:
        # Create a server with full permissions
        permissions = {
            "actions": ["read", "write", "execute"],
            "resources": [
                "terminals",
            ],
        }
        expected_codes = {200, 201, 204, None}  # Websockets don't return a code
    else:
        permissions = {"actions": [], "resources": []}
        expected_codes = {403}
    jp_serverapp.authorizer.permissions = permissions

    while True:
        code = await send_request(url, body=body, method=method)
        if code == 404:
            await asyncio.sleep(1)
            continue
        assert code in expected_codes
        break
