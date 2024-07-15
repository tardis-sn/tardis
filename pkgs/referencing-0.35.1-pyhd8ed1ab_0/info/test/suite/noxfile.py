from pathlib import Path

import nox

ROOT = Path(__file__).parent

nox.options.sessions = []


def session(default=True, **kwargs):  # noqa: D103
    def _session(fn):
        if default:
            nox.options.sessions.append(kwargs.get("name", fn.__name__))
        return nox.session(**kwargs)(fn)

    return _session


@session()
def tests(session):
    """
    Run the sanity test suite to check the tests themselves.
    """
    session.install("jsonschema", "pytest")
    session.run("pytest", *session.posargs)


@session(tags=["style"])
def style(session):
    """
    Check Python code style in the sanity test suite.
    """
    session.install("ruff")
    session.run("ruff", "check", ROOT, __file__)
