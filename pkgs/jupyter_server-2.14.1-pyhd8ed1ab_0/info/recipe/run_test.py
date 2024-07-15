import sys
import subprocess

# set in script_env
COV_THRESHOLD = 70

# only really relevant when tested by upstreams, e.g. jupyter_core
PYPY = "__pypy__" in sys.builtin_module_names

# base opinions
PYTEST_ARGS = ["pytest", "-vv", "--tb=long", "--color=yes"]

# known to fail in conda-forge CI
SKIPS = [
    # delete tests fail if /tmp and /home/.../trash are on different devices
    "delete",
    # ???
    "merge_config",
]

# -k doesn't like singletons in parentheses
SKIP_JOINED = SKIPS[0] if len(SKIPS) == 1 else f"""({" or ".join(SKIPS)})"""
PYTEST_ARGS += ["-k", f"not {SKIP_JOINED}"]

# coverage works poorly on pypy
if not PYPY:
    print("not on pypy, testing with coverage")
    PYTEST_ARGS += [
        "--cov=jupyter_server",
        "--cov-report=term-missing:skip-covered",
        "--no-cov-on-fail",
    ]
    if COV_THRESHOLD:
        PYTEST_ARGS += [f"--cov-fail-under={COV_THRESHOLD}"]

if __name__ == "__main__":
    print(">>>", "\t".join(PYTEST_ARGS))
    # finally run
    sys.exit(subprocess.call(PYTEST_ARGS))
