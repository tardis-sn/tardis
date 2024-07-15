#  tests for terminado-0.15.0-py310hff52083_0 (this is a generated file);
print('===== testing package: terminado-0.15.0-py310hff52083_0 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
""" run terminado tests with pytest, including platform- and python-based skips

    this is needed because `--pyargs` is not compatible with `-k` for
    function/method-based names
"""
import os
import sys
import pkgutil
import subprocess

platform = sys.platform
target_platform = os.environ["target_platform"]
py_major = sys.version_info[:2]
pypy = "__pypy__" in sys.builtin_module_names

loader = pkgutil.get_loader("terminado.tests")
test_path = os.path.dirname(loader.path)
pytest = [sys.executable, "-m", "pytest"]
pytest_args = [test_path, "-vv"]

if not pypy:
    pytest_args += ["--cov", "terminado", "--no-cov-on-fail"]


skips = []

# flaky tests
if platform != "linux":
    skips += [
        "basic_command",
        "max_terminals",
        "namespace",
        "single_process",
        "unique_processes"
    ]

if "aarch64" in target_platform:
    skips += ["max_terminals"]

if not skips:
    print("all tests will be run", flush=True)

elif len(skips) == 1:
    pytest_args += ["-k", "not {}".format(*skips)]
else:
    pytest_args += ["-k", "not ({})".format(" or ".join(skips))]

print("Final pytest args for", platform, target_platform, py_major)
print(" ".join([*pytest, *pytest_args]), flush=True)

# actually run the tests
sys.exit(subprocess.call([*pytest, *pytest_args]))
#  --- run_test.py (end) ---

print('===== terminado-0.15.0-py310hff52083_0 OK =====');
print("import: 'terminado'")
import terminado

