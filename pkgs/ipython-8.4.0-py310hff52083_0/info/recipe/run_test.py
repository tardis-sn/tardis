import subprocess
import platform
import os
import sys

WIN = platform.system() == "Windows"
LINUX = platform.system() == "Linux"
PYPY = "__pypy__" in sys.builtin_module_names
PPC = "ppc" in platform.machine()

COV_THRESHOLD = os.environ.get("COV_THRESHOLD")

# Environment variable should be set in the meta.yaml
MIGRATING = eval(os.environ.get("MIGRATING", "None"))

PYTEST_SKIPS = []
PYTEST_ARGS = [sys.executable, "-m", "pytest", "--pyargs", "IPython", "-vv"]

if WIN:
    pass
else:
    pass

if LINUX:
    PYTEST_SKIPS += ["system_interrupt"]

if PPC:
    PYTEST_SKIPS += ["ipython_dir_8", "audio_data"]

if COV_THRESHOLD is not None:
    PYTEST_ARGS += [
        "--cov", "IPython", "--no-cov-on-fail", "--cov-fail-under", COV_THRESHOLD,
        "--cov-report", "term-missing:skip-covered"
    ]

if len(PYTEST_SKIPS) == 1:
    PYTEST_ARGS += ["-k", f"not {PYTEST_SKIPS[0]}"]
elif PYTEST_SKIPS:
    PYTEST_ARGS += ["-k", f"""not ({" or ".join(PYTEST_SKIPS) })"""]

if __name__ == "__main__":
    print("Building on Windows?", WIN)
    print("Building on Linux?  ", LINUX)
    print("Building for PyPy?  ", PYPY)

    if MIGRATING:
        print("This is a migration, skipping test suite! Put it back later!", flush=True)
        sys.exit(0)
    else:
        print("Running pytest with args")
        print(PYTEST_ARGS, flush=True)
        sys.exit(subprocess.call(PYTEST_ARGS))
