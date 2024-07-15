#  tests for ipython-8.25.0-py312h06a4308_0 (this is a generated file);
print('===== testing package: ipython-8.25.0-py312h06a4308_0 =====');
print('running run_test.py');
#  --- run_test.py (begin) ---
import subprocess
import platform
import os
import sys
from pathlib import Path
import IPython

WIN = platform.system() == "Windows"
LINUX = platform.system() == "Linux"
PYPY = "__pypy__" in sys.builtin_module_names
PPC = "ppc" in platform.machine()

COV_THRESHOLD = os.environ.get("COV_THRESHOLD")

# Environment variable should be set in the meta.yaml
MIGRATING = eval(os.environ.get("MIGRATING", "None"))

PYTEST_SKIPS = ["decorator_skip", "pprint_heap_allocated"]
PYTEST_ARGS = [sys.executable, "-m", "pytest", "-vv"]

IGNORE_GLOBS = [
    "consoleapp.py",
    "external/*.py",
    "sphinxext/*.py",
    "terminal/console*.py",
    "terminal/pt_inputhooks/*.py",
    "utils/*.py",
]

PYTEST_ARGS += sum([[f"--ignore-glob", glob] for glob in IGNORE_GLOBS], [])

if WIN:
    pass
else:
    pass

if LINUX:
    PYTEST_SKIPS += ["system_interrupt"]

if len(PYTEST_SKIPS) == 1:
    PYTEST_ARGS += ["-k", f"not {PYTEST_SKIPS[0]}"]
elif PYTEST_SKIPS:
    PYTEST_ARGS += ["-k", f"""not ({" or ".join(PYTEST_SKIPS) })"""]

if __name__ == "__main__":
    print("Building on Windows?      ", WIN)
    print("Building on Linux?        ", LINUX)
    print("Building for PyPy?        ", PYPY)
    print("Running pytest with args")
    print(PYTEST_ARGS, flush=True)
    sys.exit(subprocess.call(PYTEST_ARGS, cwd=str(Path(IPython.__file__).parent)))#  --- run_test.py (end) ---

print('===== ipython-8.25.0-py312h06a4308_0 OK =====');
