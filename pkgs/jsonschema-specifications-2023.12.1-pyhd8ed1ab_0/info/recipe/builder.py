import tempfile, subprocess, shutil, os, sys
from pathlib import Path

PY = sys.executable
SRC = Path(os.environ["SRC_DIR"]) / "src"

with tempfile.TemporaryDirectory() as td:
    tdp = Path(td)
    dest = tdp / SRC.name
    shutil.copytree(SRC, dest)
    rc = subprocess.call(
        [PY, "-m", "pip", "install", ".", "-vv", "--no-deps", "--no-build-isolation"],
        cwd=str(dest),
    )
    sys.exit(rc)
