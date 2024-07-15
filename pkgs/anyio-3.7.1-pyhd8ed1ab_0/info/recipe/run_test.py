import os
import sys
from pathlib import Path
import subprocess

HERE = Path(__file__).parent

print("preflighting uvloop...", flush=True)
try:
    import uvloop
    print(f"... uvloop {uvloop} _appears_ to be working!")
except Exception as err:
    print(f"... uvloop import error, ignoring:\n{err}")

COV_THRESHOLD = {
    "win32": "59",
    "linux": "86",
    "darwin": "65"
}.get(sys.platform)

SKIPS = [
    "connection_refused",
    "happy_eyeballs",
    "ipv6",
    "block_device",
]

PYTEST_ARGS = [
    "-vv",
    "--cov", os.environ["PKG_NAME"],
    "--cov-fail-under", COV_THRESHOLD,
    "--cov-report", "term-missing:skip-covered",
    "--no-cov-on-fail",
    *(["-k", f"""not ({" or ".join(SKIPS)})"""] if SKIPS else [])
]

print(">>> pytest", " ".join(PYTEST_ARGS), flush=True)

sys.exit(subprocess.call(["pytest", *PYTEST_ARGS], cwd=str(HERE / "src")))
