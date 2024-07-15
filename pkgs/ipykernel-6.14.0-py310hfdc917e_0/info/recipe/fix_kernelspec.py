import os
from pathlib import Path
import sys
import json

prefix = Path(os.environ["PREFIX"])
kernel = "python{}".format(sys.version_info[0])
spec_path = prefix / "share" / "jupyter" / "kernels" / kernel / "kernel.json"
posix_exe = Path(sys.executable).as_posix()

print("Rewriting kernelspec at:\n\t{}".format(spec_path))

raw_spec = spec_path.read_text()

print(raw_spec)

spec = json.loads(raw_spec)

print("Kernel python was:\n\t{}".format(spec["argv"][0]))

if spec["argv"][0] == posix_exe:
    print("Path is fine")
else:
    print("Rewriting kernel python with:\n\t{}".format(posix_exe))
    spec["argv"][0] = posix_exe
    spec_path.write_text(json.dumps(spec, indent=2))
