print("Comparing recipe spec and nbconvert for pandoc:")
import sys
import subprocess

print("... PANDOC VERSION       ", flush=True)

print(subprocess.check_output(["pandoc", "--version"]).decode("utf-8"), flush=True)

recipe_spec = tuple(sys.argv[1:3])

print("... RECIPE PANDOC SPEC   ", recipe_spec, flush=True)

from nbconvert.utils import pandoc
pandoc_spec = (f">={pandoc._minimal_version}", f"<{pandoc._maximal_version}")

print("... NBCONVERT PANDOC SPEC", pandoc_spec, flush=True)

assert recipe_spec == pandoc_spec
