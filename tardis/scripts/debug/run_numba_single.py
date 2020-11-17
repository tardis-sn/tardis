from tardis import run_tardis
import numpy as np
from tardis.montecarlo.montecarlo_numba.base import montecarlo_main_loop
import os
import numba
import sys
import yaml


SEED = eval(sys.argv[1].split("=")[1])[0]

yaml_file, params = "tardis_example_single.yml", None

with open(yaml_file) as f:
    params = yaml.safe_load(f)

params["montecarlo"]["single_packet_seed"] = SEED

with open(yaml_file, "w") as f:
    yaml.safe_dump(params, f)

mdl = run_tardis(yaml_file)
