import tardis
from tardis import run_tardis
from tardis.io.atom_data.util import download_atom_data
import os

import tempfile

#DATA_PATH = os.path.join(tardis.__path__[0], 'plasma', 'tests', 'data')

def run_simulation(yml_file):
#   automatically download atom data
    download_atom_data('kurucz_cd23_chianti_H_He')
#   running simulation and storing in variable
    sim = run_tardis(yml_file)
    return sim

def test_write_to_dot():
    #yml_file = os.path.join(DATA_PATH, 'write_to_dot_test.yml')
    sim =  run_simulation('write_to_dot_test.yml')
    fname = "temp_graph.tex"
    try:
#       running the write_to_tex method to write a tex file from output of write_to_dot.
        sim.plasma.write_to_tex(fname)
    except:
#       the dot file couldn't be converted to tex.
        raise "dot file from plasma.write_to_dot couldn't be parsed and converted to latex."
    os.remove(fname)