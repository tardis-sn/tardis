import filecmp
import os
import pygraphviz as pgv


data_path = os.path.join('tardis', 'plasma', 'tests', 'data')


def test_write_dot(tmpdir, simulation_verysimple):
    fname = str(tmpdir.mkdir('test').join('plasma.dot'))
    simulation_verysimple.plasma.write_to_dot(fname)
    assert pgv.AGraph(fname) == pgv.AGraph(
        os.path.join(data_path, 'plasma_ref.dot'))


def test_write_tex(tmpdir, simulation_verysimple):
    fname = str(tmpdir.mkdir('test').join('plasma.tex'))
    simulation_verysimple.plasma.write_to_dot(fname)
    assert filecmp.cmp(fname, os.path.join(data_path, 'plasma_ref.tex'))
