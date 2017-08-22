import os
import pygraphviz as pgv
import pytest

data_path = os.path.join('tardis', 'plasma', 'tests', 'data')


def test_write_dot(tmpdir, simulation_verysimple):
    fname = str(tmpdir.mkdir('test_dot').join('plasma.dot'))
    simulation_verysimple.plasma.write_to_dot(fname)
    actual = pgv.AGraph(fname).to_string()
    expected = pgv.AGraph(os.path.join(data_path, 'plasma_ref.dot')).to_string()
    assert actual == expected


def test_write_tex(tmpdir, simulation_verysimple):
    fname = str(tmpdir.mkdir('test_tex').join('plasma.tex'))
    simulation_verysimple.plasma.write_to_tex(fname)
    with open(fname, 'r') as fp1, open(os.path.join(data_path, 'plasma_ref.tex'), 'r') as fp2:
        assert fp1.readline() ==fp2.readline()
