__author__ = 'wkerzend'

import os
import numpy as np
from astropy import constants, table, units

from numpy import testing
from tardis import plasma, atomic

atom_model = atomic.AtomData.from_hdf5()


def pytest_generate_tests(metafunc):
    idlist = []
    argvalues = []
    for scenario in metafunc.cls.scenarios:
        idlist.append(scenario[0])
        items = scenario[1].items()
        argnames = [x[0] for x in items]
        argvalues.append(([x[1] for x in items]))
    metafunc.parametrize(argnames, argvalues, ids=idlist, scope="class")


def convert_nist_energy(s):
    return float(s.strip().strip('?[]'))


def convert_nist_j(s):
    if '/' in s:
        numerator, denominator = s.split('/')
        return float(numerator) / float(denominator)
    else:
        return float(s)


def read_nist_data(fname):
    data = np.genfromtxt(fname, skip_header=3, delimiter='|', usecols=(2, 3),
        converters={2: convert_nist_j, 3: convert_nist_energy}, names=('j', 'energy'))
    data = data[~np.isnan(data['j'])]
    data = table.Table(data)
    data['energy'].units = units.Unit('eV')
    data['energy'].convert_units_to('erg')

    data = data[~np.isnan(data['energy'])]
    return data

tests_data_dir = os.path.join('tardis', 'tests', 'data')

nist_levels = []

for ion_number in xrange(1, 4):
    nist_levels.append(read_nist_data(os.path.join(tests_data_dir, 'nist_si%d.dat' % ion_number)))

scenarios = []
scenarios.append(('t1000', {'t_rad': 1000}))
scenarios.append(('t5000', {'t_rad': 5000}))
scenarios.append(('t10000', {'t_rad': 10000}))
scenarios.append(('t15000', {'t_rad': 15000}))
scenarios.append(('t20000', {'t_rad': 20000}))

class TestLTEPlasmaCalculations:
    scenarios = scenarios

    def test_lte_partition_functions(self, t_rad):
        plasma_model = plasma.LTEPlasma({'Si': 1.}, t_rad, 1e-8, atom_model)

        beta_rad = 1 / (constants.cgs.k_B * t_rad)

        testing.assert_almost_equal(plasma_model.beta_rad, beta_rad)
        self.nist_partitions = []
        for i, (atomic_number, ion_number, partition_function) in enumerate(plasma_model.partition_functions):
            nist_level = nist_levels[ion_number].__array__()
            nist_partition = np.sum((2 * nist_level['j'] + 1) * np.exp(-beta_rad * nist_level['energy']))
            self.nist_partitions.append(nist_partition)
            partition_delta = abs(nist_partition - partition_function) / nist_partition
            assert partition_delta < 0.05

    def test_lte_saha_calculation(self, t_rad):
        assert len(self.nist_partitions) == 0





