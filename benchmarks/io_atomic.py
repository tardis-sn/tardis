"""
Basic TARDIS Benchmark.
"""
from astropy import units as u
from astropy.tests.helper import assert_quantity_allclose
from asv_runner.benchmarks.mark import skip_benchmark

from benchmarks.benchmark_base import BenchmarkBase
from tardis import constants as const


# @skip_benchmark
class BenchmarkIoAtomic(BenchmarkBase):
    """
    Class to benchmark the atomic function.
    """

    def __init__(self):
        pass

    @property
    def basic_atom_data(self):
        return self.kurucz_atomic_data.atom_data

    @property
    def ionization_data(self):
        return self.kurucz_atomic_data.ionization_data

    @property
    def levels(self):
        return self.kurucz_atomic_data.levels

    @property
    def lines(self):
        return self.kurucz_atomic_data.lines

    def time_atom_data_basic_atom_data(self):
        assert self.basic_atom_data.loc[2, "symbol"] == "He"
        assert_quantity_allclose(
            self.basic_atom_data.at[2, "mass"] * u.Unit("g"), 4.002602 * const.u.cgs
        )

    def time_atom_data_ionization_data(self):
        assert_quantity_allclose(
            self.ionization_data.loc[(2, 1)] * u.Unit("erg"), 24.587387936 * u.Unit("eV")
        )

    def time_atom_data_levels(self):
        assert_quantity_allclose(
            u.Quantity(self.levels.at[(2, 0, 2), "energy"], u.Unit("erg")).to(
                u.Unit("cm-1"), equivalencies=u.spectral()
            ),
            166277.542 * u.Unit("cm-1"),
        )

    def time_atom_data_lines(self):
        assert_quantity_allclose(
            self.lines.loc[(2, 0, 0, 6), "wavelength_cm"].values[0] * u.Unit("cm"),
            584.335 * u.Unit("Angstrom"),
        )

    def time_atomic_reprepare(self):
        self.kurucz_atomic_data.prepare_atom_data(
            [14, 20],
            line_interaction_type="scatter",
            nlte_species=[],
            continuum_interaction_species=[],
        )
