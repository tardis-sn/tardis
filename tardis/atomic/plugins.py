from abc import ABCMeta, abstractmethod, abstractproperty
import pandas as pd


class BaseAtomicDataPlugin(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def hdf_name(self):
        raise NotImplementedError


    @abstractmethod
    def load_hdf(self, hdf_store):
        raise NotImplementedError


    def save_hdf(self, hdf_store):
        hdf_store.put(self.hdf_name, self.data, format='table',
                      data_columns=True)


    @property
    def name(self):
        return getattr(self, '_name', self.hdf_name)


class Atoms():

    hdf_name = 'atoms'





class Ions(DB2HDF):

    hdf_name = 'ions'

    def __init__(self, **kwargs):
        self.exclude_species = kwargs.get('exclude_species', [])
        self.max_ionization_energy = kwargs.get('max_ionization_energy', np.inf)


    def load_sql(self):
        ion_data = self.atomic_db.session.query(Ion, Atom).join("atom").values(
            Atom.atomic_number, Ion.ion_number,
            Ion.ionization_energy)
        self.data = self._to_data_frame(ion_data)


    def to_hdf(self, file_or_buf):
        raise NotImplementedError()


class Levels(DB2HDF):

    hdf_name = 'levels'

    def __init__(self, exclude_species=[], max_ionization_energy=np.inf):
        self.exclude_species = exclude_species
        self.max_ionization_energy = max_ionization_energy

    def load_sql(self):
        self.data = None

    def to_hdf(self, file_or_buf):
        raise NotImplementedError()

class MyData(BaseAtomicDataPlugin):
    hdf_name = 'blah'
