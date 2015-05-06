from abc import ABCMeta, abstractmethod, abstractproperty
import pandas as pd


class BaseAtomicData(object):
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


class Atoms(BaseAtomicData):
    hdf_name = 'atoms'

class Ions(BaseAtomicData):
    hdf_name = 'ions'

class Levels(BaseAtomicData):
    hdf_name = 'levels'

class Lines(BaseAtomicData):
    hdf_name = 'lines'
