from abc import ABCMeta, abstractmethod, abstractproperty
import pandas as pd


class BaseAtomicDataType(object):
    __metaclass__ = ABCMeta

    _data = None

    def hdf_name(self):
        raise NotImplementedError


    def load_hdf(self, hdf_store):
        self._data = self.hdf_name[self.hdf_name]

    @property
    def data(self):
        if self._data is None:
            raise ValueError('{0} is not available please load first'.format(
                self.__class__.__name__))
        else:
            return self._data

    @data.setter
    def data(self, value):
        self._data = value

    def save_hdf(self, hdf_store):
        hdf_store.put(self.hdf_name, self._data, format='table',
                      data_columns=True)


    @property
    def name(self):
        return getattr(self, '_name', self.hdf_name)


class Atoms(BaseAtomicDataType):
    hdf_name = 'atoms'


class Ions(BaseAtomicDataType):
    hdf_name = 'ions'

class Levels(BaseAtomicDataType):
    hdf_name = 'levels'

class Lines(BaseAtomicDataType):
    hdf_name = 'lines'
