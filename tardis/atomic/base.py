from tardis.atomic.plugins import BaseAtomicDataType
import pandas as pd

class AtomicData(object):

    all_data_types = BaseAtomicDataType.__subclasses__()

    @classmethod
    def from_hdf(cls, fname):
        return cls(pd.HDFStore(fname, mode='r'))

    def __init__(self, hdf_store, data_types=self.all_data_types):
        self.data_types = [data_type()
                           for data_type in data_types]

        self.name2data_type = {item.name:item for item in self.data_types}



    def __getattr__(self, item):
        if item in self.name2data_type:
            return self.name2date_type[item].value
        else:
            super(AtomicData, self).__getattribute__(item)


