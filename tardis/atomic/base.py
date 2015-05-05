from tardis.atomic.plugins import BaseAtomicDataPlugin
import pandas as pd

class AtomicData(object):

    plugins = BaseAtomicDataPlugin.__subclasses__()

    @classmethod
    def from_hdf(cls, fname):
        return cls(pd.HDFStore(fname, mode='r'))

    def __init__(self, hdf_store):
        self.plugins = [plugin()
                        for plugin in BaseAtomicDataPlugin.__subclasses__()]
        self.name2plugin = {item.name:item for item in self.plugins}



    def __getattr__(self, item):
        if item in self.name2plugin:
            return self.name2plugin[item].data
        else:
            super(AtomicData, self).__getattribute__(item)


