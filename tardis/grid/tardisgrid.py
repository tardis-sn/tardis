import pandas as pd
import numpy as np
import copy
import tardis

from tardis.io.config_reader import Configuration

def set_tardis_config_property(tardis_config, key, value):
    keyitems = key.split('.')
    tmp_dict = getattr(tardis_config, keyitems[0])
    for key in keyitems[1:-1]:
        tmp_dict = getattr(tmp_dict, key)
    setattr(tmp_dict, keyitems[-1], value)
    return
        

class tardisGrid:
    
    def __init__(self, configFile, gridFrame):
        try:
            tardis_config = Configuration.from_yaml(configFile)
        except TypeError:
            tardis_config = Configuration.from_config_dict(configFile)
            
        self.config = tardis_config
        self.grid = gridFrame
        
        return
    
    def grid_row_to_config(self, row_index):
        tmp_config = copy.deepcopy(self.config)
        grid_row = self.grid.iloc[row_index]
        for colname, value in zip(self.grid.columns, grid_row.values):
            set_tardis_config_property(tmp_config, colname, value)
        return tmp_config
        
    
    def run_sim_from_grid(self, row_index):
        tardis_config = self.grid_row_to_config(row_index)
        tardis_config.pop('config_dirname')
        sim = tardis.run_tardis(tardis_config)
        return sim
    
    def save_grid(self, filename):
        self.grid.to_csv(filename, index=False)
        return
    
    @classmethod
    def from_axes(cls, configFile, axesdict):
        axes = []
        dim = 1
        for key in axesdict:
            ax = axesdict[key]
            axes.append(ax)
            dim = dim * len(ax)
        axesmesh = np.meshgrid(*axes)
        tmp = np.dstack(axesmesh)
        gridpoints = tmp.reshape((dim, len(axes)), order='F')
        df = pd.DataFrame(data=gridpoints, columns=axesdict.keys())
        return cls(configFile=configFile, gridFrame=df)
    
    
