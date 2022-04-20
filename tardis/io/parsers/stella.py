import re
import pandas as pd
from astropy import units as u
import numpy as np
import yaml
from scipy import interpolate
import os
import pathlib

STELLA_CONFIG_TEMPLATE_FNAME = os.path.join(
    os.path.dirname(__file__), "tardis_stella_template.yml.j2"
)

STELLA_METADATA_MAP = {
    "days post max Lbol" : "time_explosion",
    "inner boundary mass" : "mass_inner",
    "total mass" : "mass_total",
}

STELLA_COL_MAPPER = {
    "mass of cell (g)": "mass_cell",
    "cell center m (g)": "mass_total",
    "cell center R (cm)": "r_center",
    "cell center v (cm/s)": "v_center",
    "avg density": "density",
    "radiation pressure": "pressure_rad",
    "avg temperature": "t_gas",
    "radiation temperature": "t_radiation",
    "avg opacity": "opacity",
    "outer edge m (g)": "mass_outer",
    "outer edge r (cm)": "r_outer",
}

def read_stella(fname):
    """
    Read a STELLA model into a pandas dataframe
    """
    with open(fname) as fh:
        metadata_raw = dict(
            [re.split("\s{3,}", fh.readline().strip())[:2] 
             for _ in range(4)]
        )
        for key in metadata_raw:
            metadata_raw[key] = [float(metadata_raw[key])]
        metadata_raw['zones'] = [int(metadata_raw['zones'][0])]
        metadata = pd.DataFrame(metadata_raw).rename(
            columns=STELLA_METADATA_MAP
        )
        for line in fh:
            if line.strip().startswith("mass"):
                break
    raw_columns = re.split("\s{3,}", line.strip())

    stella_model = pd.read_csv(
        fname, delim_whitespace=True, skiprows=6, names=raw_columns, index_col=0
    )
    return metadata, stella_model.rename(columns=STELLA_COL_MAPPER)

class STELLAManager:
    '''
    Manages a directory of files used to generate
    spectra using stella
    
    Parameters
    ----------
        stella_dir : str or pathlib.Path object
            Directory containing the Stella Data
    '''
    
    def __init__(self, stella_dir):
        
        self.stella_dir = pathlib.Path(stella_dir)
        self._data = None
        self._metadata = None
        self._models = None
        
    
    def get_data_files(self):
        '''
        Get all files ending in .data
        describing each epoch
        '''
        
        return self.stella_dir.glob('*.data')
    
    def load_data_files(self):
        '''
        Load all models and metadata from each data model in
        the directory
        '''
        
        self._data = {filename: read_stella(filename)
                    for filename in self.get_data_files()}
        
    @property
    def metadata(self):
        """This should be a dataframe
        multiindexed by filename then the normal 
        metadata columns
        """
        
        if self._metadata is None:
            metadata, models = zip(*self.data.values())
            metadata_df = pd.concat(metadata)
            filenames = self.data.keys()
            metadata_df['filename'] = filenames
            metadata_df.set_index('filename', inplace=True)
            self._metadata = metadata_df

        return self._metadata
    
    @property
    def models(self):
        """This should be a DataFrame
        multiindexed by filename then the normal 
        model columns"""
        if self._models is None:
            metadata, models = zip(*self.data.values())
            filenames = self.data.keys()
            self._models = pd.concat(models, 
                                     keys=filenames, 
                                     names=['filename', 'zone'])

        return self._models
    
    @property
    def data(self):
        
        if self._data is None:
            self.load_data_files()
        return self._data
   
    def query_metadata(self, query):
        '''shortcut for self.metadata.query'''
        
        return self.metadata.query(query)
    
    def query_models(self, query):
        '''shortcut for self.metadata.query'''
        
        return self.models.query(query)


def convert_stella_to_csvy(fname, out_fname=None):
    """
    Convert STELLA to TARDIS file
    """
    stella_meta, stella_model = read_stella(fname)
    v_interpolator = interpolate.interp1d(
        stella_model["r_center"], stella_model["v_center"]
    )
    # velocity needed for future check against homology
    #velocity = v_interpolator(stella_model["r_outer"].values) * u.cm / u.s
    time_explosion = stella_meta.time_explosion.values * u.day
    homologous_velocity = (
        stella_model.r_outer.values * u.cm / (time_explosion)
    )
    csvy_meta = {
        "name": "STELLA transformation",
        "model_density_time_0": str(time_explosion),
        "model_isotope_time_0": str(time_explosion),
        "description": "parsed stella file",
        "tardis_model_config_version": "v1.0",
        "datatype": {"fields": []},
    }
    csvy_table = pd.DataFrame()
    csvy_table["density"] = stella_model.density
    csvy_meta["datatype"]["fields"].append(
        {
            "name": "density",
            "unit": "g/cm^3",
            "desc": "Average density from stella",
        }
    )

    csvy_table["velocity"] = homologous_velocity.to(u.km / u.s).value
    csvy_meta["datatype"]["fields"].append(
        {
            "name": "velocity",
            "unit": "km/s",
            "desc": "WARNING - THIS IS NOT THE STELLA VELOCITY but a homologous approximation given radius and t_explosion",
        }
    )

    abundances = stella_model.iloc[:, 12:-3]
    for isotope in abundances.columns:
        csvy_meta["datatype"]["fields"].append(
            {"name": isotope, "desc": "{0} mass fraction".format(isotope)}
        )

    csvy_table = csvy_table.join(abundances)
    if out_fname is None:
        out_fname = f"{fname}_stella2tardis.csvy"

    with open(out_fname, "w") as fh:
        fh.write('---\n')
        yaml.dump(csvy_meta, fh)
        fh.write('---\n')
        csvy_table.to_csv(fh, index=False)
