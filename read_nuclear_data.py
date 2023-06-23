## A code to read hdf file and print some basic information

import pandas as pd

nuclear_data = '/Users/anirbandutta/Downloads/decay_radiation.h5'

with pd.HDFStore(nuclear_data) as hdf:
    # Print the hdf keys
    print(hdf.keys())

nuclear_df = pd.read_hdf(nuclear_data, key='decay_radiation')
print (nuclear_df.head(5))
print (nuclear_df.shape)

hdf.close()

