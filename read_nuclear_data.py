## A code to read hdf file and print some basic information

import pandas as pd

nuclear_data = '/Users/anirbandutta/Downloads/compiled_ensdf_csv.h5'

with pd.HDFStore(nuclear_data) as hdf:
    # Print the hdf keys
    print(hdf.keys())

nuclear_df = pd.read_hdf(nuclear_data, key='decay_data')
print (nuclear_df.head(5))
print (nuclear_df.shape)

hdf.close()

