import os
import sys
import argparse
import numpy as np
import pandas as pd

from tardis.atomic import AtomData


atomic_dataset = AtomData.from_hdf5()


def get_atomic_number(element):
    index = -1
    for atomic_no, row in atomic_dataset.atom_data.iterrows():
        if element in row['name']:
            index = atomic_no
            break
    return index


def extract_file_block(f):
    qty = []

    for line in f:
        items = line.split()
        if items:
            qty.extend(np.array(items).astype(np.float64))
        else:
            break

    qty = np.array(qty)

    return qty[::-1]


def convert_format(file_path):

    with open(file_path, 'r') as f:
        for line in f:
            items = line.split()
            n = len(items)

            if 'data points' in line:
                abundances_df = pd.DataFrame(columns=np.arange(int(items[n - 1])),
                                             index=pd.Index([],
                                                            name='element'),
                                             dtype=np.float64)
            if 'Time' in line:
                time_of_model = float(items[n - 1])
            if 'Velocity' in line:
                velocity = extract_file_block(f)
            if 'Density' in line:
                density = extract_file_block(f)
            if 'Electron density' in line:
                electron_density = extract_file_block(f)
            if 'Temperature' in line:
                temperature = extract_file_block(f)

            if 'mass fraction\n' in line:
                    element_string = items[0]
                    atomic_no = get_atomic_number(element_string.capitalize())
                    element_symbol = atomic_dataset.atom_data.loc[atomic_no]['symbol']

                    #Its a Isotope
                    if n == 4:
                        element_symbol += items[1]

                    abundances = extract_file_block(f)
                    abundances_df.loc[element_symbol] = abundances

        density_df = pd.DataFrame.from_records(
            [velocity, density, electron_density, temperature]).transpose()
        density_df.columns = ['Velocity', 'Densities',
                              'ElectronDensities', 'Temperature']

        return abundances_df.transpose(), density_df, time_of_model


def parse_file(args):
    abundances_df, density_df, time_of_model = convert_format(args.input_path)

    filename = os.path.splitext(os.path.basename(args.input_path))[0]
    abundances_fname = '.'.join((filename, 'csv'))
    density_fname = '.'.join((filename, 'dat'))

    abundances_file_path = os.path.join(
        args.output_path, abundances_fname)
    density_file_path = os.path.join(
        args.output_path, density_fname)

    abundances_df.to_csv(abundances_file_path, index=False, sep=' ')

    with open(density_file_path, 'w') as f:
        f.write(" ".join((str(time_of_model), "day")))
        f.write("\n")
    density_df.to_csv(density_file_path, index=True,
                      index_label='Index', sep=' ', mode='a')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', help='Path to a CMFGEN file')
    parser.add_argument(
        'output_path', help='Path to store converted TARDIS format files')
    args = parser.parse_args()
    parse_file(args)
