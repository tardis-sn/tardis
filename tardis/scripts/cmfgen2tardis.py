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


def convert_format(file_path):

    with open(file_path, 'r') as f:
        for line in f:
            items = line.split()
            n = len(items)

            if 'data points' in line:
                df = pd.DataFrame(columns=np.arange(int(items[n - 1])),
                                  index=pd.Index([],
                                                 name='element'),
                                  dtype=np.float64)

            if 'mass fraction\n' in line:
                    abundances = []
                    element_string = items[0]
                    atomic_no = get_atomic_number(element_string.capitalize())
                    element_symbol = atomic_dataset.atom_data.loc[atomic_no]['symbol']

                    #Its a Isotope
                    if n == 4:
                        element_symbol += items[1]

                    for line in f:
                        items = line.split()
                        if items:
                            abundances.extend(
                                np.array(items).astype(np.float64))
                        else:
                            break

                    df.loc[element_symbol] = abundances
        return df.transpose()


def parse_file(args):
    df = convert_format(args.input_path)
    filename = os.path.basename(args.input_path)
    save_name = '.'.join((os.path.splitext(filename)[0], 'csv'))
    df.to_csv(os.path.join(args.output_path, save_name), index=False, sep=' ')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', help='Path to a CMFGEN abundance file')
    parser.add_argument(
        'output_path', help='Path to store converted TARDIS abundance file')
    args = parser.parse_args()
    parse_file(args)
