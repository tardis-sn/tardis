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
        if element in row["name"]:
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
    # Convention in CMFGEN files is different from TARDIS files.
    # CMFGEN stores velocity values in decreasing order, while
    # TARDIS stores it in ascending order, so alongwith
    # velocity, all columns will be reversed.
    return qty[::-1]


def convert_format(file_path):
    quantities_row = []
    prop_list = ["Velocity", "Density", "Electron density", "Temperature"]
    with open(file_path, "r") as f:
        for line in f:
            items = line.replace("(", "").replace(")", "").split()
            n = len(items)

            if "data points" in line:
                abundances_df = pd.DataFrame(
                    columns=np.arange(int(items[n - 1])),
                    index=pd.Index([], name="element"),
                    dtype=np.float64,
                )
            if any(prop in line for prop in prop_list):
                quantities_row.append(items[n - 1].replace("gm", "g"))
            if "Time" in line:
                time_of_model = float(items[n - 1])
            if "Velocity" in line:
                velocity = extract_file_block(f)
            if "Density" in line:
                density = extract_file_block(f)
            if "Electron density" in line:
                electron_density = extract_file_block(f)
            if "Temperature" in line:
                temperature = extract_file_block(f)

            if "mass fraction\n" in line:
                element_string = items[0]
                atomic_no = get_atomic_number(element_string.capitalize())
                element_symbol = atomic_dataset.atom_data.loc[atomic_no][
                    "symbol"
                ]

                # Its a Isotope
                if n == 4:
                    element_symbol += items[1]

                abundances = extract_file_block(f)
                abundances_df.loc[element_symbol] = abundances

        density_df = pd.DataFrame.from_records(
            [velocity, temperature * 10 ** 4, density, electron_density]
        ).transpose()
        density_df.columns = [
            "velocity",
            "temperature",
            "densities",
            "electron_densities",
        ]
        quantities_row += abundances_df.shape[0] * [1]
        return (
            abundances_df.transpose(),
            density_df,
            time_of_model,
            quantities_row,
        )


def parse_file(args):
    abundances_df, density_df, time_of_model, quantities_row = convert_format(
        args.input_path
    )

    filename = os.path.splitext(os.path.basename(args.input_path))[0]
    save_fname = ".".join((filename, "csv"))
    resultant_df = pd.concat([density_df, abundances_df], axis=1)
    resultant_df.columns = pd.MultiIndex.from_tuples(
        zip(resultant_df.columns, quantities_row)
    )
    save_file_path = os.path.join(args.output_path, save_fname)
    with open(save_file_path, "w") as f:
        f.write(" ".join(("t0:", str(time_of_model), "day")))
        f.write("\n")

    resultant_df.to_csv(save_file_path, index=False, sep=" ", mode="a")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Path to a CMFGEN file")
    parser.add_argument(
        "output_path", help="Path to store converted TARDIS format files"
    )
    args = parser.parse_args()
    parse_file(args)
