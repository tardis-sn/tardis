def parse_atomic_data(filename):
    """
    A parser to extract data from si2_osc_kurucz to an hdf file.

    Parameters
    ----------
    filename : string
        string defining path of the atomic data file
    """
    import re
    import pandas as pd
    import numpy as np
    import os
    from pandas import HDFStore

    with open(filename) as file:

        COLUMNS_ENERGY_LEVELS = ['level', 'g', 'E(cm^-1)', '10^15 Hz', 'eV', 'Lam(A)', 'ID', 'ARAD', 'C4', 'C6']
        COLUMNS_OSCILLATOR_STRENGTHS = ['From', 'To', 'f', 'A', 'Lam(A)', 'i','j', 'Lam(obs)', '% Acc']

        meta_data = {}

        META_PATTERN = re.compile(r'^(\d{1,2}-[a-zA-Z]+-\d{4}|\d+\s|\d+\.\d+)\s*!([a-zA-Z\s]+[a-zA-Z]$)')
        PATTERN1 = re.compile(r"\d[a-z][^\s-]+")
        PATTERN2 = re.compile(r"[+-]?\d+\d+[eE][-+]?\d+|\d+\.\d+|\d+-\s+\d+|\s[+-]?\d+\s")

        energy_levels = []
        oscillator_strengths = []

        flag = 0

        line = file.readline()
        while(line):
            if(bool(re.match("\*+",line))):
                flag = flag+1
                line = file.readline()
                while(bool(re.match("\*+",line)) is False):
                    line = file.readline()
                flag = flag+1

            if(flag<=2):
                while(bool(META_PATTERN.match(line))):
                    pair = META_PATTERN.search(line)
                    meta_data[pair.group(2)] = pair.group(1)
                    line = file.readline()
                while(PATTERN1.match(line)):
                    ls = PATTERN1.findall(line)
                    ls.extend(PATTERN2.findall(line))
                    energy_levels.append(ls)
                    line = file.readline()
            else:
                while(PATTERN1.match(line)):
                    ls = PATTERN1.findall(line)
                    ls.extend(PATTERN2.findall(line))
                    oscillator_strengths.append(ls)
                    line = file.readline()

            line = file.readline()

    for tr in oscillator_strengths:
        temp = tr[-1]
        tr.pop()
        ls = temp.split('-')
        tr.append(int(ls[0].strip()))
        tr.append(int(ls[1].strip()))
        for i in range(2,len(tr)-2):
            tr[i] = float(tr[i])
        tr.extend([np.nan,np.nan])

    for e_level in energy_levels:
        for i in range(1,len(e_level)):
            e_level[i] = float(e_level[i].strip())
        e_level[6] = int(e_level[6])

    df1 = pd.DataFrame(energy_levels,columns = COLUMNS_ENERGY_LEVELS)
    df1.index.name = "energyLevels"
    df2 = pd.DataFrame(oscillator_strengths,columns = COLUMNS_OSCILLATOR_STRENGTHS)
    df2.index.name = "oscillator_strengths"

    df1.to_hdf('atomic_data.h5',df1.index.name)
    df2.to_hdf('atomic_data.h5',df2.index.name)

import pandas as pd

parse_atomic_data("si2_osc_kurucz")
df1 = pd.read_hdf("atomic_data.h5",key="energyLevels")
df2 = pd.read_hdf("atomic_data.h5",key="oscillator_strengths")
print(df1.head())
print(df2.head())
