#Tiny deduplication script to fix duplicate issues with the aoife_atomic_data file

import pandas as pd

src = "aoife_atomic_data.h5"
dst = "aoife_atomic_data_levels_dedup.h5"

with pd.HDFStore(src, "r") as inp, pd.HDFStore(dst, "w") as out:
    for key in inp.keys():
        df = inp[key]

        if key == "/levels":
            tmp = df.reset_index()
            n_before = len(tmp)
            tmp = tmp.loc[~tmp.duplicated()].copy()
            n_after = len(tmp)
            print(f"/levels: removed {n_before - n_after} exact duplicate rows")
            df = tmp.set_index(["atomic_number", "ion_number", "level_number"]).sort_index()

        out.put(key.lstrip("/"), df)

print(f"Wrote {dst}")