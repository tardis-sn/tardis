## Contents

`Universal_SDec_LIV_Plotter.py` - Takes any config + atom data combination (--config + --atomic_data) and plots the spectrum, SDec plot and LIV plot. 

`atomic_cleanup.py` - A bit of a hotfix for the aoife_atomic_data.h5 file (used by helium_HM_12_recomb.yml). Takes an atomic data file and deduplicates the levels so that they don't feed too many rows to tardis, which otherwise fails to load data

/atomic-data - atom data .h5 files used by the tested configs. Not included in PR because of size.

/configs - tardis YAML configs tested
  /dependencies - dependencies for the above YAML configs

/plots - default outdir for the spectra, SDec and LIV plots, mostly to avoid cluttering main folder. 

## Usage

`Universal_SDec_LIV_Plotter.py`

Currently mostly proof of concept, runs with CLI from terminal. 

`--config` - Input TARDIS config file for the simulation. Accepts PATH and Configuration object.   

`--atomic_data` - Input atomic data for the simulation. Accepts PATH, filename and AtomData object.   

`--plot_pdf` - Simple bool to decide whether or not to plot pdf in addition to pngs. Default is True. Accepts 0 or 1, True or False.   

`--outdir` - Specify the out directory for the plots. Default is the current folder. Accepts path.

### Examples

Example 1 - uses defaults (outdir is /plots and it plots pdfs)
```{bash}
python Universal_SDec_LIV_plotter.py --atomic_data /atomic-data/kurucz_cd23_chianti_H_He_latest.h5 --config /configs/SNIax.yml
```

Example 2 - No pdfs, and with outdir specified

```{bash}
python Universal_SDec_LIV_plotter.py --atomic_data /atomic-data/kurucz_cd23_chianti_H_He_latest.h5 --config /configs/SNIax.yml --plot_pdf False --outdir ~/plot_gallery
```