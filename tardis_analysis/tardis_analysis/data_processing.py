import pandas as pd
import numpy as np

def load_h5_data(file_path, keys):
    data = {}
    print(f"\nInspecting HDF5 file: {file_path}")
    try:
        with pd.HDFStore(file_path, mode='r') as store:
            for key in keys:
                base_path = f'/simulation/spectrum_solver/{key}'
                wavelength_path = f'{base_path}/wavelength'
                luminosity_path = f'{base_path}/luminosity'
                
                print(f"Checking {base_path}")
                
                try:
                    wavelength = pd.read_hdf(file_path, wavelength_path).to_numpy()
                    luminosity = pd.read_hdf(file_path, luminosity_path).to_numpy()
                    
                    data[key] = {
                        'wavelength': wavelength,
                        'luminosity': luminosity
                    }
                    print(f"Loaded data for {key}: wavelength shape={wavelength.shape}, luminosity shape={luminosity.shape}")
                except (KeyError, ValueError) as e:
                    print(f"Warning: Missing required paths for {key}")
                    continue
                
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return data

def calculate_residuals(wavelength, luminosity, ref_wavelength, ref_luminosity):
    if not np.array_equal(wavelength, ref_wavelength):
        return None, None, False
        
    with np.errstate(divide='ignore', invalid='ignore'):
        residuals = np.where(ref_luminosity != 0, 
                           (luminosity - ref_luminosity) / ref_luminosity, 
                           0)
    return wavelength, residuals, True 