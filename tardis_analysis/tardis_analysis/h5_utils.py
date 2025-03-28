# tardis_analysis/h5_utils.py

import h5py

def load_h5_data(file_path, keys):
    data = {}
    print(f"\nInspecting HDF5 file: {file_path}")
    try:
        with h5py.File(file_path, 'r') as f:
            print("Top-level keys:", list(f.keys()))
            if 'simulation' not in f:
                print("Error: 'simulation' group not found in the HDF5 file.")
                return data
            if 'spectrum_solver' not in f['simulation']:
                print("Error: 'spectrum_solver' group not found under 'simulation'.")
                return data
            
            for key in keys:
                group_path = f'simulation/spectrum_solver/{key}'
                print(f"Checking {group_path}")
                if group_path in f:
                    wavelength_path = f'{group_path}/wavelength/values'
                    luminosity_path = f'{group_path}/luminosity/values'
                    
                    if wavelength_path in f and luminosity_path in f:
                        wavelength = f[wavelength_path][()]
                        luminosity = f[luminosity_path][()]
                        data[key] = {
                            'wavelength': wavelength,
                            'luminosity': luminosity
                        }
                        print(f"Loaded data for {key}: wavelength shape={wavelength.shape}, luminosity shape={luminosity.shape}")
                    else:
                        print(f"Warning: Missing required paths for {key}")
                else:
                    print(f"Warning: {group_path} not found")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return data



