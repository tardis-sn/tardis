import numpy as np

def calculate_residuals(wavelength, luminosity, ref_wavelength, ref_luminosity):
    
    if not np.array_equal(wavelength, ref_wavelength):
        return None, None, False
        
    with np.errstate(divide='ignore', invalid='ignore'):
        residuals = np.where(ref_luminosity != 0, 
                             (luminosity - ref_luminosity) / ref_luminosity, 
                             0)
    return wavelength, residuals, True