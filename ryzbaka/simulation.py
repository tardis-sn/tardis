'''
this is a script for running the example tardis simulation
'''
import time
import random
import pickle
from tardis import run_tardis
from tardis.io.atom_data.util import download_atom_data
from matplotlib.pyplot import plot,show,xlabel,ylabel,legend,show,title,xlim
import os

def get_config():
    #download example config file if not in directory
    os.system('wget https://raw.githubusercontent.com/tardis-sn/tardis/master/docs/models/examples/tardis_example.yml')

def plot_spectrum(sim):
    spectrum = sim.runner.spectrum
    spectrum_virtual = sim.runner.spectrum_virtual
    spectrum_integrated = sim.runner.spectrum_integrated

    title('Spectrum Plot')
    plot(spectrum.wavelength, spectrum.luminosity_density_lambda, label='normal packets')
    plot(spectrum.wavelength, spectrum_virtual.luminosity_density_lambda, label='virtual packets')
    plot(spectrum.wavelength, spectrum_integrated.luminosity_density_lambda, label='formal integral')
    xlabel('Wavelength [$\AA$]')
    ylabel('Luminosity [erg/s/$\AA$]')
    legend()
    xlim(3000, 9000)
    show()

def save(simulation):
    '''
    Saves the simulation object as a pickle file that can be loaded later for use.
    '''
    filename='Simulation_'+str(time.time())+str(random.randint(1,1000))+'.pickle'
    with open(filename,'wb') as f:
        pickle.dump(simulation,f)

def main():
    files=os.listdir()
    if not 'tardis_example.yml' in files:
        get_config()

    download_atom_data('kurucz_cd23_chianti_H_He')
    sim=run_tardis('tardis_example.yml')
    save(sim)
if __name__=='__main__':
    main()
