#reading different model files
from numpy import recfromtxt, genfromtxt
from astropy import units as u

def read_tardis_simple_ascii_density(fname):
    """
    Reading a density file of the following structure (example; lines starting with a hash will be ignored):
    5 s
    #index velocity [km/s] density [g/cm^3]
    0 53e5 1.6e-

    Parameters
    ----------

    fname: str
        filename or path with filename


    Returns
    -------


    """

    with open(fname) as fh:
        time_of_model = u.Quantity(fh.readline().strip())

    data = recfromtxt(fname, skip_header=1, names=('index', 'velocity', 'density'), dtype=(int, float, float))




def read_tardis_simple_ascii_abundances(fname):
    pass
