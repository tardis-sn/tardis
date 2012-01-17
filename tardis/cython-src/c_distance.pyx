#distance calculations in cython
from libc.math cimport sqrt

def distance2outer(double r_outer, double r, double mu):
    d = sqrt(r_outer**2 + ((mu**2 - 1.) * r**2)) - (r * mu)
    return d