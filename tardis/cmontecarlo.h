#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

#include <Python.h>
#include <numpy/arrayobject.h>

/** Look for a place to insert a value in an inversely sorted float array.
 *
 * @param x an inversely (largest to lowest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
npy_int64 binary_search(npy_float64 *x, npy_float64 x_insert, npy_int64 imin, npy_int64 imax);

#endif // TARDIS_CMONTECARLO_H
