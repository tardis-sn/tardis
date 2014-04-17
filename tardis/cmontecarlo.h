#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

/** Look for a place to insert a value in an inversely sorted float array.
 *
 * @param x an inversely (largest to lowest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
int binary_search(double *x, double x_insert, int imin, int imax);

#endif // TARDIS_CMONTECARLO_H
