#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "randomkit/randomkit.h"

#include "errors.h"


#ifdef __clang__
#define INLINE extern inline
#else
#define INLINE inline
#endif

#define MISS_DISTANCE 1e99
#define C 29979245800.0
#define INVERSE_C 3.33564095198152e-11


/** Look for a place to insert a value in an inversely sorted float array.
 *
 * @param x an inversely (largest to lowest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
inline tardis_error_t reverse_binary_search (double *x, double x_insert,
					     int64_t imin, int64_t imax,
					     int64_t * result);

/** Look for a place to insert a value in a sorted float array.
 *
 * @param x a (lowest to largest) sorted float array
 * @param x_insert a value to insert
 * @param imin lower bound
 * @param imax upper bound
 *
 * @return index of the next boundary to the left
 */
inline tardis_error_t binary_search (double *x, double x_insert, int64_t imin,
				     int64_t imax, int64_t * result);

/** Insert a value in to an array of line frequencies
 *
 * @param nu array of line frequencies
 * @param nu_insert value of nu key
 * @param number_of_lines number of lines in the line list
 *
 * @return index of the next line ot the red. If the key value is redder than the reddest line returns number_of_lines.
 */
inline tardis_error_t line_search (double *nu, double nu_insert,
				   int64_t number_of_lines, int64_t * result);

