#ifndef TARDIS_CMONTECARLO1_H
#define TARDIS_CMONTECARLO1_H

/** Insert a value in to an array of line frequencies
 *
 * @param nu array of line frequencies
 * @param nu_insert value of nu key
 * @param number_of_lines number of lines in the line list
 *
 * @return index of the next line ot the red. If the key value is redder than the reddest line returns number_of_lines.
 */
tardis_error_t line_search (const double *nu, double nu_insert,
                            int64_t number_of_lines, int64_t * result);

#endif
