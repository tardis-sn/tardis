#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "randomkit/randomkit.h"
#include "rpacket.h"
#include "status.h"

#ifdef __clang__
#define INLINE extern inline
#else
#define INLINE inline
#endif

typedef void (*montecarlo_event_handler_t) (rpacket_t * packet,
					    storage_model_t * storage,
					    double distance);

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

/** Calculate the distance to shell boundary.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to shell boundary
 */
inline double compute_distance2boundary (rpacket_t * packet,
					 storage_model_t * storage);

/** Calculate the distance the packet has to travel until it redshifts to the first spectral line.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the next spectral line
 */
extern inline tardis_error_t compute_distance2line (rpacket_t * packet,
						    storage_model_t * storage,
						    double *result);

/** Calculate the distance to the Thomson scatter event.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the Thomson scatter event in centimeters
 */
inline double compute_distance2electron (rpacket_t * packet,
					 storage_model_t * storage);

/** Calculate the distance to the next continuum event, which can be a Thomson scattering, bound-free absorption or
    free-free transition.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * sets distance to the next continuum event (in centimeters) in packet rpacket structure
 */
inline void compute_distance2continuum (rpacket_t * packet, storage_model_t * storage);

inline int64_t macro_atom (rpacket_t * packet, storage_model_t * storage);

inline double move_packet (rpacket_t * packet, storage_model_t * storage,
			   double distance);

inline void increment_j_blue_estimator (rpacket_t * packet,
					storage_model_t * storage,
					double d_line, int64_t j_blue_idx);

int64_t montecarlo_one_packet (storage_model_t * storage, rpacket_t * packet,
			       int64_t virtual_mode);

int64_t montecarlo_one_packet_loop (storage_model_t * storage,
				    rpacket_t * packet,
				    int64_t virtual_packet);

void montecarlo_main_loop(storage_model_t * storage, 
			  int64_t virtual_packet_flag,
			  int nthreads, 
			  unsigned long seed);

/* New handlers for continuum implementation */

inline montecarlo_event_handler_t montecarlo_continuum_event_handler(rpacket_t * packet, storage_model_t * storage);

void montecarlo_free_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance);

void montecarlo_bound_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance);

/* Other new stuff */

inline void macro_atom_new (rpacket_t * packet, storage_model_t * storage, next_interaction2process * macro_atom_deactivation_type,
 int activation2level_or_cont);

void e_packet(rpacket_t * packet, storage_model_t * storage, e_packet_type etype);

inline void line_emission(rpacket_t * packet, storage_model_t * storage);

double sample_nu_free_bound(rpacket_t * packet, storage_model_t * storage, int64_t continuum_id);

inline void bf_emission(rpacket_t * packet, storage_model_t * storage);

inline void test_for_close_line(rpacket_t * packet, storage_model_t * storage);

#endif // TARDIS_CMONTECARLO_H
