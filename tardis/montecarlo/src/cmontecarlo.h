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
#include "cmontecarlo1.h"

typedef void (*montecarlo_event_handler_t) (rpacket_t * packet,
					    storage_model_t * storage,
					    double distance);

double rpacket_doppler_factor(const rpacket_t *packet, const storage_model_t *storage);

/** Calculate the distance to shell boundary.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to shell boundary
 */
double compute_distance2boundary (rpacket_t * packet,
					 const storage_model_t * storage);

/** Calculate the distance the packet has to travel until it redshifts to the first spectral line.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the next spectral line
 */
tardis_error_t compute_distance2line (const rpacket_t * packet,
						    const storage_model_t * storage,
						    double *result);

/** Calculate the distance to the next continuum event, which can be a Thomson scattering, bound-free absorption or
    free-free transition.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * sets distance to the next continuum event (in centimeters) in packet rpacket structure
 */
void compute_distance2continuum (rpacket_t * packet, storage_model_t * storage);

int64_t macro_atom (const rpacket_t * packet, const storage_model_t * storage);

double move_packet (rpacket_t * packet, storage_model_t * storage,
			   double distance);

void increment_j_blue_estimator (const rpacket_t * packet,
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

montecarlo_event_handler_t montecarlo_continuum_event_handler(rpacket_t * packet, storage_model_t * storage);

void montecarlo_free_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance);

void montecarlo_bound_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance);

double
bf_cross_section(const storage_model_t * storage, int64_t continuum_id, double comov_nu);

void calculate_chi_bf(rpacket_t * packet, storage_model_t * storage);

void
move_packet_across_shell_boundary (rpacket_t * packet,
                                   storage_model_t * storage, double distance);

void
montecarlo_thomson_scatter (rpacket_t * packet, storage_model_t * storage,
                            double distance);

void
montecarlo_line_scatter (rpacket_t * packet, storage_model_t * storage,
                         double distance);

#endif // TARDIS_CMONTECARLO_H
