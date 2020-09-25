#ifndef TARDIS_CMONTECARLO_H
#define TARDIS_CMONTECARLO_H

#include <stdint.h>
#include "randomkit/randomkit.h"
#include "rpacket.h"
#include "status.h"
#include "storage.h"

#ifdef WITH_VPACKET_LOGGING
#define LOG_VPACKETS 1
#else
#define LOG_VPACKETS 0
#endif

typedef void (*montecarlo_event_handler_t) (rpacket_t *packet,
                                            storage_model_t *storage,
                                            double distance, rk_state *mt_state);

typedef double (*pt2sample_nu) (const rpacket_t *packet,
                                       const storage_model_t *storage,
                                       rk_state *mt_state);

void initialize_random_kit (unsigned long seed);

tardis_error_t line_search (const double *nu, double nu_insert,
                            int64_t number_of_lines, int64_t * result);

tardis_error_t
reverse_binary_search (const double *x, double x_insert,
                       int64_t imin, int64_t imax, int64_t * result);

tardis_error_t
binary_search (const double *x, double x_insert,
               int64_t imin, int64_t imax, int64_t * result);

double rpacket_doppler_factor (const rpacket_t *packet, const storage_model_t *storage);

double rpacket_inverse_doppler_factor (const rpacket_t *packet, const storage_model_t *storage);

void
angle_aberration_CMF_to_LF (rpacket_t * packet, const storage_model_t * storage);

double
angle_aberration_LF_to_CMF (rpacket_t *packet, const storage_model_t *storage, double mu);

/** Calculate the distance to shell boundary.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to shell boundary
 */
void compute_distance2boundary (rpacket_t * packet,
                                  const storage_model_t * storage);

/** Calculate the distance the packet has to travel until it redshifts to the first spectral line.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * @return distance to the next spectral line
 */
tardis_error_t compute_distance2line (rpacket_t * packet,
                                      const storage_model_t * storage);

/** Calculate the distance to the next continuum event, which can be a Thomson scattering, bound-free absorption or
  free-free transition.
 *
 * @param packet rpacket structure with packet information
 * @param storage storage model data
 *
 * sets distance to the next continuum event (in centimeters) in packet rpacket structure
 */
void compute_distance2continuum (rpacket_t *packet, storage_model_t *storage);

void macro_atom (rpacket_t * packet, const storage_model_t * storage, rk_state *mt_state);

void move_packet (rpacket_t * packet, storage_model_t * storage,
                    double distance);

void increment_j_blue_estimator (const rpacket_t * packet,
                                 storage_model_t * storage,
                                 double d_line,
                                 int64_t j_blue_idx);

void increment_Edotlu_estimator (const rpacket_t * packet,
                                 storage_model_t * storage,
                                 double d_line,
                                 int64_t j_blue_idx);

double get_increment_j_blue_estimator_energy (const rpacket_t * packet,
                                              const storage_model_t * storage,
                                              double d_line);

void
increment_continuum_estimators (const rpacket_t * packet, storage_model_t * storage, double distance,
                                double comov_nu, double comov_energy);

int64_t montecarlo_one_packet (storage_model_t * storage, rpacket_t * packet,
                               int64_t virtual_mode, rk_state *mt_state);

int64_t montecarlo_one_packet_loop (storage_model_t * storage,
                                    rpacket_t * packet,
                                    int64_t virtual_packet, rk_state *mt_state);

void montecarlo_main_loop(storage_model_t * storage,
                          int64_t virtual_packet_flag,
                          int nthreads,
                          unsigned long seed);

montecarlo_event_handler_t get_event_handler (rpacket_t * packet, storage_model_t * storage,
                                              double *distance, rk_state *mt_state);

void test_for_close_line(rpacket_t * packet, const storage_model_t * storage);

/* New handlers for continuum implementation */

montecarlo_event_handler_t montecarlo_continuum_event_handler(rpacket_t * packet, storage_model_t * storage, rk_state *mt_state);

void montecarlo_free_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance, rk_state *mt_state);

void montecarlo_bound_free_scatter (rpacket_t * packet, storage_model_t * storage, double distance, rk_state *mt_state);

double
bf_cross_section(const storage_model_t * storage, int64_t continuum_id, double comov_nu);

void calculate_chi_bf(rpacket_t * packet, storage_model_t * storage);

void calculate_chi_ff(rpacket_t * packet, const storage_model_t * storage);

void
move_packet_across_shell_boundary (rpacket_t * packet,
                                   storage_model_t * storage, double distance, rk_state *mt_state);

void
montecarlo_thomson_scatter (rpacket_t * packet, storage_model_t * storage,
                            double distance, rk_state *mt_state);

void
montecarlo_line_scatter (rpacket_t * packet, storage_model_t * storage,
                         double distance, rk_state *mt_state);

void
line_emission(rpacket_t * packet, storage_model_t * storage,
              int64_t emission_line_id, rk_state *mt_state);

double sample_nu_free_free(const rpacket_t * packet, const storage_model_t * storage, rk_state *mt_state);

double sample_nu_free_bound(const rpacket_t * packet, const storage_model_t * storage, rk_state *mt_state);

void continuum_emission(rpacket_t * packet, storage_model_t * storage, rk_state *mt_state,
                        pt2sample_nu sample_nu_continuum, int64_t emission_type_id);

void
create_vpacket (storage_model_t * storage, rpacket_t * packet,
                rk_state *mt_state);

#endif // TARDIS_CMONTECARLO_H
