#ifndef TARDIS_RPACKET_H
#define TARDIS_RPACKET_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "randomkit/randomkit.h"
#include "status.h"
#include "storage.h"

#ifdef __clang__
#define INLINE extern inline
#else
#define INLINE inline
#endif

#define MISS_DISTANCE 1e99
#define C 29979245800.0
#define INVERSE_C 3.33564095198152e-11
#define H 6.6260755e-27		// erg * s, converted to CGS units from the NIST Constant Index
#define KB 1.3806488e-16	// erg / K converted to CGS units from the NIST Constant Index

/**
 * @brief A photon packet.
 */
typedef struct RPacket
{
  double nu; /**< Frequency of the packet in Hz. */
  double mu; /**< Cosine of the angle of the packet. */
  double energy; /**< Energy of the packet in erg. */
  double r; /**< Distance from center in cm. */
  double tau_event;
  double nu_line;
  int64_t current_shell_id; /**< ID of the current shell. */
  int64_t next_line_id;	/**< The index of the next line that the packet will encounter. */
  /**
   * @brief The packet has a nu red-ward of the last line.
   * It will not encounter any lines anymore.
   */
  int64_t last_line;
  /** 
   * @brief The packet just encountered a line that is very close to the next line.
   * The next iteration will automatically make an interaction with the next line 
   * (avoiding numerical problems).
   */
  int64_t close_line;
  /** 
   * @brief The packet has recently crossed the boundary and is now sitting on the boundary.
   * To avoid numerical errors, make sure that d_inner is not calculated. The value is -1
   * if the packed moved inwards, 1 if the packet moved outwards and 0 otherwise.
   */
  int64_t recently_crossed_boundary;
  /**
   * @brief packet is a virtual packet and will ignore any d_line or d_electron checks.
   * It now whenever a d_line is calculated only adds the tau_line to an
   * internal float.
   */
  int64_t current_continuum_id; /* Packet can interact with bf-continua with an index equal or bigger than this */
  int64_t virtual_packet_flag;
  int64_t virtual_packet;
  double d_line; /**< Distance to line event. */
  double d_electron; /**< Distance to electron scattering event. */
  double d_boundary; /**< Distance to shell boundary. */
  double d_cont; /**< Distance to continuum event */
  int64_t next_shell_id; /**< ID of the next shell packet visits. */
  rpacket_status_t status; /**< Packet status (in process, emitted or reabsorbed). */
  int64_t id;
  double chi_th; /**< Opacity due to electron scattering */
  double chi_cont; /**< Opacity due to continuum processes */
  double chi_ff; /**< Opacity due to free-free processes */
  double chi_bf; /**< Opacity due to bound-free processes */
} rpacket_t;

inline double rpacket_get_nu (rpacket_t * packet);

inline void rpacket_set_nu (rpacket_t * packet, double nu);

inline double rpacket_get_mu (rpacket_t * packet);

inline void rpacket_set_mu (rpacket_t * packet, double mu);

inline double rpacket_get_energy (rpacket_t * packet);

inline void rpacket_set_energy (rpacket_t * packet, double energy);

inline double rpacket_get_r (rpacket_t * packet);

inline void rpacket_set_r (rpacket_t * packet, double r);

inline double rpacket_get_tau_event (rpacket_t * packet);

inline void rpacket_set_tau_event (rpacket_t * packet, double tau_event);

inline double rpacket_get_nu_line (rpacket_t * packet);

inline void rpacket_set_nu_line (rpacket_t * packet, double nu_line);

inline unsigned int rpacket_get_current_shell_id (rpacket_t * packet);

inline void rpacket_set_current_shell_id (rpacket_t * packet,
            unsigned int current_shell_id);

inline unsigned int rpacket_get_next_line_id (rpacket_t * packet);

inline void rpacket_set_next_line_id (rpacket_t * packet,
              unsigned int next_line_id);

inline bool rpacket_get_last_line (rpacket_t * packet);

inline void rpacket_set_last_line (rpacket_t * packet, bool last_line);

inline bool rpacket_get_close_line (rpacket_t * packet);

inline void rpacket_set_close_line (rpacket_t * packet, bool close_line);

inline int rpacket_get_recently_crossed_boundary (rpacket_t * packet);

inline void rpacket_set_recently_crossed_boundary (rpacket_t * packet,
               int
               recently_crossed_boundary);

inline int rpacket_get_virtual_packet_flag (rpacket_t * packet);

inline void rpacket_set_virtual_packet_flag (rpacket_t * packet,
               int virtual_packet_flag);

inline int rpacket_get_virtual_packet (rpacket_t * packet);

inline void rpacket_set_virtual_packet (rpacket_t * packet,
          int virtual_packet);

inline double rpacket_get_d_boundary (rpacket_t * packet);

inline void rpacket_set_d_boundary (rpacket_t * packet, double d_boundary);

inline double rpacket_get_d_electron (rpacket_t * packet);

inline void rpacket_set_d_electron (rpacket_t * packet, double d_electron);

inline double rpacket_get_d_line (rpacket_t * packet);

inline void rpacket_set_d_line (rpacket_t * packet, double d_line);

inline int rpacket_get_next_shell_id (rpacket_t * packet);

inline void rpacket_set_next_shell_id (rpacket_t * packet, int next_shell_id);

inline rpacket_status_t rpacket_get_status (rpacket_t * packet);

inline void rpacket_set_status (rpacket_t * packet, rpacket_status_t status);

inline int rpacket_get_id (rpacket_t * packet);

inline void rpacket_set_id (rpacket_t * packet, int id);

inline void rpacket_reset_tau_event (rpacket_t * packet);

tardis_error_t rpacket_init (rpacket_t * packet, storage_model_t * storage,
           int packet_index, int virtual_packet_flag);

void initialize_random_kit (unsigned long seed);

/* New getter and setter methods for continuum implementation */

inline void rpacket_set_d_continuum (rpacket_t * packet, double d_continuum);

inline double rpacket_get_d_continuum (rpacket_t * packet);

inline void rpacket_set_chi_electron (rpacket_t * packet, double chi_electron);

inline double rpacket_get_chi_electron (rpacket_t * packet);

inline void rpacket_set_chi_continuum (rpacket_t * packet, double chi_continuum);

inline double rpacket_get_chi_continuum (rpacket_t * packet);

inline void rpacket_set_chi_freefree (rpacket_t * packet, double chi_freefree);

inline double rpacket_get_chi_freefree (rpacket_t * packet);

inline void rpacket_set_chi_boundfree (rpacket_t * packet, double chi_boundfree);

inline double rpacket_get_chi_boundfree (rpacket_t * packet);

inline unsigned int rpacket_get_current_continuum_id (rpacket_t * packet);

inline void rpacket_set_current_continuum_id (rpacket_t * packet, unsigned int current_continuum_id);

#endif // TARDIS_RPACKET_H
