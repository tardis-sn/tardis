#ifndef TARDIS_RPACKET_H
#define TARDIS_RPACKET_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include "randomkit/randomkit.h"
#include "status.h"
#include "storage.h"

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
   * @brief packet is a virtual packet and will ignore any d_line or d_electron checks.
   * It now whenever a d_line is calculated only adds the tau_line to an
   * internal float.
   */
  int64_t current_continuum_id; /* Packet can interact with bf-continua with an index equal or bigger than this */
  int64_t virtual_packet_flag;
  int64_t virtual_packet;
  double d_line; /**< Distance to electron event. */
  double d_electron; /**< Distance to line event. */
  double d_boundary; /**< Distance to shell boundary. */
  double d_cont; /**< Distance to continuum event */
  int64_t next_shell_id; /**< ID of the next shell packet visits. */
  rpacket_status_t status; /**< Packet status (in process, emitted or reabsorbed). */
  int64_t id;
  double chi_th; /**< Opacity due to electron scattering */
  double chi_cont; /**< Opacity due to continuum processes */
  double chi_ff; /**< Opacity due to free-free processes */
  double chi_bf; /**< Opacity due to bound-free processes */
  double *chi_bf_tmp_partial;
  int64_t macro_atom_activation_level;
  bool compute_chi_bf;
  double vpacket_weight;
} rpacket_t;

static inline double rpacket_get_nu (const rpacket_t * packet)
{
  return packet->nu;
}

static inline void rpacket_set_nu (rpacket_t * packet, double nu)
{
  packet->nu = nu;
}

static inline double rpacket_get_mu (const rpacket_t * packet)
{
  return packet->mu;
}

static inline void rpacket_set_mu (rpacket_t * packet, double mu)
{
  packet->mu = mu;
}

static inline double rpacket_get_energy (const rpacket_t * packet)
{
  return packet->energy;
}

static inline void rpacket_set_energy (rpacket_t * packet, double energy)
{
  packet->energy = energy;
}

static inline double rpacket_get_r (const rpacket_t * packet)
{
  return packet->r;
}

static inline void rpacket_set_r (rpacket_t * packet, double r)
{
  packet->r = r;
}

static inline double rpacket_get_tau_event (const rpacket_t * packet)
{
  return packet->tau_event;
}

static inline void rpacket_set_tau_event (rpacket_t * packet, double tau_event)
{
  packet->tau_event = tau_event;
}

static inline double rpacket_get_nu_line (const rpacket_t * packet)
{
  return packet->nu_line;
}

static inline void rpacket_set_nu_line (rpacket_t * packet, double nu_line)
{
  packet->nu_line = nu_line;
}

static inline unsigned int rpacket_get_current_shell_id (const rpacket_t * packet)
{
  return packet->current_shell_id;
}

static inline void rpacket_set_current_shell_id (rpacket_t * packet,
                                                 unsigned int current_shell_id)
{
  packet->current_shell_id = current_shell_id;
}

static inline unsigned int rpacket_get_next_line_id (const rpacket_t * packet)
{
  return packet->next_line_id;
}

static inline void rpacket_set_next_line_id (rpacket_t * packet,
                                             unsigned int next_line_id)
{
  packet->next_line_id = next_line_id;
}

static inline bool rpacket_get_last_line (const rpacket_t * packet)
{
  return packet->last_line;
}

static inline void rpacket_set_last_line (rpacket_t * packet, bool last_line)
{
  packet->last_line = last_line;
}

static inline bool rpacket_get_close_line (const rpacket_t * packet)
{
  return packet->close_line;
}

static inline void rpacket_set_close_line (rpacket_t * packet, bool close_line)
{
  packet->close_line = close_line;
}

static inline int rpacket_get_virtual_packet_flag (const rpacket_t * packet)
{
  return packet->virtual_packet_flag;
}

static inline void rpacket_set_virtual_packet_flag (rpacket_t * packet,
                                                    int virtual_packet_flag)
{
  packet->virtual_packet_flag = virtual_packet_flag;
}

static inline int rpacket_get_virtual_packet (const rpacket_t * packet)
{
  return packet->virtual_packet;
}

static inline void rpacket_set_virtual_packet (rpacket_t * packet,
                                               int virtual_packet)
{
  packet->virtual_packet = virtual_packet;
}

static inline double rpacket_get_d_boundary (const rpacket_t * packet)
{
  return packet->d_boundary;
}

static inline void rpacket_set_d_boundary (rpacket_t * packet, double d_boundary)
{
  packet->d_boundary = d_boundary;
}

static inline double rpacket_get_d_electron (const rpacket_t * packet)
{
  return packet->d_electron;
}

static inline void rpacket_set_d_electron (rpacket_t * packet, double d_electron)
{
  packet->d_electron = d_electron;
}

static inline double rpacket_get_d_line (const rpacket_t * packet)
{
  return packet->d_line;
}

static inline void rpacket_set_d_line (rpacket_t * packet, double d_line)
{
  packet->d_line = d_line;
}

static inline int rpacket_get_next_shell_id (const rpacket_t * packet)
{
  return packet->next_shell_id;
}

static inline void rpacket_set_next_shell_id (rpacket_t * packet, int next_shell_id)
{
  packet->next_shell_id = next_shell_id;
}

static inline rpacket_status_t rpacket_get_status (const rpacket_t * packet)
{
  return packet->status;
}

static inline void rpacket_set_status (rpacket_t * packet, rpacket_status_t status)
{
  packet->status = status;
}

static inline int rpacket_get_id (const rpacket_t * packet)
{
  return packet->id;
}

static inline void rpacket_set_id (rpacket_t * packet, int id)
{
  packet->id = id;
}

static inline void rpacket_reset_tau_event (rpacket_t * packet, rk_state *mt_state)
{
  rpacket_set_tau_event (packet, -log (rk_double (mt_state)));
}

tardis_error_t rpacket_init (rpacket_t * packet, storage_model_t * storage,
                             int packet_index, int virtual_packet_flag, double * chi_bf_tmp_partial);

/* New getter and setter methods for continuum implementation */

static inline void rpacket_set_d_continuum (rpacket_t * packet, double d_continuum)
{
  packet->d_cont = d_continuum;
}

static inline double rpacket_get_d_continuum (const rpacket_t * packet)
{
  return packet->d_cont;
}

static inline void rpacket_set_chi_electron (rpacket_t * packet, double chi_electron)
{
  packet->chi_th = chi_electron;
}

static inline double rpacket_get_chi_electron (const rpacket_t * packet)
{
  return packet->chi_th;
}

static inline void rpacket_set_chi_continuum (rpacket_t * packet, double chi_continuum)
{
  packet->chi_cont = chi_continuum;
}

static inline double rpacket_get_chi_continuum (const rpacket_t * packet)
{
  return packet->chi_cont;
}

static inline void rpacket_set_chi_freefree (rpacket_t * packet, double chi_freefree)
{
  packet->chi_ff = chi_freefree;
}

static inline double rpacket_get_chi_freefree (const rpacket_t * packet)
{
  return packet->chi_ff;
}

static inline void rpacket_set_chi_boundfree (rpacket_t * packet, double chi_boundfree)
{
  packet->chi_bf = chi_boundfree;
}

static inline double rpacket_get_chi_boundfree (const rpacket_t * packet)
{
  return packet->chi_bf;
}

static inline unsigned int rpacket_get_current_continuum_id (const rpacket_t * packet)
{
  return packet->current_continuum_id;
}

static inline void rpacket_set_current_continuum_id (rpacket_t * packet, unsigned int current_continuum_id)
{
  packet->current_continuum_id = current_continuum_id;
}

static inline void rpacket_set_macro_atom_activation_level (rpacket_t * packet, unsigned int activation_level)
{
  packet->macro_atom_activation_level = activation_level;
}

static inline unsigned int rpacket_get_macro_atom_activation_level (const rpacket_t * packet)
{
  return packet->macro_atom_activation_level;
}

#endif // TARDIS_RPACKET_H
