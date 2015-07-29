#ifndef TARDIS_STATUS_H
#define TARDIS_STATUS_H

#ifdef __clang__
#define INLINE extern inline
#else
#define INLINE inline
#endif

typedef enum
{
  TARDIS_ERROR_OK = 0,
  TARDIS_ERROR_BOUNDS_ERROR = 1,
  TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE = 2
} tardis_error_t;

typedef enum
{
  TARDIS_PACKET_STATUS_IN_PROCESS = 0,
  TARDIS_PACKET_STATUS_EMITTED = 1,
  TARDIS_PACKET_STATUS_REABSORBED = 2
} rpacket_status_t;

typedef enum
{
  CONTINUUM_OFF = 0,
  CONTINUUM_ON = 1,
} ContinuumProcessesStatus;

typedef enum
{
  FREE_FREE_OFF = 0,
  FREE_FREE_ON = 1,
} FreeFreeStatus;

typedef enum
{
  BB_EMISSION = -1,
  BF_EMISSION = -2,
  FF_EMISSION = -3,
  COLL_EXCITATION = 0,
  COLL_IONIZATION = 1,
  KPACKET_CREATION = 2
} next_interaction2process;

typedef enum
{
  EXCITATION_ENERGY = 0,
  IONIZATION_ENERGY = 1,
  THERMAL_ENERGY = 2
} e_packet_type;

#endif // TARDIS_STATUS_H