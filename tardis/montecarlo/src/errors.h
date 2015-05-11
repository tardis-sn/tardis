#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

typedef enum
{
  TARDIS_ERROR_OK = 0,
  TARDIS_ERROR_BOUNDS_ERROR = 1,
  TARDIS_ERROR_COMOV_NU_LESS_THAN_NU_LINE = 2
} tardis_error_t;

