#include <inttypes.h>
#include <stdio.h>
#include <unistd.h>

#define STATUS_FORMAT "\r\033[2K\t[%" PRId64 "%%] Packets(finished/total): %" PRId64 "/%" PRId64

static inline void
print_progress (const int64_t current, const int64_t total)
{
  if (isatty(fileno(stderr)))
    {
      fprintf(stderr, STATUS_FORMAT,
              current * 100 / total,
              current,
              total);
    }
}
