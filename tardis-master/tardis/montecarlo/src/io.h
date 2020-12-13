#define _POSIX_C_SOURCE 1

#include <unistd.h>
#include <inttypes.h>
#include <stdio.h>

#define STATUS_FORMAT "\r\033[2K\t[%" PRId64 "%%] Packets(finished/total): %" PRId64 "/%" PRId64
#define STATUS_FORMAT_FI "\r\033[2K\t[%" PRId64 "%%] Bins(finished/total): %" PRId64 "/%" PRId64

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

static inline void
print_progress_fi (const int64_t current, const int64_t total)
{
  if (isatty(fileno(stderr)))
    {
      fprintf(stderr, STATUS_FORMAT_FI,
              current * 100 / total,
              current,
              total);
    }
}
