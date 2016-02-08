#ifndef ABBREV_H
#define ABBREV_H

#include <stdlib.h>

#define FAIL() abort()

/**
 * @brief safe malloc; checks for NULL and aborts if encountered
 */
static inline void* safe_malloc(size_t size) {
  void *mem = malloc(size);
  if (mem == NULL && size != 0) FAIL();
  return mem;
}

/**
 * @brief safe realloc; checks for NULL and aborts if encountered
 */
static inline void* safe_realloc(void *ptr, size_t size) {
  void *mem = realloc(ptr, size);
  if (mem == NULL && size != 0) FAIL();
  return mem;
}

#endif /* ABBREV_H */
