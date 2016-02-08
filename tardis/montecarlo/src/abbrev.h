#ifndef ABBREV_H
#define ABBREV_H

#define FAIL() abort()

/* safe malloc */

static inline void* safe_malloc(size_t size) {
  void *mem = malloc(size);
  if (mem == NULL && size != 0) FAIL();
  return mem;
}

static inline void* safe_realloc(void *ptr, size_t size) {
  void *mem = realloc(ptr, size);
  if (mem == NULL && size != 0) FAIL();
  return mem;
}

#endif /* ABBREV_H */
