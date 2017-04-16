#include <stdlib.h>
#include <string.h>


#define ASSERT_EQUAL_FLOAT(exp, got)                                         \
    do {                                                                     \
        if ((exp) == (got)) {                                                \
            printf("\n[PASS] %s:%s():%d\n", __FILE__, __func__, __LINE__);   \
        }                                                                    \
        else {                                                               \
            printf("\n[FAIL] %s:%s():%d\n", __FILE__, __func__, __LINE__);   \
            printf("  [TST] ASSERT_EQUAL_STR(" #exp ", " #got ")\n");        \
            printf("  [EXP] %lf\n", (exp));                                   \
            printf("  [GOT] %lf\n", (got));                                   \
        }                                                                    \
    } while (0);
