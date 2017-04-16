CC=gcc
CFLAGS=-shared -Wl,-soname,cmontecarlo -lm -fPIC
SOURCES=cmontecarlo.c cmontecarlo_wrappers.c randomkit/rk_isaac.c randomkit/rk_mt.c randomkit/rk_sobol.c randomkit/rk_primitive.c
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=cmontecarlo.so

all: 
	$(CC) $(CFLAGS) $(SOURCES) -o cmontecarlo.so
