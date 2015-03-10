Q=@
ifdef V
Q=
endif

# Directory variables.
INC_DIR = .

CFLAGS += -g -Wall -Wextra -Wno-nonnull

SRCS = $(wildcard *.c)
EXES = $(patsubst %.c,%,$(SRCS))

all : $(EXES)

clean :
	$(Q)rm -f $(EXES)

%: %.c test_macros.h
	$(Q)$(CC) -I$(INC_DIR) $(CFLAGS) $< -o $@ $(LIBS)

