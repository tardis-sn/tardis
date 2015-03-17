Q=@
ifdef V
Q=
endif
# Directory variables.
INC_DIR = .
# Flags passed to the C++ compiler.
CFLAGS += -g -Wall -Wextra -Wno-nonnull
# Source files to compile.
EXCLUDES = cmontecarlo.c
SRCS = $(filter-out $(EXCLUDES), $(wildcard *.c))
EXES = $(patsubst %.c,%,$(SRCS))
all : $(EXES)
clean :
	$(Q)rm -f $(EXES)
# Executable targets.
%: %.c test_macros.h
	$(Q)$(CC) -I$(INC_DIR) $(CFLAGS) $< -o $@ $(LIBS)
