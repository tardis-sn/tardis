#ifndef TARDIS_CMONTECARLO_MAINLOOP_H
#define TARDIS_CMONTECARLO_MAINLOOP_H

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include "cmontecarlo.h"

#ifdef __clang__
#define INLINE extern inline
#else
#define INLINE inline
#endif

#define MAX_TODO_LEN 10000

#define INFO(msg) \
    fprintf(stderr, "info: %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, "%s", msg);

static void *safe_malloc(size_t n);

void montecarlo_main_loop(storage_model_t * storage,
			  int64_t openmp_threads, int64_t virtual_packet_flag);

void montecarlo_parallel_loop(storage_model_t * storage, int64_t openmp_threads,
			      int64_t virtual_mode);

void montecarlo_serial_loop(storage_model_t * storage, int64_t openmp_threads,
			    int64_t virtual_mode);

void
move_data_packet2storage(storage_model_t * storage, rpacket_t * packet,
			 int64_t ip, int64_t reabsorbed);

void
assign_packet_todo(rpacket_t * packet, int64_t threads_num,
		   int64_t * big_spectrum_virt_nu_todo_index,
		   double *big_spectrum_virt_nu_todo_value,
		   int64_t * big_line_lists_j_blues_todo_index,
		   double *big_line_lists_j_blues_todo_value);

#endif				/*TARDIS_CMONTECARLO_MAINLOOP_H */
