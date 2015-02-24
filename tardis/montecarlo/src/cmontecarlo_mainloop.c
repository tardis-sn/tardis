#define OPENMP = True
#ifdef OPENMP
#include <omp.h>
#endif
#include "cmontecarlo_mainloop.h"
#include "cmontecarlo.h"

static void *safe_malloc(size_t n)
{
	void *p = malloc(n);
	if (!p) {
		perror("CRITICAL - Error in malloc");
		abort();
	}
	return p;
}

void
montecarlo_main_loop(storage_model_t * storage,
		     int64_t virtual_mode, int64_t openmp_threads,
		     int64_t virtual_packet_flag)
{
	if (openmp_threads == 0) {
		montecarlo_serial_loop(storage, openmp_threads,
				       virtual_packet_flag);
	} else {
		montecarlo_parallel_loop(storage, openmp_threads,
					 virtual_packet_flag);
	}
}

void
montecarlo_parallel_loop(storage_model_t * storage, int64_t openmp_threads,
			 int64_t virtual_packet_flag)
{
#ifdef OPENMP
	rpacket_t *big_packet = NULL;
	double *big_thread_spectrum_virt_nu = NULL;
	double *big_thread_line_lists_j_blues = NULL;
	int64_t ip;

#pragma omp parallel shared(storage,openmp_threads, virtual_packet_flag, big_packet, big_thread_spectrum_virt_nu, big_thread_line_lists_j_blues),privat(ip) ,default(none)
	{
		int64_t npacket;
		int64_t reabsorbed = 0;
		int64_t thread_num = 0;
		int64_t spectrum_virt_nu_len = 0;
		int64_t line_lists_j_blues_len = 0;
		rpacket_t *packet = NULL;

		npacket = storage->no_of_packets;
		spectrum_virt_nu_len = storage->spectrum_virt_nu_len;
		line_lists_j_blues_len = storage->line_lists_j_blues_len;

		omp_set_num_threads(openmp_threads);

#pragma single
		{
			int64_t num_threads = omp_get_num_threads();
			big_packet =
			    (rpacket_t *) safe_malloc(sizeof(rpacket_t) * num_threads);
			big_thread_spectrum_virt_nu =
			    (double *)safe_malloc(sizeof(double) *
					     (spectrum_virt_nu_len *
					      num_threads));
			big_thread_line_lists_j_blues =
			    (double *)safe_malloc(sizeof(double) *
					     (line_lists_j_blues_len *
					      num_threads));

		}
#pragma omp for
		for (ip = 0; ip <= npacket; ip++) {
			thread_num = omp_get_thread_num();
			packet = (rpacket_t *) (big_packet + thread_num);
			rpacket_set_id(packet, ip);
			packet->thread_spectrum_virt_nu =
			    (double *)(big_thread_spectrum_virt_nu +
				       thread_num);
			packet->thread_line_lists_j_blues =
			    (double *)(big_thread_line_lists_j_blues +
				       thread_num);
			rpacket_init(packet, &storage, ip, virtual_packet_flag);

			if (!(npacket % (npacket / 20))) {
				printf("%d\n", npacket);
			}
			if (virtual_packet_flag > 0) {
				reabsorbed =
				    montecarlo_one_packet(&storage, &packet,
							  -1);
			} else {
				reabsorbed =
				    montecarlo_one_packet(&storage, &packet, 0);
			}

			move_data_packet2storage(storage, packet, ip,
						 reabsorbed);
		}
#pragma single
		{
			free(big_packet);
			big_packet = NULL;
			free(big_thread_spectrum_virt_nu);
			big_thread_spectrum_virt_nu = NULL;
			free(big_thread_line_lists_j_blues);
			big_thread_line_lists_j_blues = NULL;
		}

	}
#else
	printf
	    ('CRITICAL - Tardis was built without openmp support. Fallback to none parallel run');
	montecarlo_serial_loop(storage, 0, virtual_packet_flag);

#endif

}

void
montecarlo_serial_loop(storage_model_t * storage, int64_t openmp_threads,
		       int64_t virtual_packet_flag)
{
	int64_t npacket;
	int64_t reabsorbed;
	rpacket_t *packet = NULL;
	int64_t spectrum_virt_nu_len = 0;
	int64_t line_lists_j_blues_len = 0;
	int64_t ip =0;

	npacket = storage->no_of_packets;
	spectrum_virt_nu_len = storage->spectrum_virt_nu_len;
	line_lists_j_blues_len = storage->line_lists_j_blues_len;
	packet = (rpacket_t*) safe_malloc(sizeof(rpacket_t));
	packet->spectrum_virt_nu =
	    (double*)safe_malloc(sizeof(double) * (spectrum_virt_nu_len));
	packet->line_lists_j_blues =
	    (double*)safe_malloc(sizeof(double) * (line_lists_j_blues_len));

	for (ip = 0; ip <= npacket; ip++) {
		storage->current_packet_id = ip;
		if (!(npacket % (npacket / 20))) {
			printf("%d\n", npacket);
		}
		rpacket_init(&packet, &storage, ip, virtual_packet_flag);
		if (virtual_packet_flag > 0) {
			reabsorbed =
			    montecarlo_one_packet(&storage, &packet, -1);
		}
		reabsorbed = montecarlo_one_packet(&storage, &packet, 0);
		move_data_packet2storage(storage, packet, ip, reabsorbed);
	}
	free(packet);
	packet = NULL;
}

void
move_data_packet2storage(storage_model_t * storage,rpacket_t  * packet,
			 int64_t ip, int64_t reabsorbed)
{
	int64_t svn_len = storage->spectrum_virt_nu_len;
	int64_t lljb_len = storage->line_lists_j_blues_len;
	int64_t i =0;

	storage->output_nus[ip] = rpacket_get_nu(&packet);
	if (reabsorbed == 1) {
		storage->output_energies[ip] = -rpacket_get_energy(&packet);
	} else {
		storage->output_energies[ip] = rpacket_get_energy(&packet);
	}
	storage->last_line_interaction_in_id[ip] =
	    rpacket_get_last_line_interaction_in_id(packet);
	storage->last_line_interaction_shell_id[ip] =
	    rpacket_get_last_line_interaction_shell_id(packet);
	storage->last_interaction_type[ip] =
	    rpacket_get_last_interaction_type(packet);
	storage->last_line_interaction_out_id[ip] =
	    rpacket_get_last_line_interaction_out_id(packet);

	for (i = 0; i < svn_len; i++) {
#pragma omp atomic
		storage->spectrum_virt_nu[i] += packet->spectrum_virt_nu[i];
	}
	for (i = 0; i < lljb_len; i++) {
#pragma omp atomic
		storage->line_lists_j_blues[i] +=
		    packet->line_lists_j_blues[i];
	}
}
