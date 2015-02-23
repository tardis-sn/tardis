#define OPENMP = True
#ifdef OPENMP
#include <omp.h>
#endif
#include "cmontecarlo_mainloop.h"
#include "cmontecarlo.h"

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
	int64_t npacket;
	int64_t reabsorbed = 0;
	int64_t num_threads = 0;
	int64_t threads_num = 0;
	int64_t spectrum_virt_nu_len = 0;
	int64_t line_lists_j_blues_len = 0;
	RPacket *packet = NULL;
	RPacket *big_packet = NULL;
	double *big_thread_spectrum_virt_nu = NULL;
	double *big_thread_line_lists_j_blues = NULL;

	npacket = storage->no_of_packets;
	spectrum_virt_nu_len = storage->spectrum_virt_nu_len;
	line_lists_j_blues_len = storage->line_lists_j_blues_len;

	omp_set_num_threads(openmp_threads);

	num_threads = omp_get_num_threads();
	big_packet = (RPacket *) malloc(sizeof(RPacket) * num_threads);
	big_thread_spectrum_virt_nu =
	    (double *)malloc(sizeof(double) *
			     (spectrum_virt_nu_len * num_threads));
	big_thread_line_lists_j_blues =
	    (double *)malloc(sizeof(double) *
			     (line_lists_j_blues_len * num_threads));

#pragma omp parallel for (shared(storage,openmp_threads, virtual_packet_flag), private(reabsorbed,packet))
	for (int ip = 0; ip <= npacket; ip++) {
		threads_num = omp_get_thread_num();
		packet = (RPacket *) (big_packet + threads_num);
		rpacket_set_id(packet, ip);
		packet->thread_spectrum_virt_nu =
		    (double *)(big_thread_spectrum_virt_nu + threads_num);
		packet->thread_line_lists_j_blues =
		    (double *)(big_thread_line_lists_j_blues + threads_num);
		rpacket_init(packet, &storage, ip, virtual_packet_flag);

		if (!(npacket % (npacket / 20))) {
			printf("%d\n", npacket);
		}
		if (virtual_packet_flag > 0) {
			reabsorbed =
			    montecarlo_one_packet(&storage, &packet, -1);
		} else {
			reabsorbed =
			    montecarlo_one_packet(&storage, &packet, 0);
		}

		move_data_packet2storage(storage, packet, ip, reabsorbed);
	}
	free(packet);
	packet = NULL;


#else
	printf
	    ('CRITICAL - Tardis was build without openmp support. Fallback to none parallel run');
	montecarlo_serial_loop(storage, 0, virtual_packet_flag);


#endif

}

void
montecarlo_serial_loop(storage_model_t * storage, int64_t openmp_threads,
		       int64_t virtual_packet_flag)
{
	int64_t npacket;
	int64_t reabsorbed;
	RPacket *packet = NULL;
	int64_t spectrum_virt_nu_len = 0;
	int64_t line_lists_j_blues_len = 0;

	npacket = storage->no_of_packets;
	spectrum_virt_nu_len = storage->spectrum_virt_nu_len;
	line_lists_j_blues_len = storage->line_lists_j_blues_len;
	packet = (RPacket *) malloc(sizeof(RPacket));
	packet->thread_spectrum_virt_nu =
	    (double *)malloc(sizeof(double) * (spectrum_virt_nu_len));
	packet->line_lists_j_blues =
	    (double *)malloc(sizeof(double) * (line_lists_j_blues_len));

	for (int ip = 0; ip <= npacket; ip++) {
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
move_data_packet2storage(storage_model_t * storage, RPacket * packet,
			 int64_t ip, int64_t reabsorbed)
{
	int64_t svn_len = storage->spectrum_virt_nu_len;
	int64_t lljb_len = storage->line_lists_j_blues_len;

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

	for (int i = 0; i < svn_len; i++) {
#pragma omp atomic
		storage->spectrum_virt_nu[i] += packet->spectrum_virt_nu[i];
	}
	for (int i = 0; i < lljb_len; i++) {
#pragma omp atomic
		storage->line_lists_j_blues_[i] +=
		    packet->line_lists_j_blues_[i];
	}
}
