#define OPENMP
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
		     int64_t openmp_threads,
		     int64_t virtual_packet_flag)
{
    //printf("montecarlo_main_loop - storage prt %p\n", storage);
    //printf("montecarlo_main_loop -storage->no_of_lines: %d\n",storage->no_of_lines);
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
	int64_t *big_spectrum_virt_nu_todo_index = NULL;
	double * big_spectrum_virt_nu_todo_value = NULL;
	int64_t *big_line_lists_j_blues_todo_index =NULL;
	double * big_line_lists_j_blues_todo_value=NULL;
	int64_t ip;
    //printf("montecarlo_parallel_loop - storage prt %p", storage);

	//omp_set_num_threads(20);
	omp_set_num_threads(openmp_threads);
	
#pragma omp parallel shared(storage,openmp_threads, virtual_packet_flag, big_packet, big_spectrum_virt_nu_todo_index,big_spectrum_virt_nu_todo_value, big_line_lists_j_blues_todo_index,big_line_lists_j_blues_todo_value) ,private(ip) ,default(none)
	{
		int64_t npacket;
		int64_t reabsorbed = 0;
		int64_t thread_num = 0;
		rpacket_t *packet = NULL;

		npacket = storage->no_of_packets;


#pragma single
		{
			int64_t num_threads = omp_get_num_threads();
			big_packet =
			    (rpacket_t *) safe_malloc(sizeof(rpacket_t) * num_threads);

			big_spectrum_virt_nu_todo_index = (int64_t* )safe_malloc(sizeof(int64_t) * num_threads * MAX_TODO_LEN);
			big_spectrum_virt_nu_todo_value = (double* )safe_malloc(sizeof(int64_t) * num_threads * MAX_TODO_LEN);

			big_line_lists_j_blues_todo_index = (int64_t* )safe_malloc(sizeof(int64_t) * num_threads * MAX_TODO_LEN);
			big_line_lists_j_blues_todo_value = (double* )safe_malloc(sizeof(int64_t) * num_threads * MAX_TODO_LEN);


		}
#pragma omp for
		for (ip = 0; ip <= npacket; ip++) {
			thread_num = omp_get_thread_num();
			packet = (rpacket_t *) (big_packet + thread_num);
			rpacket_set_id(packet, ip);
			packet->spectrum_virt_nu = NULL;
			packet->line_lists_j_blues = NULL;

    assign_packet_todo(packet, 0,big_spectrum_virt_nu_todo_index, big_spectrum_virt_nu_todo_value, big_line_lists_j_blues_todo_index, big_line_lists_j_blues_todo_value);

		    //printf("montecarlo_parallel_loop - storage prt %p\n", storage);
			rpacket_init(packet, storage, ip, virtual_packet_flag);

			if (!(ip % (npacket / 20))) {
				printf("%d\n", ip);
			}
			if (virtual_packet_flag > 0) {
				reabsorbed =
				    montecarlo_one_packet(storage, packet,
							  -1);
			} else {
				reabsorbed =
				    montecarlo_one_packet(storage, packet, 0);
			}

			move_data_packet2storage(storage, packet, ip,
						 reabsorbed);
		}
#pragma single
		{
			free(big_packet);
			big_packet = NULL;
			//free(big_thread_spectrum_virt_nu);
		    //big_thread_spectrum_virt_nu = NULL;
			//free(big_thread_line_lists_j_blues);
			//big_thread_line_lists_j_blues = NULL;
		}

	}
#else
	printf
	    ("CRITICAL - Tardis was built without openmp support. Fallback to none parallel run");
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
	int64_t ip =0;

	int64_t *big_spectrum_virt_nu_todo_index = NULL;
	double * big_spectrum_virt_nu_todo_value = NULL;
	int64_t *big_line_lists_j_blues_todo_index =NULL;
	double * big_line_lists_j_blues_todo_value=NULL;

	npacket = storage->no_of_packets;
	packet = (rpacket_t*) safe_malloc(sizeof(rpacket_t));
	packet->spectrum_virt_nu = NULL;
	packet->line_lists_j_blues = NULL;

    		big_spectrum_virt_nu_todo_index = (int64_t* )safe_malloc(sizeof(int64_t) *  MAX_TODO_LEN);
			big_spectrum_virt_nu_todo_value = (double* )safe_malloc(sizeof(int64_t) *  MAX_TODO_LEN);

			big_line_lists_j_blues_todo_index = (int64_t* )safe_malloc(sizeof(int64_t) *  MAX_TODO_LEN);
			big_line_lists_j_blues_todo_value = (double* )safe_malloc(sizeof(int64_t) *  MAX_TODO_LEN);



	for (ip = 0; ip <= npacket; ip++) {
		storage->current_packet_id = ip;
        assign_packet_todo(packet, 0,big_spectrum_virt_nu_todo_index, big_spectrum_virt_nu_todo_value, big_line_lists_j_blues_todo_index, big_line_lists_j_blues_todo_value);
		if (!(ip % (npacket / 20))) {
			printf("%d\n", ip);
		}
		//printf("montecarlo_serial_loop - storage prt %p\n", storage);
        //printf("montecarlo_serial_loop -storage->no_of_lines: %d\n",storage->no_of_lines);
		rpacket_init(&packet, storage, ip, virtual_packet_flag);
		if (virtual_packet_flag > 0) {
			reabsorbed =
			    montecarlo_one_packet(storage, packet, -1);
		}
		reabsorbed = montecarlo_one_packet(storage, packet, 0);
		move_data_packet2storage(storage, packet, ip, reabsorbed);
	}
	free(packet);
	packet = NULL;
}

void
move_data_packet2storage(storage_model_t * storage,rpacket_t  * packet,
			 int64_t ip, int64_t reabsorbed)
{
	int64_t i =0;
	int64_t ii = 0;

	int64_t svntd_len = packet->spectrum_virt_nu_todo_len;
	int64_t lljbtd_len = packet->line_lists_j_blues_todo_len;

	storage->output_nus[ip] = rpacket_get_nu(packet);
	if (reabsorbed == 1) {
		storage->output_energies[ip] = -rpacket_get_energy(packet);
	} else {
		storage->output_energies[ip] = rpacket_get_energy(packet);
	}
	storage->last_line_interaction_in_id[ip] =
	    rpacket_get_last_line_interaction_in_id(packet);
	storage->last_line_interaction_shell_id[ip] =
	    rpacket_get_last_line_interaction_shell_id(packet);
	storage->last_interaction_type[ip] =
	    rpacket_get_last_interaction_type(packet);
	storage->last_line_interaction_out_id[ip] =
	    rpacket_get_last_line_interaction_out_id(packet);


    if (svntd_len > 0)
    {
    for (i = 0; i < svntd_len; i++) {
        ii = packet->spectrum_virt_nu_todo_index[i];
        #pragma omp atomic
        storage->spectrum_virt_nu[ii] += packet->spectrum_virt_nu_todo_value[i];
        }
    }

    if (svntd_len > 0)
    {
    for (i = 0; i < lljbtd_len; i++) {
        ii = packet->line_lists_j_blues_todo_index[i];
        #pragma omp atomic
        storage->line_lists_j_blues[ii] += packet->line_lists_j_blues_todo_value[i];
        }
    }
}

void
assign_packet_todo(rpacket_t * packet, int64_t thread_num, int64_t * big_spectrum_virt_nu_todo_index,
double * big_spectrum_virt_nu_todo_value, int64_t * big_line_lists_j_blues_todo_index, double * big_line_lists_j_blues_todo_value
){
packet->spectrum_virt_nu_todo_len = 0;
packet->spectrum_virt_nu_todo_index = (int64_t*) big_spectrum_virt_nu_todo_index + (thread_num * MAX_TODO_LEN);
packet->spectrum_virt_nu_todo_value = (double*) big_spectrum_virt_nu_todo_value +(thread_num * MAX_TODO_LEN);
packet->line_lists_j_blues_todo_len = 0;
packet->line_lists_j_blues_todo_index = (int64_t*) big_line_lists_j_blues_todo_index +(thread_num * MAX_TODO_LEN);
packet->line_lists_j_blues_todo_value = (double*) big_line_lists_j_blues_todo_value +(thread_num * MAX_TODO_LEN);
}

