
from numba import prange, njit
from tardis.montecarlo.montecarlo_numba.rpacket import RPacket

def montecarlo_radial1d(model, plasma, runner):
    storage_model_kwargs = {'packet_nus', packet_mus, packet_energies, 
    output_nus, output_energies, no_of_packets, no_of_shells, 
    r_inner, r_outer, v_inner, time_explosion, electron_densities, line_list_nu, line_lists_tau_sobolevs, line_lists_tau_sobolevs_nd, 
    no_of_lines, no_of_edges, line_interaction_id, 
    inverse_sigma_thomson}
    pass
        #montecarlo.montecarlo_radial1d(
        #    model, plasma, self,
        #    virtual_packet_flag=no_of_virtual_packets,
        #    nthreads=nthreads,
        #    last_run=last_run)
        # Workaround so that j_blue_estimator is in the right ordering
        # They are written as an array of dimension (no_of_shells, no_of_lines)
        # but python expects (no_of_lines, no_of_shells)

@njit
def montecarlo_main_loop(storage):
    for i in prange(storage.no_of_packets):
        from


void
montecarlo_main_loop(storage_model_t * storage, int64_t virtual_packet_flag, int nthreads, unsigned long seed)
{
  int64_t finished_packets = 0;
  storage->virt_packet_count = 0;
#ifdef WITH_VPACKET_LOGGING
  storage->virt_packet_nus = (double *)safe_malloc(sizeof(double) * storage->no_of_packets);
  storage->virt_packet_energies = (double *)safe_malloc(sizeof(double) * storage->no_of_packets);
  storage->virt_packet_last_interaction_in_nu = (double *)safe_malloc(sizeof(double) * storage->no_of_packets);
  storage->virt_packet_last_interaction_type = (int64_t *)safe_malloc(sizeof(int64_t) * storage->no_of_packets);
  storage->virt_packet_last_line_interaction_in_id = (int64_t *)safe_malloc(sizeof(int64_t) * storage->no_of_packets);
  storage->virt_packet_last_line_interaction_out_id = (int64_t *)safe_malloc(sizeof(int64_t) * storage->no_of_packets);
  storage->virt_array_size = storage->no_of_packets;
#endif // WITH_VPACKET_LOGGING
#ifdef WITHOPENMP
  omp_set_dynamic(0);
  if (nthreads > 0)
    {
      omp_set_num_threads(nthreads);
    }

#pragma omp parallel firstprivate(finished_packets)
    {
      rk_state mt_state;
      rk_seed (seed + omp_get_thread_num(), &mt_state);
#pragma omp master
      {
        fprintf(stderr, "Running with OpenMP - %d threads\n", omp_get_num_threads());
        print_progress(0, storage->no_of_packets);
      }
#else
      rk_state mt_state;
      rk_seed (seed, &mt_state);
      fprintf(stderr, "Running without OpenMP\n");
#endif
      int64_t chi_bf_tmp_size = (storage->cont_status) ? storage->no_of_edges : 0;
      double *chi_bf_tmp_partial = safe_malloc(sizeof(double) * chi_bf_tmp_size);

      #pragma omp for
      for (int64_t packet_index = 0; packet_index < storage->no_of_packets; ++packet_index)
        {
          int reabsorbed = 0;
          rpacket_t packet;
          rpacket_set_id(&packet, packet_index);
          rpacket_init(&packet, storage, packet_index, virtual_packet_flag, chi_bf_tmp_partial);
          if (virtual_packet_flag > 0)
            {
              reabsorbed = montecarlo_one_packet(storage, &packet, -1, &mt_state);
            }
          reabsorbed = montecarlo_one_packet(storage, &packet, 0, &mt_state);
          storage->output_nus[packet_index] = rpacket_get_nu(&packet);
          if (reabsorbed == 1)
            {
              storage->output_energies[packet_index] = -rpacket_get_energy(&packet);
            }
          else
            {
              storage->output_energies[packet_index] = rpacket_get_energy(&packet);
            }
          if ( ++finished_packets%100 == 0 )
            {
#ifdef WITHOPENMP
              // WARNING: This only works with a static sheduler and gives an approximation of progress.
              // The alternative would be to have a shared variable but that could potentially decrease performance when using many threads.
              if (omp_get_thread_num() == 0 )
                print_progress(finished_packets * omp_get_num_threads(), storage->no_of_packets);
#else
              print_progress(finished_packets, storage->no_of_packets);
#endif
            }
        }
      free(chi_bf_tmp_partial);
#ifdef WITHOPENMP
    }
#endif
  print_progress(storage->no_of_packets, storage->no_of_packets);
  fprintf(stderr,"\n");
}

