#ifdef WITHOPENMP
  #include <omp.h>
#else
int omp_get_num_threads(){
  return 1;
}
#endif
