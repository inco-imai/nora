#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#include "mymurmurhash3.h"

int LOOP = 1000;
int N_KINT = 100000;
uint32_t WINDOW = 1000;

int main(int argc, char ** argv)
{
  int N_THREADS = omp_get_max_threads();

  uint32_t * seeds;
  seeds = (uint32_t*)malloc(sizeof(uint32_t)*LOOP);
  if(seeds == NULL){
    fprintf(stderr, "cannot allocate memory seeds\n");
    return 1;
  }
  int i;
  for(i=0; i<LOOP; ++i){
    seeds[i] = (uint32_t)i;
  }

  uint32_t * kmerints;
  kmerints = (uint32_t*)malloc(sizeof(uint32_t)*N_KINT);
  if(kmerints == NULL){
    fprintf(stderr, "cannot allocate memory kmerints\n");
    return 1;
  }
  uint32_t * fingerprints;
  fingerprints = (uint32_t*)malloc(sizeof(uint32_t)*N_THREADS);
  if(fingerprints == NULL){
    fprintf(stderr, "cannot allocate memory fingerprints\n");
    return 1;
  }

  for(i=0; i<N_KINT; ++i){
    kmerints[i] = (uint32_t)i;
  }

  int * counts;
  counts = (int*)malloc(sizeof(int)*WINDOW);
  if(counts == NULL){
    fprintf(stderr, "cannot allocate memory counts\n");
    return 1;
  }
  for(i=0; i<WINDOW; ++i){
    counts[i] = 0;
  }

  int j;
  for(i=0; i<N_KINT; ++i){
    #pragma omp parallel for
    for(j=0; j<LOOP; ++j){
      int thread_num = omp_get_thread_num();
      mymurmurhash3(kmerints[i], seeds[j], fingerprints+thread_num);
      //printf("kint %d, seed %d, fp %d\n",kmerints[i],seeds[j],fingerprints[thread_num]);
      #pragma omp atomic
      ++counts[fingerprints[thread_num] % WINDOW];
    }
  }

  for(i=0; i<WINDOW; ++i){
    printf("%d %d\n",i,counts[i]);
  }

  free(seeds);
  free(kmerints);
  free(fingerprints);
  free(counts);

  return 0;
}

