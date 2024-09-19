// originally written by stanaka2 at https://qiita.com/stanaka2/items/526102c6c56759c22b8f
// modified by Imai
// other refefence https://github.com/eduardlopez/quicksort-parallel
// (linked by stanaka2)


#define PARALLEL_QSORT_THRESHOLD (1000)

static inline signed long long median(signed long long a, signed long long b, signed long long c){
  if(a>b){
    if(b>c){
      return b;
    }
    else if(c>a){
      return a;
    }
    else{
      return c;
    }
  }
  else{
    if(a>c){
      return a;
    }
    else if(c>b){
      return b;
    }
    else{
      return c;
    }
  }
}

static inline void qsort_partitioning(record_t ** p_kmers, signed long long * left, signed long long * right, signed long long pivot){
  signed long long i,j;
  i=*left;
  j=*right;

  while(1){
    while(p_kmers[i]->kmer < pivot){
      ++i;
    }
    while(p_kmers[j]->kmer > pivot){
      --j;
    }
    if(i>=j){
      break;
    }
    
    record_t * tmp = p_kmers[i];
    p_kmers[i] = p_kmers[j];
    p_kmers[j] = tmp;
 
    ++i;
    --j;
  }

  *left = i;
  *right = j;
}

void single_thread_qsort(record_t ** p_kmers, signed long long left, signed long long right){
  if(left < right){
    signed long long i,j;
    i=left;
    j=right;

    signed long long pivot;
    pivot = median(p_kmers[i]->kmer, p_kmers[j]->kmer, p_kmers[(i+j)/2]->kmer);

    qsort_partitioning(p_kmers, &i, &j, pivot);

    single_thread_qsort(p_kmers, left, i-1);
    single_thread_qsort(p_kmers, j+1, right);
  }
}

void parallel_task_qsort_internal(record_t ** p_kmers, signed long long left, signed long long right){
  signed long long len = right-left+1;
  if(len < PARALLEL_QSORT_THRESHOLD){
    single_thread_qsort(p_kmers, left, right);
    return;
  }

  signed long long i,j;
  i=left;
  j=right;

  signed long long pivot = median(p_kmers[i]->kmer, p_kmers[j]->kmer, p_kmers[(i+j)/2]->kmer);

  qsort_partitioning(p_kmers, &i, &j, pivot);
  
  #pragma omp task
  parallel_task_qsort_internal(p_kmers,left,j);
  #pragma omp task
  parallel_task_qsort_internal(p_kmers,i,right);
}

void parallel_task_qsort(record_t ** p_kmers, signed long long left, signed long long right){
  #pragma omp parallel
  {
    #pragma omp single nowait
    {
      parallel_task_qsort_internal(p_kmers, left, right);
    }
  }
}

