#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include <limits.h>
#include <stdint.h>
//#include "b_heikki.h"
#include "b_heikki.dynamic.h"
#include "mymurmurhash3.h"

unsigned long long K = 14; // kmer's k. K <= 32
unsigned long long K_SUGAR = 0;
double T_PIDENT = 0.75; // XXX
int L_CHECK = 250; // for error rate check in b_heikki
int L_CHUNK = 2000; // chunk size used by dynamic b_heikki
int H_SIZE = 500;// XXX

int LSEQ = 65535;//XXX

int MAX_DEPTH = 200;// XXX
int MAX_CAND;// = LSEQ*MAX_DEPTH;// XXX
int THRE_J = 10;// TODO
int THRE_J_UPPER = 10000;// TODO
int L_CAN_MAX=150;//XXX
int TOO_MANY_KMER=180;//XXX
unsigned long long PM=100;
unsigned long long DM=100;
double V_ACCURACY = 0.8;

//int MIN_OL = 50;
int EPS = 0;
int DEBUG_ALI_PRINT=0;

int MIN_ALI_SCORE = 300;
double MIN_ALI_FRAC = 0.0;//XXX
double P_SIMILAR = 0.95;

unsigned long long s_kmers;

char mat[128];
unsigned long long to_code[128];

typedef struct record_t{
  unsigned long id;
  unsigned long long kmer;
  unsigned long long offset;
  char strand;
  int f_valid;
}record_t;

typedef struct kv_t{
  int key;
  int value;
}kv_t;

#include "myqsort.h"

typedef struct pos_t{
  unsigned long long stt;
  unsigned long long end;
}pos_t;

typedef struct bfmt7_t{
  //char qname[256];
  //char sname[256];
  unsigned long long qid;
  unsigned long long sid;
  char * qseq;
  char * sseq;
  char * ali;
  int sstt,send,qstt,qend;
  int aliscore;
  double evalue,pident;
  int l_data;
  int qlen;
  int slen;
  int sstt2,send2,qstt2,qend2;
  //char kmerseq[256];
  //char bbuf[256];
  //char abuf[256];
  //char alibuf[256];
  //char subbuf[256];
}bfmt7_t;

typedef struct name_t{
  char name[256];
}name_t;

void malloc_bfmt7array(bfmt7_t* array, int len){
  int i;
  for(i=0; i<len; ++i){
    array[i].qseq = (char*)malloc(sizeof(char)*LSEQ);
    if(array[i].qseq == NULL){
      fprintf(stderr, "cannot allocate memory: array[%d].qseq",i);
    }
    array[i].sseq = (char*)malloc(sizeof(char)*LSEQ);
    if(array[i].sseq == NULL){
      fprintf(stderr, "cannot allocate memory: array[%d].sseq",i);
    }
    array[i].ali = (char*)malloc(sizeof(char)*LSEQ);
    if(array[i].ali == NULL){
      fprintf(stderr, "cannot allocate memory: array[%d].ali",i);
    }
  }
}

void free_bfmt7array(bfmt7_t* array, int len){
  int i;
  for(i=0; i<len; ++i){
    free(array[i].qseq);
    free(array[i].sseq);
    free(array[i].ali);
  }
}

void chomp2(char * s){
  int len = strlen(s);
  (s[len-1] == '\n') ? s[len-1] = '\0' : fprintf(stderr, "strange str %s\nnot end with newline (too long line? > %d)", s,LSEQ); // chomp2
}

void reversestr(char * str, int len){
  //int len = strlen(str);
  int loop=len/2;
  int i;
  char tmp;
  for(i=0; i<loop; ++i){
    // swap
    tmp = str[i];
    str[i] = str[len-1-i];
    str[len-1-i] = tmp;
  }
  str[len] = '\0';
}

void reversecomplementstr(char * str, int len){
  //int len = strlen(str);
  int loop=len/2;
  int i;
  char tmp;
  for(i=0; i<loop; ++i){
    // swap
    tmp = str[i];
    str[i] = str[len-1-i];
    str[len-1-i] = tmp;
  }
  str[len] = '\0';
  for(i=0; i<len; ++i){
    str[i] = mat[(int)str[i]];
  }
}
void my_strcpy(char * dest, char * source, int len){
  int i;
  for(i=0; i<len; ++i){
    dest[i] = source[i];
  }
  dest[len] = '\0';
}

unsigned long long ntuple_code(const char * str, int stt, int n){
  if(n>32){
    fprintf(stderr, "nruplw_code: cannot handle %d(>32)-mer\n",n);
    exit(1);
  }
  unsigned long long ret=0;
  int i;
  for(i=0; i<n; ++i){
    ret |= (to_code[(int)str[stt+i]] << (2*(n-i-1)));
  }
  return ret;
}

double gettimeofday_sec(){
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec * 1e-6;
}

int opt_overhangingonly=0;
int opt_count_for_gnuplot=0;

void kmerize(char * reads, unsigned long long stt, unsigned long long end, record_t * kmers, int id, char strand);
int compare_kmer_records(const void *a, const void *b);
void make_sketch(uint32_t * ske, unsigned long long lb);
int get_Jaccard_similarity_like(unsigned long long i, unsigned long long j, uint32_t * sketches);
int f_hit(unsigned long long target, unsigned long long * candidates, int len);

int opt_old=0;
int opt_print_Ja=0;
int opt_DEBUG=0;
void checker(char mark, char * orig, char * ali);

int opt_print_Jmatrix=0;
kv_t * Jmat;
int compare_kv(const void *a, const void *b);

int opt_single=0;

int main(int argc, char * argv[]){
  mat[(int)'A'] = 'T';
  mat[(int)'a'] = 'T';
  mat[(int)'C'] = 'G';
  mat[(int)'c'] = 'G';
  mat[(int)'G'] = 'C';
  mat[(int)'g'] = 'C';
  mat[(int)'T'] = 'A';
  mat[(int)'t'] = 'A';
  mat[(int)'N'] = 'N';
  mat[(int)'n'] = 'N';
  mat[(int)'-'] = '-';
  to_code[(int)'A'] = 0ull;
  to_code[(int)'a'] = 0ull;
  to_code[(int)'C'] = 1ull;
  to_code[(int)'c'] = 1ull;
  to_code[(int)'G'] = 2ull;
  to_code[(int)'g'] = 2ull;
  to_code[(int)'T'] = 3ull;
  to_code[(int)'t'] = 3ull;
  to_code[(int)'N'] = 0ull;// XXX
  to_code[(int)'n'] = 0ull;

  int hitnum=0;
  {
    int result;
    while((result=getopt(argc,argv,"K:H:t:T:Om:M:cP:C:dS:JDsI:o")) != -1){
      switch(result){
        case 'K':
          K = (unsigned long long)atoi(optarg);
          if(K > 16 || K < 1){
            fprintf(stderr, "K must be 1 <= K <= 16 .\n");
            return 1;
          }
          hitnum+=2;
          break;
        case 'S':
          K_SUGAR = (unsigned long long)atoi(optarg);
          if(K_SUGAR > 10 || K_SUGAR < 0){
            fprintf(stderr, "K_SUGAR must be 0 <= K_SUGAR <= 10 .\n");
            return 1;
          }
          hitnum+=2;
          break;
        case 'H':
          H_SIZE = atoi(optarg);
          if(H_SIZE < 0){
            fprintf(stderr, "H must be >= 0. (0: don't use minhash)\n");
            return 1;
          }
          hitnum+=2;
          break;
        case 't':
          THRE_J = atoi(optarg);
          if(THRE_J < 1){
            fprintf(stderr, "THRE_J must be >= 1.\n");
            return 1;
          }
          hitnum+=2;
          break;
        case 'T':
          THRE_J_UPPER = atoi(optarg);
          if(THRE_J_UPPER < 1){
            fprintf(stderr, "THRE_J_UPPER must be >= 1.\n");
            return 1;
          }
          if(THRE_J_UPPER < THRE_J){
            fprintf(stderr, "THRE_J_UPPER must be larger than or equal to THRE_J.\n");
            return 1;
          }
          hitnum+=2;
          break;
        case 'O':
          opt_overhangingonly=1;
          ++hitnum;
          break;
        case 'm':
          MIN_ALI_FRAC = atof(optarg);
          if(MIN_ALI_FRAC < 0.0 || MIN_ALI_FRAC > 1.0){
            fprintf(stderr, "MIN_ALI_FRAC must be in [0.0 .. 1.0]. default: %lf\n",MIN_ALI_FRAC);
            return 1;
          }
          hitnum+=2;
          break;
        case 'M':
          TOO_MANY_KMER = atoi(optarg);
          if(TOO_MANY_KMER < 1){
            fprintf(stderr, "TOO_MANY_KMER must be >= 1. default: %d\n",TOO_MANY_KMER);
            return 1;
          }
          hitnum+=2;
          break;
        case 'c':
          opt_count_for_gnuplot=1;
          ++hitnum;
          break;
        case 'P':
          PM = (unsigned long long)atoi(optarg);
          if(PM > 1000){
            fprintf(stderr, "WARNING: large PM: %llu\n",PM);
            //return 1;
          }
          DM = PM; // TODO
          hitnum+=2;
          break;
        case 'C':
          L_CHUNK = atoi(optarg);
          hitnum+=2;
          if(L_CHUNK <= 0){
            fprintf(stderr, "-C (L_CHUNK) must be > 0\n");
            exit(1);
          }
          break;
        case 'd':
          opt_old=1;
          ++hitnum;
          break;
        case 'J':
          opt_print_Ja=1;
          ++hitnum;
          break;
        case 'D':
          opt_DEBUG=1;
          ++hitnum;
          break;
        case 's': // similarity
          opt_print_Jmatrix=1;
          ++hitnum;
          break;
        case 'I':
          T_PIDENT = atof(optarg);
          if(T_PIDENT < 0.0 || T_PIDENT > 1.0){
            fprintf(stderr, "T_PIDENT must be in [0.0 .. 1.0]. default: %lf\n",T_PIDENT);
            return 1;
          }
          hitnum+=2;
          break;
        case 'o':
          opt_single=1;
          ++hitnum;
          break;
        case '?':
          printf("humei\n");
          break;
        default:
          break;
      }
    }
  }

  if(argc != 2+hitnum){
    fprintf(stderr, "USAGE: <this> <in.fa>\n");
    fprintf(stderr, "\t-K=i: kmer's K (default: %llu)\n",K);
    fprintf(stderr, "\t-q: use fastq format\n");
    fprintf(stderr, "\t-t=float: percent of identity\n");
    fprintf(stderr, "\t-H=positive integer: H size of minhash\n");
    fprintf(stderr, "\t-T=positive integer: threshold of Jaccard similarity*H_SIZE valude of minhash\n");
    fprintf(stderr, "\t-o: output overhanging only\n");
    fprintf(stderr, "\t-m=float: specify minimum fraction per overlap against read length\n");
    fprintf(stderr, "\t-M=positive integer: discard the kmers if the number of the kmerint becomes more than or equal to this value\n");
    fprintf(stderr, "\t-P=integer: discard the radius P kmers around kmers discarded by the -M\n");
    fprintf(stderr, "\t-c: output kmer freq data only. use this for plotting\n");
    return 1;
  }

  if(opt_DEBUG){
    H_SIZE=0;
    MIN_ALI_FRAC=0.0;
    //P_SIMILAR = 1.0;
    T_PIDENT = 0.65; // XXX
  }

  char * in_fa = argv[1+hitnum];

  FILE * fp = fopen(in_fa,"r");
  if(fp == NULL){
    fprintf(stderr, "cannot open %s\n", in_fa);
    exit(1);
  }

  unsigned long long max_n_bases=800llu*(1llu<<20);
  char * allreads = (char*)malloc(sizeof(char)*max_n_bases);
  if (allreads == NULL){
    fprintf(stderr, "cannot allocate memory: allreads\n");
    exit(1);
  }
  //char * rev = (char*)malloc(sizeof(char)*max_n_bases);
  //if(rev == NULL){
    //fprintf(stderr,"cannot allocate memory: rev\n");
    //exit(1);
  //}
  //char ** reads[2];
  //reads[0] = &fwd;
  //reads[1] = &rev;

  unsigned long long max_n_records = (1llu<<20);
  unsigned long long * bounds = (unsigned long long*)malloc(sizeof(unsigned long long)*max_n_records*2);
  if(bounds == NULL){
    fprintf(stderr, "cannot allocate memory: bounds\n");
    exit(1);
  }
  unsigned long long * len_kmers = (unsigned long long*)malloc(sizeof(unsigned long long)*(max_n_records<<1));
  if(len_kmers == NULL){
    fprintf(stderr, "cannot allocate memory: len_kmers\n");
    exit(1);
  }

  char * buf = (char*)malloc(sizeof(char)*LSEQ);
  if(buf == NULL){
    fprintf(stderr, "cannot allocate memory: buf\n");
    exit(1);
  }
  name_t * names = (name_t*)malloc(sizeof(name_t)*(max_n_records<<1));
  if(names == NULL){
    fprintf(stderr, "cannot allocate memory: names\n");
    return 1;
  }
  //fprintf(stderr,"fgets 1 start\n");

  unsigned long long l_bounds=0;// l_records, l_reads
  bounds[l_bounds] = 0;
  len_kmers[l_bounds] = 0;
  allreads[0] = '\0';
  //rev[0] = '\0';
  unsigned long long l_kmers=0;
  while(fgets(buf,LSEQ,fp)!=NULL){
    chomp2(buf);
    strcpy(names[l_bounds].name,&buf[1]);
    //strcpy(names[l_bounds+1].name,&buf[1]);
    fgets(buf,LSEQ,fp);
    int line_length;
    line_length = strlen(buf);
    if(line_length >= LSEQ-1){
      while(line_length >= LSEQ-1 && buf[line_length-1] != '\n'){
        char * newbuf;
        newbuf = (char*)realloc(buf, LSEQ*2);
        if(newbuf == NULL){
          fprintf(stderr,"cannot allocate memory: newbuf\n");
          exit(1);
        }
        buf = newbuf;
        fgets(buf+line_length, LSEQ+1,fp);
        line_length = strlen(buf);
        LSEQ *= 2;
      }
    }
    chomp2(buf);
    line_length = strlen(buf);
    /*
    if(strlen(buf) < K){
      continue;
    }
    */
    if(line_length < K){
      continue;
    }
    //strcat(allreads,buf);
    strcpy(allreads+bounds[l_bounds],buf);
    ++l_bounds;
    bounds[l_bounds] = bounds[l_bounds-1]+line_length;
    len_kmers[l_bounds] = len_kmers[l_bounds-1]+(unsigned long long)(line_length-K+1);

    /*
    reversecomplementstr(buf,line_length);

    strcpy(allreads+bounds[l_bounds],buf);
    ++l_bounds;
    bounds[l_bounds] = bounds[l_bounds-1]+line_length;
    len_kmers[l_bounds] = len_kmers[l_bounds-1]+(unsigned long long)(line_length-K+1);
    */

    //bounds[l_bounds] = (unsigned long long)strlen(allreads);
    //bounds[l_bounds] = bounds[l_bounds-1]+strlen(buf);
    if(bounds[l_bounds] > max_n_bases-(LSEQ<<1)){
      fprintf(stderr, "max_n_bases exceeded. set it larger\n");
      exit(1);
    }
    //len_kmers[l_bounds] = len_kmers[l_bounds-1]+(unsigned long long)(strlen(buf)-K+1);
    if(l_bounds >= max_n_records*2){
      fprintf(stderr,"l_bounds >= max_n_records\n");
      exit(1);
    }
  }

  if(l_bounds == 0){
    fprintf(stderr, "no valid record. maybe K is too large.\n");
    exit(1);
  }

/*
  fseek(fp,0,0);
  rev_first_index = l_bounds;
  while(fgets(buf,LSEQ,fp)!=NULL){
    chomp2(buf);
    fgets(buf,LSEQ,fp);
    chomp2(buf);
    if(strlen(buf) < K){
      continue;
    }

    reversecomplement(buf);

    strcpy(names[l_bounds].name,&buf[1]);
    //strcat(allreads,buf);
    strcpy(allreads+bounds[l_bounds],buf);
    ++l_bounds;
    //bounds[l_bounds] = (unsigned long long)strlen(allreads);
    bounds[l_bounds] = bounds[l_bounds-1]+strlen(buf);
    if(bounds[l_bounds] > max_n_bases-LSEQ){
      fprintf(stderr, "max_n_bases exceeded. set it larger 2\n");
      exit(1);
    }
    len_kmers[l_bounds] = len_kmers[l_bounds-1]+(strlen(buf)-K+1);
    if(l_bounds >= max_n_records*2){
      fprintf(stderr,"l_bounds >= max_n_records\n");
      exit(1);
    }
    //reversecomplement(buf);
    //strcat(rev,buf);
  }
  */
  l_kmers=len_kmers[l_bounds];
//	printf("lb1 %llu\n",l_bounds);
  //fprintf(stderr,"fgets 2 done\n");

  MAX_CAND = LSEQ*MAX_DEPTH;// XXX

  unsigned long long kmerkinds = 1ull;
  kmerkinds <<=2*K;
  /*
  unsigned long long * htable = (unsigned long long*)malloc(sizeof(unsigned long long)*kmerkinds);
  if(htable == NULL){
    fprintf(stderr, "cannot allocate memory: htable\n");
    exit(1);
  }
  {
    unsigned long long i;
    #pragma omp parallel for
    for(i=0; i<kmerkinds; ++i){
      htable[i] = ULLONG_MAX;
    }
  }
  */
  pos_t * kmerrange = (pos_t *)malloc(sizeof(pos_t)*kmerkinds);
  if(kmerrange == NULL){
    fprintf(stderr, "cannot allocate memory: kmerrange\n");
    exit(1);
  }
  {
    unsigned long long i;
    #pragma omp parallel for
    for(i=0; i<kmerkinds; ++i){
      kmerrange[i].stt = ULLONG_MAX;
      kmerrange[i].end = ULLONG_MAX;
    }
  }

  record_t * kmers = (record_t*)malloc(sizeof(record_t)*(max_n_bases*2));
  if(kmers == NULL){
    fprintf(stderr, "cannot allocate memory: kmers\n");
    exit(1);
  }

  uint32_t * sketches = (uint32_t*)malloc(sizeof(uint32_t)*H_SIZE*max_n_records*2);
  if(sketches == NULL){
    fprintf(stderr, "cannot allocate memory: sketches\n");
    exit(1);
  }

  int NUM_THREADS = omp_get_max_threads();
  if(NUM_THREADS <= 0){
    NUM_THREADS=1;
    //fprintf(stderr, "set OMP_NUM_THREADS environment variable. othewise, this program uses only 1 thread.\n");
    fprintf(stderr, "set OMP_NUM_THREADS environment variable.\n");
    return 1;
  }
  else{
    fprintf(stderr,"omp_get_max_threads %d\n",NUM_THREADS);
    //fprintf(stderr, "NUM_THREADS: %d\n",NUM_THREADS);
  }
  int NRECORD_MAX = 400;
  bfmt7_t * b_array = (bfmt7_t*)malloc(sizeof(bfmt7_t)*NRECORD_MAX*NUM_THREADS);
  if(b_array == NULL){
    fprintf(stderr, "cannot allocate memory: b_array\n");
    return 1;
  }
  malloc_bfmt7array(b_array,NRECORD_MAX*NUM_THREADS);

  char * qseq;
  char * sseq;
  char * sbj_buf;
  char * que_buf;
  char * ali_buf;
  char * tmp_qbuf;
  char * tmp_sbuf;
  char * tmp_abuf;
  char * res_qbuf;
  char * res_sbuf;
  char * res_abuf;
  uint32_t * Gx;
  qseq = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(qseq == NULL){
    fprintf(stderr,"cannot allocate memory: qseq\n");
    return 1;
  }
  sseq = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(sseq == NULL){
    fprintf(stderr,"cannot allocate memory: sseq\n");
    return 1;
  }
  tmp_qbuf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(tmp_qbuf == NULL){
    fprintf(stderr,"cannot allocate memory: tmp_qbuf\n");
    return 1;
  }
  tmp_sbuf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(tmp_sbuf == NULL){
    fprintf(stderr,"cannot allocate memory: tmp_sbuf\n");
    return 1;
  }
  tmp_abuf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(tmp_abuf == NULL){
    fprintf(stderr,"cannot allocate memory: tmp_abuf\n");
    return 1;
  }
  res_qbuf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(res_qbuf == NULL){
    fprintf(stderr,"cannot allocate memory: res_qbuf\n");
    return 1;
  }
  res_sbuf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(res_sbuf == NULL){
    fprintf(stderr,"cannot allocate memory: res_sbuf\n");
    return 1;
  }
  res_abuf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(res_abuf == NULL){
    fprintf(stderr,"cannot allocate memory: res_abuf\n");
    return 1;
  }
  sbj_buf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(sbj_buf == NULL){
    fprintf(stderr,"cannot allocate memory: sbj_buf\n");
    return 1;
  }
  que_buf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(que_buf == NULL){
    fprintf(stderr,"cannot allocate memory: que_buf\n");
    return 1;
  }
  ali_buf = (char*)malloc(sizeof(char)*(LSEQ+1)*NUM_THREADS);
  if(ali_buf == NULL){
    fprintf(stderr,"cannot allocate memory: ali_buf\n");
    return 1;
  }
  unsigned long long ** D0ss = (unsigned long long**)malloc(sizeof(unsigned long long*)*NUM_THREADS);
  if(D0ss == NULL){
    fprintf(stderr, "cannot allocate memory: D0ss\n");
    return 1;
  }
  unsigned long long ** HPss = (unsigned long long**)malloc(sizeof(unsigned long long*)*NUM_THREADS);
  if(HPss == NULL){
    fprintf(stderr, "cannot allocate memory: HPss\n");
    return 1;
  }
  unsigned long long ** VPss = (unsigned long long**)malloc(sizeof(unsigned long long*)*NUM_THREADS);
  if(VPss == NULL){
    fprintf(stderr, "cannot allocate memory: VPss\n");
    return 1;
  }
  unsigned long long ** HNss = (unsigned long long**)malloc(sizeof(unsigned long long*)*NUM_THREADS);
  if(HNss == NULL){
    fprintf(stderr, "cannot allocate memory: HNss\n");
    return 1;
  }
  unsigned long long ** VNss = (unsigned long long**)malloc(sizeof(unsigned long long*)*NUM_THREADS);
  if(VNss == NULL){
    fprintf(stderr, "cannot allocate memory: VNss\n");
    return 1;
  }
  {
    int i;
    for(i=0; i<NUM_THREADS; ++i){
      malloc_h64(&D0ss[i],&HPss[i],&VPss[i],&HNss[i],&VNss[i]);
    }
  }
  //record_t ** p_kmers = (record_t**)malloc(sizeof(record_t*)*(l_kmers));
  record_t ** p_kmers = (record_t**)malloc(sizeof(record_t*)*max_n_bases*2);
  if(p_kmers == NULL){
    fprintf(stderr,"cannot allocate memory: p_kmers\n");
    exit(1);
  }

  unsigned long long * candidates;
  if(!opt_print_Jmatrix){
    candidates = (unsigned long long*)malloc(sizeof(unsigned long long)*MAX_CAND*NUM_THREADS);
  }
  else{
    candidates = (unsigned long long*)malloc(sizeof(unsigned long long)*LSEQ*NUM_THREADS); // candidates variable will NOT be used under this option. no need to malloc
  }
  if(candidates == NULL){
    fprintf(stderr,"cannot allocate memory: candidates\n");
    exit(1);
  }

  Gx = (uint32_t*)malloc(sizeof(uint32_t)*(LSEQ+1)*NUM_THREADS);
  if(Gx == NULL){
    fprintf(stderr,"cannot allocate memory: Gx\n");
    return 1;
  }

  int MAXFREQ = 200;
  int * freq = (int*)malloc(sizeof(int)*(MAXFREQ+1)*NUM_THREADS);
  if(freq == NULL){
    fprintf(stderr,"cannot allocate memory: freq\n");
    exit(1);
  }
  {
    int i;
    for(i=0; i<(MAXFREQ+1)*NUM_THREADS; ++i){
      freq[i] = 0;
    }
  }

  if(opt_print_Jmatrix){
    Jmat = (kv_t*)malloc(sizeof(kv_t)*l_bounds*NUM_THREADS);
    if(Jmat == NULL){
      fprintf(stderr, "cannot allocate memory: Jmat\n");
      exit(1);
    }
  }

  FILE ** out_bfmt7 = (FILE**)malloc(sizeof(FILE*)*NUM_THREADS);
  if(out_bfmt7 == NULL){
    fprintf(stderr, "cannot allocate memory: out_bfmt7\n");
    exit(1);
  }
  if(!opt_print_Jmatrix)
  {
    int i;
    char buffer[256];
    char prefix[256];
    strcpy(prefix, in_fa);
    //fprintf(stderr,"bef: %s\n",prefix);
    for(i=strlen(prefix)-1; i>0; --i){
      if(prefix[i] == '.'){
        prefix[i] = '\0';
        //fprintf(stderr,"aft: %s\n",prefix);
        //exit(1);
        break;
      }
    }
    for(i=0; i<NUM_THREADS; ++i){
      sprintf(buffer, "%s.%04d.bfmt7",prefix,i);
      out_bfmt7[i] = fopen(buffer, "w");
      if(out_bfmt7[i] == NULL){
        fprintf(stderr, "cannot open %s\n", buffer);
        exit(1);
      }
    }
  }

  {
    unsigned long long i;
    #pragma omp parallel for
    for(i=0; i<l_bounds; ++i){
      //if(i < rev_first_index){
      // TODO
      if(i & 1){
        kmerize(allreads,bounds[i],bounds[i+1],kmers+len_kmers[i],i+1,'-');
      }
      else{
        kmerize(allreads,bounds[i],bounds[i+1],kmers+len_kmers[i],i+1,'+');
      }
      //}
      //else if(i >= rev_first_index){
        //kmerize(allreads,bounds[i],bounds[i+1],&kmers[len_kmers[i]],i+1,'-');
      //}
    }
  }
  fprintf(stderr,"kmerize done\n");
  // if rev, then stt = stt+sth;
  /*
  {
    unsigned long long i;
    for(i=0; i<l_kmers; ++i){
      if(kmers[i].id & 1ull){
        kmers[i].pos += strlen(allreads);
      }
    }
  }
  */

  //record_t ** p_kmers = (record_t**)malloc(sizeof(record_t*)*(l_kmers));
  //if(p_kmers == NULL){
    //fprintf(stderr,"cannot allocate memory: p_kmers\n");
    //exit(1);
  //}
  {
    unsigned long long i;
    #pragma omp parallel for private(i)
    for(i=0; i<l_kmers; ++i){
      p_kmers[i] = &kmers[i];
    }
  }
fprintf(stderr,"qsort started. %llu\n",l_kmers);
  //qsort(p_kmers,l_kmers,sizeof(record_t*),compare_kmer_records);
  parallel_task_qsort(p_kmers,0,(signed long long)l_kmers-1);
fprintf(stderr,"qsort done\n");
  {
    unsigned long long i;
    for(i=0; i<l_kmers-1; ++i){
      if(p_kmers[i]->kmer > p_kmers[i+1]->kmer){
        fprintf(stderr, "strange %llu > %llu\n",p_kmers[i]->kmer, p_kmers[i+1]->kmer);
        exit(1);
      }
    }
    fprintf(stderr, "sort is ok\n");
    //return 0;
  }
  //printf("%s\n",allreads);
  /*
  {
    unsigned long long i;
    for(i=0; i<l_kmers; ++i){
      //record_t * p = p_kmers[i];
      //printf("%llu\t%lu\t%llu\t%c\n",p->kmer,p->id,p->offset,p->strand);
    }
  }
  */

fprintf(stderr,"kmer range started\n");
  {
     unsigned long long i,j1,j2,foo;
     unsigned long long prev;
     prev = p_kmers[0]->kmer;
     kmerrange[prev].stt = 0;
     for(i=1; i<l_kmers; ++i){
       if(p_kmers[i]->kmer != prev){
         kmerrange[prev].end = i-1;
         prev = p_kmers[i]->kmer;
         kmerrange[prev].stt = i;
       }
     }
     kmerrange[prev].end = l_kmers-1;
     #pragma omp parallel for private(i)
     for(i=0; i<kmerkinds; ++i){
       j1 = kmerrange[i].stt;
       j2 = kmerrange[i].end;
       if(j1 == j2){ // not shared
         if(j1 != ULLONG_MAX){// i is used
           p_kmers[j1]->f_valid = 0;// not shared => no need to use
         }
         else{// i is not used
         }
         //kmerrange[i].stt = ULLONG_MAX;
         //kmerrange[i].end = ULLONG_MAX;
       }
       else if(j2-j1+1 >= TOO_MANY_KMER){
         for(foo=j1; foo<=j2; ++foo){
           p_kmers[foo]->f_valid = 0;
         }
       }
     }
  }
fprintf(stderr,"kmer range ended\n");

  if(opt_count_for_gnuplot)
  {
fprintf(stderr,"kmer freq counting started\n");
     unsigned long long i;
     int thread_num;
     int count;
     #pragma omp parallel for private(i,thread_num,count)
     for(i=0; i<kmerkinds; ++i){
       thread_num = omp_get_thread_num();
       if(kmerrange[i].stt != ULLONG_MAX){// i is used
         count = kmerrange[i].end - kmerrange[i].stt + 1;
         //count /= 10;
         if(count<MAXFREQ){
           ++freq[(MAXFREQ+1)*thread_num+count];
         }
         else{
           ++freq[(MAXFREQ+1)*thread_num+MAXFREQ];
         }
       }
     }
     int i2;
     int i3;
     for(i2=0; i2<=MAXFREQ; ++i2){
       for(i3=1; i3<NUM_THREADS; ++i3){
         freq[i2] += freq[i2+(MAXFREQ+1)*i3];
       }
     }
     for(i2=0; i2<=MAXFREQ; ++i2){
       printf("%d %d\n",i2,freq[i2]);
     }
    /*
    record_t * p;
    p = p_kmers[0];
    unsigned long long prevkint = p->kmer;
    int count = 1;
    while(p != p_kmers[l_kmers-1]){
      count = 1;
      while(prevkint == p->kmer && p != p_kmers[l_kmers-1]){
        ++count;
        ++p;
      }
      //printf("%llu %d\n",prevkint,count);
      if(count<MAXFREQ){
        ++freq[count];
      }
      else{
        ++freq[MAXFREQ];
      }
      prevkint = p->kmer;
    }
    if(prevkint != p->kmer){
      //printf("%llu %d\n",prevkint,count);
      if(count<MAXFREQ){
        ++freq[count];
      }
      else{
        ++freq[MAXFREQ];
      }
    }
    int i;
    for(i=0; i<=MAXFREQ; ++i){
      printf("%d %d\n",i,freq[i]);
    }
    */
    return 0;
  }

//fprintf(stderr,"PM started\n");
  //if(PM>0)
  if(0)
  {
    // invalidate high-freq kmers and its surroundings XXX
    unsigned long long i,j,j1,j2,kint,kint2;
    unsigned long long qlen;
    //unsigned long long slen;
    unsigned long long k,k1,k2;
    unsigned long long m,m1,m2;
    #pragma omp parallel for private(i,j,j1,j2,kint,kint2,qlen,k,k1,k2,m,m1,m2)
    for(i=0; i<l_bounds; ++i){
      qlen = bounds[i+1]-bounds[i];
      for(j=0; j<=qlen-K; ++j){
        kint = kmers[len_kmers[i]+j].kmer;
        j1 = kmerrange[kint].stt;
        j2 = kmerrange[kint].end;
        if(j2-j1+1 >= TOO_MANY_KMER){
          if(j>PM){
            k1=j-PM;
          }
          else{
            k1=0;
          }
          k2=j+DM;
          if(k2>qlen-K){
            k2=qlen-K;
          }
          for(k=k1; k<=k2; ++k){
            kint2 = kmers[len_kmers[i]+k].kmer;
            m1 = kmerrange[kint2].stt;
            m2 = kmerrange[kint2].end;
            for(m=m1; m<=m2; ++m){
              p_kmers[m]->f_valid = 0;
            }
          }
        }
      }
    }
  }
//fprintf(stderr,"PM ended\n");

  // make sketches
  if(H_SIZE>0)
  {
  fprintf(stderr, "sketch started\n");
    unsigned long long i;
    unsigned long long j;
    int s;
    uint32_t minfp;
    uint32_t buf;
    unsigned long long qlen;
    unsigned long long kint;
    //unsigned long long dummystt;
    unsigned long long f_valid2;
    int thread_num;
    uint32_t tmpv;
    if(!opt_old){
    #pragma omp parallel for private(j,s,minfp,buf,qlen,kint,f_valid2,thread_num,tmpv)
    for(i=0; i<l_bounds; ++i){
      thread_num = omp_get_thread_num();
      for(s=0; s<1; ++s){
        minfp = UINT32_MAX;
        qlen = bounds[i+1]-bounds[i];
        for(j=0; j<=qlen-K; ++j){
          kint = kmers[len_kmers[i]+j].kmer;
          f_valid2 = p_kmers[kmerrange[kint].stt]->f_valid;//kmerrange[kint].stt ~ end 's f_valid are the same
          if(f_valid2){
            mymurmurhash3(kint,(uint32_t)s,&buf);
            Gx[(LSEQ+1)*thread_num+j] = buf;
            if(minfp>buf){
              minfp = buf;
            }
          }
          else{
            Gx[(LSEQ+1)*thread_num+j] = 0;// this means to discard AAAAA...
          }
        }
        sketches[H_SIZE*i+s] = minfp;
      }
      for(s=1; s<H_SIZE; ++s){
        minfp = UINT32_MAX;
        for(j=0; j<=qlen-K; ++j){
          tmpv = Gx[(LSEQ+1)*thread_num+j];
          if(tmpv != 0){
            // XORShift
            tmpv ^= tmpv << 13;
            tmpv ^= tmpv >> 17;
            tmpv ^= tmpv << 5;
            Gx[(LSEQ+1)*thread_num+j] = tmpv;
            if(minfp>tmpv){
              minfp = tmpv;
            }
          }
          else{
          }
        }
        sketches[H_SIZE*i+s] = minfp;
      }
    }
    }
    else{
    #pragma omp parallel for private(j,s,minfp,buf,qlen,kint,f_valid2,thread_num,tmpv)
    for(i=0; i<l_bounds; ++i){
      for(s=0; s<H_SIZE; ++s){
        minfp = UINT32_MAX;
        qlen = bounds[i+1]-bounds[i];
        for(j=0; j<=qlen-K; ++j){
          kint = kmers[len_kmers[i]+j].kmer;
          //dummystt = kmerrange[kint].stt;
          f_valid2 = p_kmers[kmerrange[kint].stt]->f_valid;//kmerrange[kint].stt ~ end 's f_valid are the same
          //if(dummystt != ULLONG_MAX)
          if(f_valid2){
            mymurmurhash3(kint,(uint32_t)s,&buf);
            if(minfp>buf){
              minfp = buf;
            }
          }
          else{
            // discard thie kmer
          }
        }
        sketches[H_SIZE*i+s] = minfp;
      }
    }
    }
  fprintf(stderr, "sketch done\n");
  }
  //int n_threads = omp_get_num_threads();
  //printf("n_threads: %d\n",n_threads);
  //return 0;
  //int n_threads = 2;
  //int NUM_THREADS = atoi(getenv("OMP_NUM_THREADS"));
  //NUM_THREADS=1;
//printf("nexted?: %d\n",omp_get_nested());
  fprintf(stderr, "ol started\n");
  {
    unsigned long long i;
    omp_lock_t omp_lock;
    omp_init_lock(&omp_lock);
    //#pragma omp parallel for num_threads(NUM_THREADS)
    int thread_num;
    int idx;
    unsigned long long qlen;
    unsigned long long kint;
    unsigned long long j;
    unsigned long long stt;
    unsigned long long end;
    unsigned long long k;
    unsigned long long id;
    unsigned long long slen;
    unsigned long long offset;
    int n_qleftbases;
    int n_qrightbases;
    int n_sleftbases;
    int n_srightbases;
    //int n_qleftbases2;
    //int n_qrightbases2;
    //int n_sleftbases2;
    //int n_srightbases2;
    int itmp,alilen;
    int tmpi;
    int debug_heikki_all=0;
    unsigned long long tmpid;
    int tmp;
    int f_alreadyhit = 0;
    int index;
    int b_offset;
    int tmpalilen;
    int alii;
    int match;
    int l_offset;
    int foo;
    int l_can;
    int tmpthre_j;
    int f_discard_this_record;
    int min_ali_score_for_this;
    int tmp_jsl;

    if(opt_print_Ja){
      int val;
      #pragma omp parallel for private (val)
      for(i=0; i<l_bounds; ++i){
        for(j=0; j<l_bounds; ++j){
          if(i==j){
            continue;
          }
          val = get_Jaccard_similarity_like(i,j,sketches);
          if(val){
            printf("%d\n",val);
          }
        }
      }
      return 0;
    }
    if(opt_print_Jmatrix){
      thread_num = omp_get_thread_num();
      //#pragma omp parallel for private (thread_num,j)
      for(i=0; i<l_bounds; ++i){
        for(j=0; j<l_bounds; ++j){
          Jmat[l_bounds*thread_num+j].key = (int)j;
          Jmat[l_bounds*thread_num+j].value = get_Jaccard_similarity_like(i,j,sketches);
        }
        qsort(Jmat+l_bounds*thread_num, l_bounds, sizeof(kv_t), compare_kv);
        if(!(i & 1)){
          omp_set_lock(&omp_lock);
          for(j=0; j<l_bounds-1; ++j){
            printf("%llu,%d,%d\t",i+1,Jmat[l_bounds*thread_num+j].key+1,Jmat[l_bounds*thread_num+j].value);// ids are 1-based
          }
          printf("%llu,%d,%d\n",i+1,Jmat[l_bounds*thread_num+j].key+1,Jmat[l_bounds*thread_num+j].value);
          omp_unset_lock(&omp_lock);
        }
      }
      return 0;
    }
    #pragma omp parallel for private (thread_num,idx,qlen,kint,j,stt,end,k,id,slen,offset,n_qleftbases,n_qrightbases,n_sleftbases,n_srightbases,itmp,alilen,tmpi,debug_heikki_all,tmpid,tmp,f_alreadyhit,index,b_offset,tmpalilen,alii,match,l_offset,foo,l_can,tmpthre_j,f_discard_this_record,min_ali_score_for_this,tmp_jsl)
    for(i=0; i<l_bounds; i+=2){// fwd only
      thread_num = omp_get_thread_num();
      //int thread_num = 0;
      //printf("thread_num: %d\n",thread_num);
      idx = (int)((LSEQ+1)*thread_num);

      //fprintf(stderr,"kita0 %llu\n",i);


      l_can = 0;
      f_discard_this_record=0;
      if(H_SIZE>0){
        //l_can = 0;
    //fprintf(stderr,"before get_Js %llu\n",i);
        tmpthre_j = THRE_J;
        for(j=0; j<l_bounds; ++j){
          if(i==j){
            continue;
          }
          tmp_jsl = get_Jaccard_similarity_like(i,j,sketches);
          if(tmp_jsl >= tmpthre_j && tmp_jsl <= THRE_J_UPPER){
            candidates[MAX_CAND*thread_num+l_can] = j;
            if(l_can <= L_CAN_MAX && l_can < MAX_CAND){
              ++l_can;
            }
            else{
              // TODO
              break;
            }
          }
        }
        if(l_can > L_CAN_MAX){
          f_discard_this_record = 1;
        }
        //fprintf(stderr,"%d\n",l_can);
    //fprintf(stderr,"after get_Js %llu\n",i);
      }
      if(f_discard_this_record){
        continue;
      }

      //printf("init before heikki: idx: %d\n",idx);
      //int thread_num = 0;
      sseq[idx] = '\0';
      qseq[idx] = '\0';
      b_array[(NRECORD_MAX)*thread_num].l_data = 0;
      //unsigned long long MAXHIT = 1;
      qlen = bounds[i+1]-bounds[i];
      if(qlen>0){
        strncpy(&qseq[idx],&allreads[bounds[i]],qlen);
        qseq[idx+qlen] = '\0';
      }
      else{
        fprintf(stderr,"bug? qlen=%llu\n",qlen);
      }

      //min_ali_score_for_this = qlen * MIN_ALI_FRAC;

      for(j=0; j<=qlen-K; ++j){
        kint = kmers[len_kmers[i]+j].kmer;
        stt = kmerrange[kint].stt;
        end = kmerrange[kint].end;
        //if(stt == ULLONG_MAX)// not shared

        if(p_kmers[stt]->f_valid==0){
          //printf("kita ULLONG_MAX\n");
          continue;
        }

        for(k=stt; k<=end && k<stt+NRECORD_MAX; ++k){
          id = p_kmers[k]->id;
          if(i+1 == id){// self hit
            continue;
          }
          if(i+2 == id){// rev.c hist
            continue;
          }
          //unsigned long long slen = bounds[(id-1)+1]-bounds[id-1];
          if(H_SIZE>0){
            if(!f_hit(id-1,candidates+MAX_CAND*thread_num,l_can)){
              //fprintf(stderr,"kita f_hit\n");
              continue;
            }
          }

          index = NRECORD_MAX*thread_num;
          b_offset = b_array[index].l_data;
          if(b_offset >= NRECORD_MAX){
            //fprintf(stderr,"koreka b_offset: %d\n",b_offset);
            //discard XXX
            continue;
          }

          f_alreadyhit = 0;
          for(tmp=0; tmp<b_array[NRECORD_MAX*thread_num].l_data; ++tmp){
            tmpid = id;
            if(b_array[NRECORD_MAX*thread_num+tmp].sid == tmpid){
              f_alreadyhit = 1;
              break;
            }
          }
          if(f_alreadyhit){
            // I don't do my best.
            continue;
          }

          slen = bounds[id]-bounds[id-1];
          if(slen>0){
            strncpy(sseq+idx,allreads+bounds[id-1],slen);
            sseq[idx+slen] = '\0';
          }
          else{
            fprintf(stderr,"bug? slen=%llu\n",slen);
          }

          // adjust min_ali_score_for_this to the shorter read
          if(qlen <= slen){
            min_ali_score_for_this = qlen * MIN_ALI_FRAC;
          }
          else{
            min_ali_score_for_this = slen * MIN_ALI_FRAC;
          }

          if(j>=qlen-min_ali_score_for_this){
            // I don't think this search is promising.
            continue;
          }

          //printf("qseq-1: %s thread_num: %d, stt: %llu, end: %llu\n",qseq+idx,thread_num,stt,end);
          //printf("sseq-1: %s thread_num: %d, stt: %llu, end: %llu\n",sseq+idx,thread_num,stt,end);
          offset = p_kmers[k]->offset;
          if(offset >= slen-min_ali_score_for_this){
            // I don't think this search is promising too.
            continue;
          }
          res_sbuf[idx] = '\0';
          res_qbuf[idx] = '\0';
          res_abuf[idx] = '\0';
          que_buf[idx] = '\0';
          ali_buf[idx] = '\0';
          sbj_buf[idx] = '\0';
          n_qleftbases = 0;
          n_qrightbases = 0;
          n_sleftbases = 0;
          n_srightbases = 0;

          // against hash collision + k-mer extension
          strncpy(tmp_qbuf+idx,qseq+idx+j,K+K_SUGAR);
          tmp_qbuf[idx+K+K_SUGAR] = '\0';
          strncpy(tmp_sbuf+idx,sseq+idx+offset,K+K_SUGAR);
          tmp_sbuf[idx+K+K_SUGAR] = '\0';
          f_discard_this_record = 0;
          for(itmp=K+K_SUGAR-1; itmp>=0 ; --itmp){
            if(tmp_qbuf[idx+itmp] != tmp_sbuf[idx+itmp]){
              f_discard_this_record = 1;
              break;
            }
          }
          if(f_discard_this_record){
            continue;
          }

          //printf("check: thread_num: %d, idx: %d\n",thread_num,idx);
          // left
          //printf("thread_num: %d, offset: %llu, j: %llu,qlen-K: %llu\n",thread_num,offset,j,qlen-K);
          if(offset>=1 && j>=1){
            //printf("kita\n");
            strncpy(tmp_qbuf+idx,qseq+idx,j);
            tmp_qbuf[idx+j] = '\0';
            strncpy(tmp_sbuf+idx,sseq+idx,offset);
            tmp_sbuf[idx+offset] = '\0';
            //printf("left q: %s\tj     : %llu\n",tmp_qbuf+idx,j);
            //printf("left s: %s\toffset: %llu\n",tmp_sbuf+idx,offset);
            reversestr(tmp_sbuf+idx,offset);
            reversestr(tmp_qbuf+idx,j);
          res_sbuf[idx+0] = '\0';
          res_qbuf[idx+0] = '\0';
          res_abuf[idx+0] = '\0';


//fprintf(stderr,"bef s: %s\n",tmp_sbuf+idx);
//fprintf(stderr,"bef q: %s\n",tmp_qbuf+idx);
          //fprintf(stderr, "%s\n",tmp_qbuf+idx);
          //fprintf(stderr, "%s\n",tmp_sbuf+idx);
            //printf("before heikki len res_sbuf: %d\n",(int)strlen(res_sbuf+idx));
            //printf("before heikki len res_qbuf: %d\n",(int)strlen(res_qbuf+idx));
            //printf("before heikki len res_abuf: %d\n",(int)strlen(res_abuf+idx));
            //printf("before heikki: res-q: %s que_buf: %s res-q-len: %d que_buf_len: %d idx: %d\n",res_qbuf+idx,que_buf+idx,(int)strlen(res_qbuf+idx),(int)strlen(que_buf+idx),idx);
            //printf("before heikki: res-a: %s ali_buf: %s res-a-len: %d ali_buf_len: %d idx: %d\n",res_abuf+idx,ali_buf+idx,(int)strlen(res_abuf+idx),(int)strlen(ali_buf+idx),idx);
            //printf("before heikki: res-s: %s sbf_buf: %s res-s-len: %d sbf_buf_len: %d idx: %d\n",res_sbuf+idx,sbj_buf+idx,(int)strlen(res_sbuf+idx),(int)strlen(sbj_buf+idx),idx);
          //fprintf(stdout,"##left stt q: %llu, s: %llu\n",i+1,id);
            b_heikki_re(tmp_sbuf+idx, tmp_qbuf+idx, sbj_buf+idx, que_buf+idx, ali_buf+idx,D0ss[thread_num],HPss[thread_num],VPss[thread_num],HNss[thread_num],VNss[thread_num],LSEQ,1.0-T_PIDENT,L_CHECK,L_CHUNK);
          //fprintf(stdout,"##left end q: %llu, s: %llu\n",i+1,id);
          //fprintf(stderr,"# b_heikki done\n");
            if(que_buf[idx] == '\0'){// out of the band
          //    fprintf(stdout,"left out of the band %llu, %llu\n",i+1,id);
              continue;
            }
            //checker('s', tmp_sbuf+idx, sbj_buf+idx);
            //checker('q', tmp_qbuf+idx, que_buf+idx);
            {
              int tmpalilen2 = strlen(ali_buf+idx);
              int match2 = 0;
              int i2;
              for(i2=0; i2<tmpalilen2; ++i2){
                if(ali_buf[idx+i2] == '|'){
                  ++match2;
                }
              }
              double pident2 = (double)match2/(double)tmpalilen2;
              if(pident2 < T_PIDENT){
              //fprintf(stderr,"kita T_PIDNET %llu, %llu\n",i+1,id);
              //fprintf(stderr,"%s\n",que_buf+idx);
              //fprintf(stderr,"%s\n",ali_buf+idx);
              //fprintf(stderr,"%s\n",sbj_buf+idx);
                continue;
              }
              else{
                //fprintf(stdout,"T_PIDENT ok %f\n",pident2);
              }
            }
            //printf("before strcpy: res-q: %s que_buf: %s len: %d\n",res_qbuf+idx,que_buf+idx,(int)strlen(que_buf+idx));
            //printf("before strcpy: res-a: %s ali_buf: %s len: %d\n",res_abuf+idx,ali_buf+idx,(int)strlen(ali_buf+idx));
            //printf("before strcpy: res-s: %s sbf_buf: %s len: %d\n",res_sbuf+idx,sbj_buf+idx,(int)strlen(sbj_buf+idx));
            strcpy(res_qbuf+idx,que_buf+idx);
            strcpy(res_abuf+idx,ali_buf+idx);
            strcpy(res_sbuf+idx,sbj_buf+idx);
//fprintf(stderr,"ori %s\n",res_abuf+idx);
            //printf("thread_num: %d\n",thread_num);
            //printf("strcpy: res-q: %s que_buf: %s len: %d\n",res_qbuf+idx,que_buf+idx,(int)strlen(que_buf+idx));
            //printf("strcpy: res-a: %s ali_buf: %s len: %d\n",res_abuf+idx,ali_buf+idx,(int)strlen(ali_buf+idx));
            //printf("strcpy: res-s: %s sbf_buf: %s len: %d\n",res_sbuf+idx,sbj_buf+idx,(int)strlen(sbj_buf+idx));
            //my_strcpy(res_sbuf+idx,sbj_buf+idx,(int)strlen(sbj_buf+idx));
            //my_strcpy(res_qbuf+idx,que_buf+idx,(int)strlen(que_buf+idx));
            //my_strcpy(res_abuf+idx,ali_buf+idx,(int)strlen(ali_buf+idx));
            //printf("mystrcpy: res-q: %s que_buf: %s len: %d\n",res_qbuf+idx,que_buf+idx,(int)strlen(que_buf+idx));
            //printf("mystrcpy: res-a: %s ali_buf: %s len: %d\n",res_abuf+idx,ali_buf+idx,(int)strlen(ali_buf+idx));
            //printf("mystrcpy: res-s: %s sbf_buf: %s len: %d\n",res_sbuf+idx,sbj_buf+idx,(int)strlen(sbj_buf+idx));
            reversestr(res_sbuf+idx,strlen(res_sbuf+idx));
            reversestr(res_qbuf+idx,strlen(res_qbuf+idx));
            reversestr(res_abuf+idx,strlen(res_abuf+idx));
//fprintf(stderr,"rev %s\n",res_abuf+idx);
            alilen = strlen(res_abuf+idx);
            //int qbuflen = strlen(res_qbuf+idx);
            //if(alilen != qbuflen){
            //  fprintf(stderr,"kitaaaa\n");
            //fprintf(stderr, "revstr: res-q: %s que_buf: %s len: %d\n",res_qbuf+idx,que_buf+idx,(int)strlen(que_buf+idx));
            //fprintf(stderr, "revstr: res-a: %s ali_buf: %s len: %d\n",res_abuf+idx,ali_buf+idx,(int)strlen(ali_buf+idx));
            //fprintf(stderr, "revstr: res-s: %s sbf_buf: %s len: %d\n",res_sbuf+idx,sbj_buf+idx,(int)strlen(sbj_buf+idx));
            //  exit(1);
            //}
            for(itmp=0;itmp<alilen;++itmp){
              if(res_sbuf[idx+itmp] != '-'){
                ++n_sleftbases;
                if(n_sleftbases > offset){
                  //printf("#### n_sleftbases\n");
                  //printf("qid: %llu\n",i+1);
                  //printf("sid: %llu\n",id);
                  //printf("j(kmer pos): %llu\n",j);
                  //printf("qseq: %s\n",qseq+idx);
                  //printf("sseq: %s\n",sseq+idx);
                  //printf("left-qseq: %s\n",tmp_qbuf+idx);
                  //printf("left-sseq: %s\n",tmp_sbuf+idx);
                  //printf("res-q: %s\n",res_qbuf+idx);
                  //printf("res-a: %s\n",res_abuf+idx);
                  //printf("res-s: %s\n",res_sbuf+idx);
                  //printf("que_buf: %s\n",que_buf+idx);
                  //printf("ali_buf: %s\n",ali_buf+idx);
                  //printf("sbj_buf: %s\n",sbj_buf+idx);
                  //printf("tmpq: %s\n",tmp_qbuf+idx);
                  //printf("tmps: %s\n",tmp_sbuf+idx);
                  //printf("offset: %llu\n",offset);
                  //printf("n_sleftbases: %d\n",n_sleftbases);
                  //printf("####\n");
                }
              }
              if(res_qbuf[idx+itmp] != '-'){
                ++n_qleftbases;
                if(n_qleftbases > j){
                  //printf("#### n_qleftbases\n");
                  //printf("qid: %llu\n",i+1);
                  //printf("sid: %llu\n",id);
                  //printf("j(kmer pos): %llu\n",j);
                  //printf("qseq: %s\n",qseq+idx);
                  //printf("sseq: %s\n",sseq+idx);
                  //printf("q: %s\n",res_qbuf+idx);
                  //printf("a: %s\n",res_abuf+idx);
                  //printf("s: %s\n",res_sbuf+idx);
                  //printf("tmpq: %s\n",tmp_qbuf+idx);
                  //printf("tmps: %s\n",tmp_sbuf+idx);
                  //printf("####\n");
                }
              }
            }
            //if((int)strlen(tmp_sbuf+idx) != n_sleftbases){
            //  fprintf(stderr, "diff1 %d, %d\n",(int)strlen(tmp_sbuf+idx),n_sleftbases);
            //  exit(1);
            //}
            //if((int)strlen(tmp_qbuf+idx) != n_qleftbases){
            //  fprintf(stderr, "diff2 %d, %d\n",(int)strlen(tmp_qbuf+idx),n_qleftbases);
            //  exit(1);
            //}
          }
            /*
            printf("l after  qseq: %s\n",qseq+idx);
            printf("l after  sseq: %s\n",sseq+idx);
            printf("l q: %s\n",res_qbuf+idx);
            printf("l a: %s\n",res_abuf+idx);
            printf("l s: %s\n",res_sbuf+idx);
            printf("####l \n");
            */

          // kmer
          tmp_qbuf[idx] = '\0';

          strncpy(tmp_qbuf+idx,qseq+idx+j,K);
          tmp_qbuf[idx+K] = '\0';

          //fprintf(stderr,"kmer: %s\n",tmp_qbuf+idx);
//strcpy(b_array[NRECORD_MAX*thread_num+b_array[NRECORD_MAX*thread_num].l_data].kmerseq,tmp_qbuf+idx);

          strcat(res_sbuf+idx,tmp_qbuf+idx);
          strcat(res_qbuf+idx,tmp_qbuf+idx);

          for(tmpi=0; tmpi<K; ++tmpi){
            tmp_qbuf[idx+tmpi] = '|';
          }
          tmp_qbuf[idx+tmpi] = '\0';
          strcat(res_abuf+idx,tmp_qbuf+idx);

          // right
          if(offset+K <= slen-1 && j+K <= qlen-1){
            strncpy(tmp_sbuf+idx,sseq+idx+offset+K,slen-offset-K);
            tmp_sbuf[idx+slen-offset-K] = '\0';
            strncpy(tmp_qbuf+idx,qseq+idx+j+K,qlen-j-K);
            tmp_qbuf[idx+qlen-j-K] = '\0';
          //  fprintf(stdout,"##right stt q: %llu, s: %llu\n",i+1,id);
//strncpy(b_array[NRECORD_MAX*thread_num+b_array[NRECORD_MAX*thread_num].l_data].bbuf,tmp_qbuf+idx,K);
            b_heikki_re(tmp_sbuf+idx, tmp_qbuf+idx, sbj_buf+idx, que_buf+idx, ali_buf+idx,D0ss[thread_num],HPss[thread_num],VPss[thread_num],HNss[thread_num],VNss[thread_num],LSEQ,1.0-T_PIDENT,L_CHECK,L_CHUNK);
//strncpy(b_array[NRECORD_MAX*thread_num+b_array[NRECORD_MAX*thread_num].l_data].abuf,que_buf+idx,K);
//strncpy(b_array[NRECORD_MAX*thread_num+b_array[NRECORD_MAX*thread_num].l_data].alibuf,ali_buf+idx,K);
//strncpy(b_array[NRECORD_MAX*thread_num+b_array[NRECORD_MAX*thread_num].l_data].subbuf,sbj_buf+idx,K);
//if(strlen(tmp_sbuf+idx)<strlen(tmp_qbuf+idx)){
//  strcat(b_array[NRECORD_MAX*thread_num+b_array[NRECORD_MAX*thread_num].l_data].subbuf,"sw");
//}
//else{
//  strcat(b_array[NRECORD_MAX*thread_num+b_array[NRECORD_MAX*thread_num].l_data].subbuf,"ns");
//}
          //  fprintf(stdout,"##right end q: %llu, s: %llu\n",i+1,id);
            if(que_buf[idx] == '\0'){// out of the band
              //fprintf(stdout,"right out of the band %llu, %llu\n",i+1,id);
              //fprintf(stdout,"%s\n",tmp_qbuf+idx);
              //fprintf(stdout,"%s\n",tmp_sbuf+idx);
              continue;
            }
            {
              int tmpalilen2 = strlen(ali_buf+idx);
              int match2 = 0;
              int i2;
              for(i2=0; i2<tmpalilen2; ++i2){
                if(ali_buf[idx+i2] == '|'){
                  ++match2;
                }
              }
              double pident2 = (double)match2/(double)tmpalilen2;
              if(pident2 < T_PIDENT){
              //fprintf(stderr,"kita2 T_PIDNET %llu, %llu\n",i+1,id);
              //fprintf(stderr,"%s\n",que_buf+idx);
              //fprintf(stderr,"%s\n",ali_buf+idx);
              //fprintf(stderr,"%s\n",sbj_buf+idx);
                continue;
                //fprintf(stdout,"T_PIDENT2 bad %f\n",pident2);
              }
              else{
                //fprintf(stdout,"T_PIDENT2 ok %f\n",pident2);
              }
            }
            strcat(res_sbuf+idx,sbj_buf+idx);
            strcat(res_qbuf+idx,que_buf+idx);
            strcat(res_abuf+idx,ali_buf+idx);
            alilen = strlen(ali_buf+idx);
            for(itmp=0;itmp<alilen;++itmp){
              if(sbj_buf[idx+itmp] != '-'){
                ++n_srightbases;
              }
              if(que_buf[idx+itmp] != '-'){
                ++n_qrightbases;
              }
            }
          }
// TODO
//fprintf(stderr,"done %s\n",res_abuf+idx);

      //fprintf(stderr,"kita1.5 j: %llu\n",j);
          //if(strlen(qseq) < strlen(res_sbuf))
          //if(strlen(qseq) < n_qleftbases+K+n_qrightbases){
          //  printf("qid: %llu\n",i+1);
          //  printf("sid: %llu\n",id);
          //  printf("j(kmer pos): %llu\n",j);
          //  printf("qseq: %s\n",qseq+idx);
          //  printf("sseq: %s\n",sseq+idx);
          //  printf("q: %s\n",res_qbuf+idx);
          //  printf("a: %s\n",res_abuf+idx);
          //  printf("s: %s\n",res_sbuf+idx);
          //  printf("tmpq: %s\n",tmp_qbuf+idx);
          //  printf("tmps: %s\n",tmp_sbuf+idx);
          //  printf("####\n");
          //}

          //b_heikki_re(&sseq[(LSEQ+1)*thread_num], &qseq[(LSEQ+1)*thread_num], &sbj_buf[(LSEQ+1)*thread_num], &que_buf[(LSEQ+1)*thread_num], &ali_buf[(LSEQ+1)*thread_num],D0ss[thread_num],HPss[thread_num],VPss[thread_num],HNss[thread_num],VNss[thread_num]);
          debug_heikki_all=0;
          if(debug_heikki_all){
            printf("before qseq: %s\n",qseq+idx);
            printf("before sseq: %s\n",sseq+idx);
            b_heikki_re(sseq+idx, qseq+idx, sbj_buf+idx, que_buf+idx, ali_buf+idx,D0ss[thread_num],HPss[thread_num],VPss[thread_num],HNss[thread_num],VNss[thread_num],LSEQ,1.0-T_PIDENT,L_CHECK,L_CHUNK);
            printf("after  qseq: %s\n",qseq+idx);
            printf("after  sseq: %s\n",sseq+idx);
            printf("q: %s\n",que_buf+idx);
            printf("a: %s\n",ali_buf+idx);
            printf("s: %s\n",sbj_buf+idx);
            printf("####\n");
          }
          //printf("b_heikki_re done\n");
          // move data to b_array
          /*
          f_alreadyhit = 0;
          for(tmp=0; tmp<b_array[NRECORD_MAX*thread_num].l_data; ++tmp){
            tmpid = id;
            if(tmpid-1 >= rev_first_index){
              tmpid -= rev_first_index;
            }
            if(b_array[NRECORD_MAX*thread_num+tmp].sid == tmpid){
              f_alreadyhit = 1;
              break;
            }
          }
          */
          //printf("flag done. flag:%d\n",f_alreadyhit);
          //if(n_qleftbases+K+n_qrightbases < MIN_OL){
            //continue;
          //}
          //n_qleftbases2 = 0;
          //n_qrightbases2 = 0;
          //n_sleftbases2 = 0;
          //n_srightbases2 = 0;
          if(!f_alreadyhit){
            index = NRECORD_MAX*thread_num;
            b_offset = b_array[index].l_data;
            if(b_offset >= NRECORD_MAX){
              //fprintf(stderr,"koreka b_offset: %d\n",b_offset);
              //discard XXX
              continue;
            }
            if(j > offset+EPS && slen-offset > qlen-j+EPS){
              // right
              //q ==>
              //s  ==>
              //fprintf(stderr,"kita1 stt:%llu\n",j-n_qleftbases+1);
              b_array[index+b_offset].qstt = j-n_qleftbases+1;//1-based
              b_array[index+b_offset].qend = qlen;
              b_array[index+b_offset].sstt = 1;
              b_array[index+b_offset].send = offset+K+n_srightbases;
            }
            else if(offset > j+EPS && qlen-j > slen-offset+EPS){
              // left
              //q  ==>
              //s ==>
              //printf("kita2\n");
              //if(offset-n_sleftbases == 0){
                //fprintf(stderr,"kita\n");
              //}
              b_array[index+b_offset].qstt = 1;
              b_array[index+b_offset].qend = j+K+n_qrightbases;
              b_array[index+b_offset].sstt = offset-n_sleftbases+1;
              b_array[index+b_offset].send = slen;
            }
            else if(j > offset+EPS && qlen-j > slen-offset+EPS){
              if(opt_overhangingonly){
                continue;
              }
              // containing
              //q =====>
              //s   =>
              b_array[index+b_offset].qstt = j-n_qleftbases+1;
              b_array[index+b_offset].qend = j+K+n_qrightbases;
              //printf("kita3\n");
              b_array[index+b_offset].sstt = 1;
              b_array[index+b_offset].send = slen;
            }
            else if(j == offset && qlen-j == slen-offset){
              if(opt_overhangingonly){
                continue;
              }
              // equal
              //q =====>
              //s =====>
              b_array[index+b_offset].qstt = 1;
              b_array[index+b_offset].qend = qlen;
              //printf("kita3\n");
              b_array[index+b_offset].sstt = 1;
              b_array[index+b_offset].send = slen;
            }
            else if(slen >= qlen && ((double)(slen-qlen)/(double)qlen) < 1.0-P_SIMILAR){ // pident filter was already done
              if(opt_overhangingonly){
                continue;
              }
              // contained but similar
              //q  ===>
              //s =====>
              b_array[index+b_offset].qstt = 1;
              b_array[index+b_offset].qend = qlen;
              b_array[index+b_offset].sstt = offset-n_sleftbases+1;
              b_array[index+b_offset].send = offset+K+n_srightbases;
            }
            else{
              // contained, self or similar
              // discard
              if(opt_DEBUG){
                // force to output
                b_array[index+b_offset].qstt = j-n_qleftbases+1;
                b_array[index+b_offset].qend = j+K+n_qrightbases;
                b_array[index+b_offset].sstt = offset-n_sleftbases+1;
                b_array[index+b_offset].send = offset+K+n_srightbases;
              }
              else{
//fprintf(stderr,"kita\n");
                continue;
              }
            }

            b_array[index+b_offset].qid = i+1;// +1: index to id
            b_array[index+b_offset].sid = id;
            strcpy(b_array[index+b_offset].qseq, res_qbuf+idx); 
            strcpy(b_array[index+b_offset].sseq, res_sbuf+idx); 
            strcpy(b_array[index+b_offset].ali, res_abuf+idx); 
            tmpalilen = strlen(b_array[index+b_offset].ali);
            match = 0;
            for(alii=0; alii<tmpalilen; ++alii){
              if(b_array[index+b_offset].ali[alii] == '|'){
                ++match;
              }
            }
      //fprintf(stderr,"kita1.9 j: %llu\n",j);
            b_array[index+b_offset].pident = (double)match/(double)tmpalilen;
            //b_array[index+b_offset].aliscore = 100.0*b_array[index+b_offset].pident;//TODO
            b_array[index+b_offset].aliscore = match;//TODO
            b_array[index+b_offset].evalue = 0.0;//TODO
            b_array[index+b_offset].qlen = (int)qlen;
            b_array[index+b_offset].slen = (int)slen;
//fprintf(stderr, "prev %llu, %llu\n",i+1,id);
            if(b_array[index+b_offset].aliscore < min_ali_score_for_this){// TODO
//fprintf(stderr, "kita min_ali %llu, %llu\n",i+1,id);
              continue;
            }
//fprintf(stderr, "pident %llu, %llu, %f (%f)\n",i+1,id,b_array[index+b_offset].pident,T_PIDENT);
//fprintf(stderr, "%s\n",b_array[index+b_offset].ali);
//fprintf(stderr, "%s\n",b_array[index+b_offset].qseq);
//fprintf(stderr, "%s\n",b_array[index+b_offset].sseq);
            if(b_array[index+b_offset].pident < T_PIDENT){
              continue;
            }
            else{
              b_array[index].l_data += 1;
//fprintf(stderr, "hit %llu, %llu\n",i+1,id);
            }
          }
          else{
            // do nothing
//fprintf(stderr,"kita already hit\n");
          }
      //fprintf(stderr,"kita2 %llu\n",k);
        }
      }
      //printf("i=%llu loop done\n",i);
      index = NRECORD_MAX*thread_num;
      l_offset = b_array[index].l_data;
      // output in bfmt7
    if(opt_single){
      omp_set_lock(&omp_lock);
      //printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%.2lf\t%.2lf\t%.2lf\t%s\t%s\n",
      //  i+1,1,(int)qlen,i+1,1,(int)qlen,100.0,0.0,1.0,&qseq[(LSEQ+1)*thread_num],&qseq[(LSEQ+1)*thread_num]);
      //printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\n",
        //i+1,1,(int)qlen,i+1,1,(int)qlen,(int)qlen,0.0,1.0,&qseq[(LSEQ+1)*thread_num],&qseq[(LSEQ+1)*thread_num],(int)qlen,(int)qlen);
      //printf("%s\t%d\t%d\t%s\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\n",
      //  names[i].name,1,(int)qlen,names[i].name,1,(int)qlen,(int)qlen,0.0,1.0,&qseq[(LSEQ+1)*thread_num],&qseq[(LSEQ+1)*thread_num],(int)qlen,(int)qlen);
      printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t++\n",
        i+1,1,(int)qlen,i+1,1,(int)qlen,(int)qlen,0.0,1.0,&qseq[(LSEQ+1)*thread_num],&qseq[(LSEQ+1)*thread_num],(int)qlen,(int)qlen);
      for(foo=0; foo<l_offset; ++foo){
        DEBUG_ALI_PRINT = 0;
        if(DEBUG_ALI_PRINT){
          printf("qid: %llu\n",b_array[index+foo].qid);
          printf("sid: %llu\n",b_array[index+foo].sid);
          printf("%s\n",b_array[index+foo].qseq);
          printf("%s\n",b_array[index+foo].ali);
          printf("%s\n",b_array[index+foo].sseq);
          printf("####\toffset=%d\n",foo);
        }
        else{
          if((b_array[index+foo].sid-1) & 1){// sseq is rev.c
            // +,-
            printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t+-\n",
              b_array[index+foo].qid,
              //names[b_array[index+foo].qid-1].name,
              b_array[index+foo].qstt,
              b_array[index+foo].qend,
              b_array[index+foo].sid,
              //names[b_array[index+foo].sid-1].name,
              b_array[index+foo].sstt,
              b_array[index+foo].send,
              b_array[index+foo].aliscore,
              b_array[index+foo].evalue,
              b_array[index+foo].pident,
              b_array[index+foo].qseq,
              b_array[index+foo].sseq,
              b_array[index+foo].qlen,
              b_array[index+foo].slen
            );
          }
          else{
            // +,+
            printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t++\n",
              b_array[index+foo].qid,
              //names[b_array[index+foo].qid-1].name,
              b_array[index+foo].qstt,
              b_array[index+foo].qend,
              b_array[index+foo].sid,
              //names[b_array[index+foo].sid-1].name,
              b_array[index+foo].sstt,
              b_array[index+foo].send,
              b_array[index+foo].aliscore,
              b_array[index+foo].evalue,
              b_array[index+foo].pident,
              b_array[index+foo].qseq,
              b_array[index+foo].sseq,
              b_array[index+foo].qlen,
              b_array[index+foo].slen
            );
          }
        }
      }
      reversecomplementstr(&qseq[(LSEQ+1)*thread_num],qlen);
      printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t--\n",
        i+2,1,(int)qlen,i+2,1,(int)qlen,(int)qlen,0.0,1.0,&qseq[(LSEQ+1)*thread_num],&qseq[(LSEQ+1)*thread_num],(int)qlen,(int)qlen);
      for(foo=0; foo<l_offset; ++foo){
        DEBUG_ALI_PRINT = 0;
        if(DEBUG_ALI_PRINT){
          printf("qid: %llu\n",b_array[index+foo].qid);
          printf("sid: %llu\n",b_array[index+foo].sid);
          printf("%s\n",b_array[index+foo].qseq);
          printf("%s\n",b_array[index+foo].ali);
          printf("%s\n",b_array[index+foo].sseq);
          printf("####\toffset=%d\n",foo);
        }
        else{
          if((b_array[index+foo].sid-1) & 1){// sseq is rev.c
            // +,- -> -,+
            reversecomplementstr(b_array[index+foo].qseq,strlen(b_array[index+foo].qseq));
            reversecomplementstr(b_array[index+foo].sseq,strlen(b_array[index+foo].sseq));
            printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t-+\n",
              b_array[index+foo].qid+1,
              //names[b_array[index+foo].qid-1].name,
              b_array[index+foo].qlen - b_array[index+foo].qend + 1,
              b_array[index+foo].qlen - b_array[index+foo].qstt + 1,
              b_array[index+foo].sid-1,
              //names[b_array[index+foo].sid-1].name,
              b_array[index+foo].slen - b_array[index+foo].send + 1,
              b_array[index+foo].slen - b_array[index+foo].sstt + 1,
              b_array[index+foo].aliscore,
              b_array[index+foo].evalue,
              b_array[index+foo].pident,
              b_array[index+foo].qseq,
              b_array[index+foo].sseq,
              b_array[index+foo].qlen,
              b_array[index+foo].slen
            );
          }
          else{
            // +,+ -> -,-
            reversecomplementstr(b_array[index+foo].qseq,strlen(b_array[index+foo].qseq));
            reversecomplementstr(b_array[index+foo].sseq,strlen(b_array[index+foo].sseq));
            printf("%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t--\n",
              b_array[index+foo].qid+1,
              //names[b_array[index+foo].qid-1].name,
              b_array[index+foo].qlen - b_array[index+foo].qend + 1,
              b_array[index+foo].qlen - b_array[index+foo].qstt + 1,
              b_array[index+foo].sid+1,
              //names[b_array[index+foo].sid-1].name,
              b_array[index+foo].slen - b_array[index+foo].send + 1,
              b_array[index+foo].slen - b_array[index+foo].sstt + 1,
              b_array[index+foo].aliscore,
              b_array[index+foo].evalue,
              b_array[index+foo].pident,
              b_array[index+foo].qseq,
              b_array[index+foo].sseq,
              b_array[index+foo].qlen,
              b_array[index+foo].slen
            );
          }
        }
      }
      //printf("i=%llu\n",i);
      omp_unset_lock(&omp_lock);
    }
    else{
      fprintf(out_bfmt7[thread_num],"%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t++\n",
        i+1,1,(int)qlen,i+1,1,(int)qlen,(int)qlen,0.0,1.0,&qseq[(LSEQ+1)*thread_num],&qseq[(LSEQ+1)*thread_num],(int)qlen,(int)qlen);
      for(foo=0; foo<l_offset; ++foo){
        if((b_array[index+foo].sid-1) & 1){// sseq is rev.c
          // +,-
          fprintf(out_bfmt7[thread_num],"%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t+-\n",
            b_array[index+foo].qid,
            //names[b_array[index+foo].qid-1].name,
            b_array[index+foo].qstt,
            b_array[index+foo].qend,
            b_array[index+foo].sid,
            //names[b_array[index+foo].sid-1].name,
            b_array[index+foo].sstt,
            b_array[index+foo].send,
            b_array[index+foo].aliscore,
            b_array[index+foo].evalue,
            b_array[index+foo].pident,
            b_array[index+foo].qseq,
            b_array[index+foo].sseq,
            b_array[index+foo].qlen,
            b_array[index+foo].slen
          );
        }
        else{
          // +,+
          fprintf(out_bfmt7[thread_num],"%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t++\n",
            b_array[index+foo].qid,
            //names[b_array[index+foo].qid-1].name,
            b_array[index+foo].qstt,
            b_array[index+foo].qend,
            b_array[index+foo].sid,
            //names[b_array[index+foo].sid-1].name,
            b_array[index+foo].sstt,
            b_array[index+foo].send,
            b_array[index+foo].aliscore,
            b_array[index+foo].evalue,
            b_array[index+foo].pident,
            b_array[index+foo].qseq,
            b_array[index+foo].sseq,
            b_array[index+foo].qlen,
            b_array[index+foo].slen
          );
        }
      }
      reversecomplementstr(&qseq[(LSEQ+1)*thread_num],qlen);
      fprintf(out_bfmt7[thread_num],"%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t--\n",
        i+2,1,(int)qlen,i+2,1,(int)qlen,(int)qlen,0.0,1.0,&qseq[(LSEQ+1)*thread_num],&qseq[(LSEQ+1)*thread_num],(int)qlen,(int)qlen);
      for(foo=0; foo<l_offset; ++foo){
        if((b_array[index+foo].sid-1) & 1){// sseq is rev.c
          // +,- -> -,+
          reversecomplementstr(b_array[index+foo].qseq,strlen(b_array[index+foo].qseq));
          reversecomplementstr(b_array[index+foo].sseq,strlen(b_array[index+foo].sseq));
          fprintf(out_bfmt7[thread_num],"%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t-+\n",
            b_array[index+foo].qid+1,
            //names[b_array[index+foo].qid-1].name,
            b_array[index+foo].qlen - b_array[index+foo].qend + 1,
            b_array[index+foo].qlen - b_array[index+foo].qstt + 1,
            b_array[index+foo].sid-1,
            //names[b_array[index+foo].sid-1].name,
            b_array[index+foo].slen - b_array[index+foo].send + 1,
            b_array[index+foo].slen - b_array[index+foo].sstt + 1,
            b_array[index+foo].aliscore,
            b_array[index+foo].evalue,
            b_array[index+foo].pident,
            b_array[index+foo].qseq,
            b_array[index+foo].sseq,
            b_array[index+foo].qlen,
            b_array[index+foo].slen
          );
        }
        else{
          // +,+ -> -,-
          reversecomplementstr(b_array[index+foo].qseq,strlen(b_array[index+foo].qseq));
          reversecomplementstr(b_array[index+foo].sseq,strlen(b_array[index+foo].sseq));
          fprintf(out_bfmt7[thread_num],"%llu\t%d\t%d\t%llu\t%d\t%d\t%d\t%.2lf\t%.2lf\t%s\t%s\t%d\t%d\t--\n",
            b_array[index+foo].qid+1,
            //names[b_array[index+foo].qid-1].name,
            b_array[index+foo].qlen - b_array[index+foo].qend + 1,
            b_array[index+foo].qlen - b_array[index+foo].qstt + 1,
            b_array[index+foo].sid+1,
            //names[b_array[index+foo].sid-1].name,
            b_array[index+foo].slen - b_array[index+foo].send + 1,
            b_array[index+foo].slen - b_array[index+foo].sstt + 1,
            b_array[index+foo].aliscore,
            b_array[index+foo].evalue,
            b_array[index+foo].pident,
            b_array[index+foo].qseq,
            b_array[index+foo].sseq,
            b_array[index+foo].qlen,
            b_array[index+foo].slen
          );
        }
      }
    }
  }

  }
  fprintf(stderr, "ol done\n");
//  fprintf(stderr,"%d\n",conf);


  fclose(fp);

  free(allreads);
  //free(rev);
  free(bounds);
  free(len_kmers);
  free(buf);
  free(kmers);
  free(sketches);
  free_bfmt7array(b_array,NUM_THREADS);
  free(qseq);
  free(sseq);
  free(tmp_qbuf);
  free(tmp_sbuf);
  free(tmp_abuf);
  free(res_qbuf);
  free(res_sbuf);
  free(res_abuf);
  free(sbj_buf);
  free(que_buf);
  free(ali_buf);
  {
    int i;
    for(i=0; i<NUM_THREADS; ++i){
      free_h64(&D0ss[i],&HPss[i],&VPss[i],&HNss[i],&VNss[i]);
    }
  }
  free(p_kmers);
  free(candidates);
  free(freq);
  free(names);
  free(Gx);
  if(opt_print_Ja){
    free(Jmat);
  }
  if(!opt_print_Jmatrix)
  {
    int i;
    char buffer[256];
    char prefix[256];
    strcpy(prefix, in_fa);
    for(i=strlen(prefix)-1; i>0; --i){
      if(prefix[i] == '.'){
        prefix[i] = '\0';
        break;
      }
    }
    for(i=0; i<NUM_THREADS; ++i){
      sprintf(buffer, "%s.%04d.bfmt7",prefix,i);
      fclose(out_bfmt7[i]);
    }
  }
  free(out_bfmt7);
  return 0;
}

void kmerize(char * reads, unsigned long long stt, unsigned long long end, record_t * kmers, int id, char strand){
  // kmer <- stt.. end-1 or [stt, end)
  unsigned long long kalph = ntuple_code(reads,stt,K);
  kmers[0].kmer = kalph;
  kmers[0].id = id;
  kmers[0].offset = 0;
  kmers[0].strand = strand;
  kmers[0].f_valid = 1;

  //fprintf(stderr, "strand: %c\n",strand);

  unsigned long long mask=1ull;
  if(K<32){
    mask<<=2*K;
    mask -= 1ull;
  }
  else{
    mask = ~0ull;
  }
  unsigned long long i;
  unsigned long long n_kmers = end-stt-K+1;
  if(n_kmers > end){
    fprintf(stderr,"uf stt:%llu, end:%llu\n",stt,end);
    exit(1);
  }
  for(i=1; i<n_kmers; ++i){
    kalph <<= 2;
    kalph &= mask;
    kalph |= to_code[(int)reads[stt+(unsigned long long)K-1ull+i]];
  //fprintf(stderr,"aiu\n");
//    if(*l_kmers+i >= s_kmers){
//      fprintf(stderr, "%llu %llu\n",*l_kmers+i,s_kmers);
//      exit(1);
//    }
    kmers[i].kmer = kalph;
    kmers[i].id = id;
    kmers[i].offset = i;
    kmers[i].strand = strand;
    kmers[i].f_valid = 1;
  //fprintf(stderr,"kkk\n");
  }
//  *l_kmers = *l_kmers+imax+1ull;
}

int compare_kmer_records(const void *a, const void *b){
    record_t * foo = *(record_t**)a;
    record_t * bar = *(record_t**)b;
//    printf("%llu\n",foo);
//    printf("%llu\n",bar);
    if(foo->kmer < bar->kmer){
      return -1;
    }
    else if(foo->kmer == bar->kmer){
      return 0;
    }
    else{
      return 1;
    }
}

void make_sketches(uint32_t * ske, unsigned long long lb){
}
int get_Jaccard_similarity_like(unsigned long long i, unsigned long long j, uint32_t * sketches){
  int n_matches = 0;
  int s;
  for(s=0; s<H_SIZE; ++s){
    if(sketches[H_SIZE*i+s] == sketches[H_SIZE*j+s]){
      ++n_matches;
    }
  }
  return n_matches;
}

int f_hit(unsigned long long target, unsigned long long * candidates, int len){
  if(len <= 0){
    return 0;// false
  }
  int index = len/2;
  if(candidates[index]>target){
    return f_hit(target,candidates,index);
  }
  else if(candidates[index]<target){
    return f_hit(target,candidates+index+1,len-index-1);
  }
  else if(candidates[index] == target){
    return 1;//true
  }
  else{
    fprintf(stderr, "buggy: f_hit\n");
    exit(1);
  }
}

void checker(char mark, char * orig, char * ali){
  int lo = strlen(orig);
  int la = strlen(ali);
  int i;
  int pos=0;
  fprintf(stderr,"\nali: %s\n",ali);
  for(i=0; i<la; ++i){
    if(ali[i] != '-'){
      ali[pos] = ali[i];
      ++pos;
    }
  }
  ali[pos] = '\0';
  if(pos > lo){
    pos = lo;
    ali[pos] = '\0';
    fprintf(stderr,"pos > lo: %d, %d\n",pos,lo);
  }
  for(i=0; i<pos; ++i){
    if(ali[i] != orig[i]){
      fprintf(stderr,"original\n%s\n%s\n\n",orig,ali);
      orig[i+1] = '\0';
      fprintf(stderr,"mark: %c\n",mark);
      fprintf(stderr,"org: %s\nali: %s\n",orig,ali);
      fprintf(stderr,"pos: %d\n",pos);
      fprintf(stderr,"i: %d\n",i);
      exit(1);
    }
  }
  fprintf(stderr,"ok: %c\n",mark);
}

int compare_kv(const void *a, const void *b){
  kv_t * a1 = (kv_t*)a;
  kv_t * b1 = (kv_t*)b;
  return (b1->value - a1->value);
}



