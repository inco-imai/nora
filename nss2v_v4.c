#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h> /* toupper */
#include <assert.h>

int LSEQ = 65536;
#define LAUX 32

/* ~hit num */
#define NROWS 512
#define LISTSIZE (1024)
#define FoS 16

void* 
malloc_or_die(size_t size, const char* message)
{
  void *pval;
  pval = malloc(size); 
  if(pval == NULL){ 
    fputs("Fatal error!\n", stderr);
    fputs("cannot allocate memory: ", stderr);
    fputs(message, stderr);
    fputs("\n", stderr);
    exit(1); 
  } 
  return pval;
}
void* 
realloc_or_die(void*buf, size_t size, const char* message)
{
  void *pval;
  pval = realloc(buf, size); 
  if(pval == NULL){ 
    fputs("Fatal error!\n", stderr);
    fputs("cannot reallocate memory: ", stderr);
    fputs(message, stderr);
    fputs("\n", stderr);
    exit(1); 
  } 
  return pval;
}


void* 
calloc_or_die(size_t nmemb, size_t size, const char* message)
{
  void *pval;
  pval = calloc(nmemb, size); 
  if(pval == NULL){ 
    fputs("Fatal error!\n", stderr);
    fputs("cannot allocate memory: ", stderr);
    fputs(message, stderr);
    fputs("\n", stderr);
    exit(1); 
  } 
  return pval;
}
void
reset_i_list(int *i_list, size_t nrows)
{
  int i;
  for(i = 0; i < nrows; ++i){
    i_list[i]=0;
  }
}

const unsigned INTERVAL=16;/* >0ull */
const int max_ins = 50;/* if a record has more than 31 insertions, this program will discard the record. */
const int max_depth = FoS*32;

typedef struct regions_t{
  unsigned long long stt;
  unsigned long long end;
}regions_t;

typedef struct cigar_t{
  int num;
  char sym;
}cigar_t;

typedef struct cmaux_t{
  char sym;
  int num;
}cmaux_t;

typedef struct sam_t{
  char * qname;
  int flag;
  char * rname;
  int pos;
  int mapq;
  char *cigar_string;
  struct cigar_t * cigar;
  int cigar_length;
  char * rnext;
  int pnext;
  int tlen;
  char * seq;
  char * qual;
  char as[LAUX];
  char ev[LAUX];
  int bit_score;
  double evalue;
  char * buf;
  size_t buffer_size;
}sam_t;

typedef struct base_t{
  char base;
  char qv;
}base_t;

typedef struct el_t{
  struct base_t el;
  struct el_t * down;
}el_t;
typedef struct el_pool_t{
  size_t index;
  size_t size; /* currently not exceeding 2^31 */
  el_t * free_el_list;
  struct el_pool_t *ext;
}el_pool_t;

el_pool_t*
init_el_pool(el_pool_t * pool_ptr, size_t size)
{
  pool_ptr->size = size;
  pool_ptr->index = 0;
  pool_ptr->ext = NULL;
  pool_ptr->free_el_list = 
   (el_t*)malloc_or_die(sizeof(el_t)*size, "free_el_list");
  return pool_ptr;
}

int free_el_pool(el_pool_t * pool_ptr)
{
  int size;
  if(pool_ptr == NULL)
    return 0;
  size = pool_ptr->size; 
  size += free_el_pool(pool_ptr->ext);
  pool_ptr->ext = NULL;
  free(pool_ptr->free_el_list);
  return size;
}

el_pool_t*
flush_el_pool(el_pool_t * pool_ptr)
{
  int size;
  size = free_el_pool(pool_ptr->ext);
  pool_ptr->ext = NULL;
  if(size){
    pool_ptr->size += size;
    pool_ptr->free_el_list = 
      realloc_or_die(pool_ptr->free_el_list, 
      sizeof(el_t)*(pool_ptr->size), "free_el_list");
  }
  pool_ptr->index = 0;
  return pool_ptr;
}

typedef struct column_t{
  int score;
  int s_depth;
  struct column_t * next;
  el_t * to_el;
}column_t;

typedef struct b_col_t{
  size_t size;
  int end;
  column_t * array;
}b_col_t;

b_col_t*
init_expand_b_col(b_col_t * basic_column, size_t end)
{
  int i;
  if(basic_column->size <= end){
    basic_column->array = realloc_or_die(basic_column->array, sizeof(column_t)*(end+1), "expand_b_col");
  }
  for(i=basic_column->end + 1; i <= end; ++i){
    basic_column->array[i].to_el = NULL;
    basic_column->array[i].next = NULL;
    basic_column->array[i].score = 0;
    basic_column->array[i].s_depth = 0;
  }
  basic_column->end = end;
  return basic_column;
}
typedef struct col_pool_t{
  unsigned long index;
  unsigned long size;/* = (LSEQ*(max_ins+1)); */
  column_t* free_col_list;
  struct col_pool_t * ext; 
}col_pool_t;

el_t*
init_el(el_t* el)
{
  el->el.base = ' ';
  el->el.qv = ' ';
  el->down = NULL;
  return el;
}

col_pool_t*
init_col_pool(col_pool_t* pool_ptr, size_t size)
{
  pool_ptr->size = size;
  pool_ptr->index = 0;
  pool_ptr->ext = NULL;
  pool_ptr->free_col_list = 
   (column_t*)malloc_or_die(sizeof(column_t)*size, "free_col_list");
  return pool_ptr;
}
int
free_col_pool(col_pool_t* pool_ptr)
{
  int size;
  if(pool_ptr == NULL)
    return 0;
  size = pool_ptr->size; 
  size += free_col_pool(pool_ptr->ext);
  pool_ptr->ext = NULL;
  free(pool_ptr->free_col_list);
  return size;
}
col_pool_t*
flush_col_pool(col_pool_t* pool_ptr)
{
  int size;
  size = free_col_pool(pool_ptr->ext);
  pool_ptr->ext = NULL;
  if(size){
    pool_ptr->size += size;
    pool_ptr->free_col_list =
      realloc_or_die(pool_ptr->free_col_list, 
        sizeof(column_t)*(pool_ptr->size), "free_col_list");
  }
  pool_ptr->index = 0;
  return pool_ptr;
}

el_t * 
get_el(el_pool_t *el_pool)
{
  if(el_pool->index >= el_pool->size){
    if(el_pool->ext == NULL){
      el_pool->ext = (el_pool_t*)malloc_or_die(sizeof(el_pool_t),"el_pool_expand");
      init_el_pool(el_pool->ext, el_pool->size * 1.4);
    }
    return get_el(el_pool->ext);
  }
  if(el_pool->free_el_list == NULL){
    fprintf(stderr, "free_el_list is NULL\n");
    exit(1);
  }
  return init_el(el_pool->free_el_list+el_pool->index++);
}

static inline int is_nt_char(char c)
{
  if(c == 'A' || c == 'C' || c == 'G' || c== 'T' || c == 'N')
    return 1;
  return 0;
}


column_t *
init_col(column_t *col, el_pool_t *el_pool)
{
  col->score = 0;
  col->s_depth = 0;
  col->next = NULL;
  col->to_el = get_el(el_pool);
  return col;
}
column_t * 
get_col(col_pool_t *col_pool, el_pool_t*el_pool){
  if(col_pool->index >= col_pool->size){
    if(col_pool->ext == NULL){
      col_pool->ext = (col_pool_t*)malloc_or_die(sizeof(col_pool_t),"col_pool_expand");
      init_col_pool(col_pool->ext, col_pool->size * 1.4);
    }
    return get_col(col_pool->ext, el_pool);
  }
  if(col_pool->free_col_list == NULL){
    fprintf(stderr, "free_col_list is NULL\n");
    abort();
  }
  return init_col(col_pool->free_col_list+col_pool->index++, el_pool);
}

regions_t ** alloc_lists_or_die(size_t nrows, size_t listsize)
{
  int i;
  char buf[32];
  regions_t ** lists;
  lists = (regions_t**)malloc_or_die(sizeof(regions_t*)*nrows, "lists");
  for(i=0;i<NROWS;++i){
    sprintf(buf, "lists[%i]", i);
    lists[i] = (regions_t*)malloc_or_die(sizeof(regions_t)*listsize, buf);
  }
  return lists;
}
void init_sam(sam_t * s, size_t LBUF){
  s->buffer_size = LBUF;
  s->qname = (char*)malloc_or_die(LBUF, "s->qname");
  s->rname = (char*)malloc_or_die(LBUF, "s->rname");
  s->cigar_string = (char*)malloc_or_die(LBUF, "s->rname");
  s->cigar = (cigar_t*)malloc_or_die(sizeof(cigar_t)*LSEQ, "s->cigar");
  s->rnext = (char*)malloc_or_die(LBUF, "s->rnext");
  s->seq = (char*)malloc_or_die(LBUF, "s->seq");
  s->qual = (char*)malloc_or_die(LBUF, "s->qual");
  s->buf = (char*)malloc_or_die(LBUF, "s->buf");
  return;
}

void realloc_sam(sam_t * s, size_t LBUF){
  s->buffer_size = LBUF;
  s->qname = (char*)realloc_or_die(s->qname, LBUF, "s->qname");
  s->rname = (char*)realloc_or_die(s->rname, LBUF, "s->rname");
  s->cigar_string = (char*)realloc_or_die(s->cigar_string, LBUF, "s->rname");
  s->cigar = (cigar_t*)realloc_or_die(s->cigar, sizeof(cigar_t)*LBUF/2, "s->cigar");
  s->rnext = (char*)realloc_or_die(s->rnext, LBUF, "s->rnext");
  s->seq = (char*)realloc_or_die(s->seq, LBUF, "s->seq");
  s->qual = (char*)realloc_or_die(s->qual, LBUF, "s->qual");
  s->buf = (char*)realloc_or_die(s->buf, LBUF, "s->buf");
  return;
}

void reset_sam(sam_t * s){
  s->qname[0] = '\0';
  s->rname[0] = '\0';
  s->cigar_length=0;
  s->rnext[0] = '\0';
  s->seq[0] = '\0';
  s->qual[0] = '\0';
  return;
}

void free_sam(sam_t * s){
  free(s->qname);
  free(s->rname);
  free(s->cigar_string);
  free(s->cigar);
  free(s->rnext);
  free(s->seq);
  free(s->qual);
  return;
}

int opt_separate=0;

void print_vertical(char *, b_col_t *, char* ,char *, el_pool_t*);
int ref_count = 0;

void
parse_cigar(sam_t*t_sam)
{
  int stt_cigar_index;
  int end_cigar_index;
  int offset=0;
  int j;
  int ret;
  int tmpoffset;
  int len = strlen(t_sam->cigar_string);
  for(j=0; offset < len && j<t_sam->buffer_size/2; ++j){
    ret = sscanf((t_sam->cigar_string+offset),"%d%c%n",&t_sam->cigar[t_sam->cigar_length].num,&t_sam->cigar[t_sam->cigar_length].sym,&tmpoffset);
    if(ret != 2){
      printf("%s\n",t_sam->cigar_string);
      fprintf(stderr, "souteigai 1 ret: %d %s\n",ret, t_sam->cigar_string+offset);
      fprintf(stderr, "offset: %d \n",offset);
      abort();
    }
    offset += tmpoffset;
    ++t_sam->cigar_length;
  }
  stt_cigar_index=0;
  end_cigar_index=t_sam->cigar_length-1;
  if(offset != len){
    fprintf(stderr, "souteigai 2 %d %d\n", offset, len);
    abort();
  }
  if(t_sam->cigar[stt_cigar_index].sym == 'H'){
    fprintf(stderr,"err: heading H\n");
    exit(1);
  }
  if(t_sam->cigar[end_cigar_index].sym == 'H'){
    fprintf(stderr,"err: tailing H\n");
    exit(1);
  }
  if(t_sam->cigar[stt_cigar_index].sym == 'I'){
    fprintf(stderr,"err: heading I\n");
    exit(1);
  }
  if(t_sam->cigar[end_cigar_index].sym == 'I'){
    fprintf(stderr,"err: tailing I\n");
    exit(1);
  }
}

void
normalize_naseq(char * seq)
{
  /* small -> capital */
  /* later, small letters will be ignored */
  int j;
  int loop = strlen(seq);
  for(j=0; j<loop; ++j){
    seq[j] = toupper(seq[j]);
    if(!is_nt_char(seq[j])){
      fputs("Fatal error!\n", stderr);
      fputs("unexpected character:", stderr);
      fprintf(stderr, "'%c' in\n%s\n",seq[j],seq);
      abort();
    }
  }
}

char*
gen_inserted_seq(char *sbuf, sam_t * t_sam)
{
  int cp=0;
  int si=0;
  int j, k;
  for(j=0; j<t_sam->cigar_length; ++j){
    char sym = t_sam->cigar[j].sym;
    if(sym == 'D'){
      for(k=0; k<t_sam->cigar[j].num; ++k){
        sbuf[si++] = '-';
      }
    }else if(sym == 'I'){
      for(k=0; k<t_sam->cigar[j].num; ++k){
        sbuf[si++] = t_sam->seq[cp++];
      }
    }else if(sym == 'M'){
      for(k=0; k<t_sam->cigar[j].num; ++k){
        sbuf[si++] = t_sam->seq[cp++];
      }
    }else{
      fprintf(stderr, "strange cigar %c\n",sym);
      exit(1);
    }
  }
  sbuf[si] = '\0';
  return sbuf;
}

int
get_sam_record(sam_t*t_sam, FILE* in_sam)
{
  char * retp;
  int sretval;
  int line_length;
  char * buf = t_sam->buf;
  size_t buffer_size = t_sam->buffer_size;
  retp = fgets(buf, buffer_size, in_sam);
  if(!retp) return 0;
  line_length = strlen(buf);
  if(line_length >= t_sam->buffer_size - 1){
    /* The line is incompletely read if the buffer 
     * is fully occuppied and not ending with newline */
    while(line_length >= buffer_size - 1 && buf[line_length - 1] != '\n'){
      buf = (char*)realloc_or_die(buf, buffer_size*2, "buf");
      fgets(buf+line_length, buffer_size + 1, in_sam);
      line_length = strlen(buf);
      buffer_size *= 2;
    }
    t_sam->buf=buf;
/* now that the length of the buffer was measured, 
 * the line should contain qseq and sseq of the same length.
 * Thus the sequence should be less than half of the line length */ 
    realloc_sam(t_sam, buffer_size);
  }
  reset_sam(t_sam);
  t_sam->cigar_length=0;
  sretval = sscanf(t_sam->buf,"%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
        t_sam->qname,
        &t_sam->flag,
        t_sam->rname,
        &t_sam->pos,  /* 1-origin */
        &t_sam->mapq,
        t_sam->cigar_string,
        t_sam->rnext,
        &t_sam->pnext,
        &t_sam->tlen,
        t_sam->seq,
        t_sam->qual,
        t_sam->as,
        t_sam->ev);
        t_sam->tlen = strlen(t_sam->seq);
  if(t_sam->pos > 0){
    t_sam->pos--;/* 0-origin */
  }else{
    fprintf(stderr,"strange pos: %d . this must be > 0\n",t_sam->pos);
    exit(1);
  }
  if(sretval != 13){/* XXX 13 may be changed */
    fprintf(stderr, "souteigai sam: retval %d\n",sretval);
    fprintf(stderr, "%s\n", t_sam->buf);
    abort();
  }
  if(strcmp("*",t_sam->qual) == 0){
    int len = strlen(t_sam->seq);
    int j;
    for(j=0; j<len; ++j){
      t_sam->qual[j] = '!';
    }
  }
  sretval = sscanf(t_sam->as, "AS:i:%d", &(t_sam->bit_score));
  if(sretval != 1){
    fprintf(stderr, "no alignment score ? : retval %d\n",sretval);
    fprintf(stderr, "%s\n",t_sam->as);
    fprintf(stderr, "%s\n",t_sam->qname);
    fprintf(stderr, "%s\n",t_sam->rname);
    abort();
  }
  sretval = sscanf(t_sam->ev, "EV:Z:%lf", &(t_sam->evalue));
  if(sretval != 1){
    fprintf(stderr, "no evalue ? : retval %d\n",sretval);
    fprintf(stderr, "%s\n",t_sam->ev);
    abort();
  }
  parse_cigar(t_sam);
  normalize_naseq(t_sam->seq);
  return 1;
}
int
calculate_region_length(sam_t* t_sam)
{
  int j;
  int length_of_region = 0;
  for(j = 0; j < t_sam->cigar_length; j++){
    char sym = t_sam->cigar[j].sym;
    if(sym == 'M' || sym == 'D'){
      length_of_region += t_sam->cigar[j].num;
    }else if(sym == 'I' && t_sam->cigar[j].num > max_ins){
      /* max_ins is a global variable */
      fprintf(stderr, "WARNING: too large ins: %d . skipped\n",t_sam->cigar[j].num);
      return 0;
    }
  }
  return length_of_region;
}
int
find_row_with_space(regions_t**lists, int* i_list, int listsize, int nrows, int opt_separate, int stt, int end, int interval)
{
  int row;
  for(row=0; row < nrows; ++row){
    int f_free = 0;
    /* find a free region */
    if(!opt_separate){
      int l;
      int not_collide = 1;
      for(l=0; l<i_list[row] && l<listsize; ++l){/* check the rgn does not overlap any reads */
        if(end+interval < lists[row][l].stt || stt > lists[row][l].end+interval){
          continue;
        }else{
          not_collide=0;
        }
      }
      if(not_collide==1 && i_list[row] <= listsize){
        f_free = 1;
      }
    }
    if(opt_separate){
      if(i_list[row] == 0){
        f_free = 1;
      }
    }
    if(f_free){
      return row;
    }
  }
  return row;
}
int
paste_match(b_col_t * basic_column, el_pool_t * el_pool, int *pos, int *si, const char *sbuf, int num, int row, int stt)
{
  /* basic_column is the major operand that will be changed */
  /* el_pool supplies element cells */

  /* num is referred only once */
  /* pos is incremented for num times*/
  /* k is referred once */
  /* row is not changed */
  /* si is incremented */
  /* stt is not changed */
  /* sbuf is read only */

  int t = *pos+num;
  int l;
  int t_si = *si;
  for(l=*pos; l<t; ++l,++*pos){
    int score = 0;
    el_t * tmp;
    int m;
    if(basic_column->array[l].to_el == NULL){
      /* initialize */
      basic_column->array[l].to_el = get_el(el_pool);
      basic_column->array[l].next = NULL;
      basic_column->array[l].score = 0;
    }
    tmp = basic_column->array[l].to_el;
    for(m=0; m<row; ++m){
      if(tmp->down == NULL){
        tmp->down = get_el(el_pool);
      }
      if(tmp->el.base != ' '){
        ++score;
      }
      tmp = tmp->down;
    }
    if(sbuf[t_si] != ' '){
      assert(tmp->el.base == ' ');
      tmp->el.base = sbuf[t_si];
      ++score;
    }else{
      if(tmp->el.base != ' '){
        ++score;
      }
    }
    ++t_si;
    while(tmp->down){
      tmp = tmp->down;
      if(tmp->el.base != ' '){
        ++score;
      }
    }
    basic_column->array[l].score = score;

    if(*pos != stt && basic_column->array[l].to_el->el.base != ' ' && basic_column->array[l-1].to_el->el.base != ' ')
    {
      column_t * cur_col;
      column_t * l_col;
      column_t * r_col;
      assert(l >= 1);
      cur_col = basic_column->array[l-1].next;/* ins column */
      l_col = &(basic_column->array[l-1]);
      r_col = &(basic_column->array[l]);
      for(; cur_col; cur_col=cur_col->next){
        {
          el_t * l_el = l_col->to_el;
          el_t * r_el = r_col->to_el;
          el_t * c_el = cur_col->to_el;
          int depth = 0;
          int score = 0;
          while(l_el != NULL && r_el != NULL)
          {
            if(c_el->el.base == ' '){
              if(l_el->el.base != ' ' && r_el->el.base != ' '){
                c_el->el.base = '-';
                ++score;
              }
            }else{
              ++score;
            }
            ++depth;

            l_el = l_el->down;
            r_el = r_el->down;
            if(depth > row && (l_el == NULL || r_el == NULL)){
              break;
            }
            if(c_el->down == NULL){
              c_el->down = get_el(el_pool);
            }
            c_el = c_el->down;
          }
          while(depth <= row){
            if(depth == row){
              c_el->el.base = '-';
              ++score;
              break;
            }else if(c_el->el.base == ' '){
            }else{
              ++score;
            }
            if(c_el->down == NULL){
              c_el->down = get_el(el_pool);
            }
            c_el = c_el->down;
            ++depth;
          }
          cur_col->score = score;
        }
      }
    }
  }
  *si = t_si;
  return *pos;
}
int
paste_insert(b_col_t * basic_column, el_pool_t * el_pool, col_pool_t *col_pool,  int pos, int *si, const char *sbuf, int num, int row, int stt)
{
  int p;
  column_t *l_col, *r_col, *cur_col;
  assert(pos >=1);
  l_col = &(basic_column->array[pos-1]);
  r_col = &(basic_column->array[pos]);
  cur_col = l_col;

  for(p=0; p< num; ++p){
    if(cur_col->next == NULL){
      cur_col->next = get_col(col_pool, el_pool);
    }
    cur_col = cur_col->next;
    {
      el_t * l_el = l_col->to_el;
      el_t * r_el = r_col->to_el;
      el_t * c_el = cur_col->to_el;
      int depth = 0;
      int score = 0;
      while(l_el != NULL && r_el != NULL)
      {
        if(depth == row){
          if(sbuf[*si] != ' '){
            assert(c_el->el.base == ' ');
            c_el->el.base = sbuf[*si];
            ++score;
          }else{
            if(c_el->el.base != ' '){
              ++score;
            }
          }
          ++*si;
        }else if(c_el->el.base == ' '){
          if(l_el->el.base != ' ' && r_el->el.base != ' '){
            c_el->el.base = '-';
            ++score;
          }
        }else{
          /* keep */
          ++score;
        }
        ++depth;

        l_el = l_el->down;
        r_el = r_el->down;
        if(depth > row && (l_el == NULL || r_el == NULL)){
          break;
        }
        if(c_el->down == NULL){
          c_el->down = get_el(el_pool);
        }
        c_el = c_el->down;
      }
      while(depth <= row){
        if(depth == row){
          if(sbuf[*si] != ' '){
            c_el->el.base = sbuf[*si];
            ++score;
          }
          ++*si;
          break;
        }else if(c_el->el.base == ' '){
        }else{
          ++score;
        }
        if(c_el->down == NULL){
          c_el->down = get_el(el_pool);
        }
        c_el = c_el->down;
        ++depth;
      }
      cur_col->score = score;
    }
  }
  return *si;
}
int main(int argc, char ** argv)
{
  int hitnum=0;
  int LBUF;
  FILE * IN_SAM = NULL;
  int valid_voters = (NROWS-1);
  char *prevchr;
  b_col_t basic_column;
  int * i_list;
  char *pbuf;
  char *sbuf;
  int max_end=0; /* only int variable end2 is assigned */
  el_pool_t el_pool;
  col_pool_t col_pool;
  regions_t ** lists;
  sam_t * t_sam;
  {
    int result;
    int tmp;
    while((result=getopt(argc,argv,"r:s")) != -1){
      switch(result){
        case 's':
          opt_separate=1;
          ++hitnum;
          break;
        case 'r':
          tmp = atoi(optarg);
          if(tmp <= 0){
            fprintf(stderr, "r(max_read_length) must be > 0: %s\n", optarg);
            exit(1);
          }
          LSEQ = tmp;
          if(tmp>=65536){
            ++LSEQ;
          }
          hitnum += 2;
          break;
        case '?':
          printf("humei\n");
          exit(1);
          break;
        default:
          break;
      }
    }
  }

  LBUF = 2*LSEQ;
    if(argc != 2+hitnum)
    {
      fprintf(stderr, "USAGE: <this> <in.nss>\n");
      fprintf(stderr, "\tinput MUST be sorted by 1st chromosome name (must), 2nd evalue and 3rd bitScore\n");
      fprintf(stderr, "\t-r <integer>: max read length. default %d\n",LSEQ);
      return 1;
    }
  if(LSEQ <= 65536){
    el_pool.size = LSEQ*(1+max_ins);
    el_pool.size *= max_depth/FoS; /* 32 */
  }else{
    el_pool.size = LSEQ*100;
  }
  col_pool.size = LSEQ*(1+max_ins);

  if(strcmp("-",argv[1+hitnum])==0){
    IN_SAM = stdin;
  }else{
    IN_SAM = fopen(argv[1+hitnum],"r");
    if(IN_SAM == NULL){
      fprintf(stderr, "cannot open the file %s\n", argv[1+hitnum]);
      abort();
    }
  }
  prevchr = (char*)malloc_or_die(LBUF, "prevchr");
  basic_column.array = (column_t*)malloc_or_die(sizeof(column_t)*LSEQ, "basic_column");
  basic_column.size = LSEQ;
  basic_column.end=-1;
  init_el_pool(&el_pool, el_pool.size);
  init_col_pool(&col_pool, col_pool.size);
  lists = alloc_lists_or_die(NROWS, LISTSIZE);
  i_list = (int*)calloc_or_die(NROWS, sizeof(int), "i_list");
  pbuf = (char*)malloc_or_die(LSEQ, "pbuf");
  sbuf = (char*)malloc_or_die(LSEQ*2, "sbuf");
  t_sam = (sam_t*)malloc_or_die(sizeof(sam_t), "t_sam");
  init_sam(t_sam, LBUF);

  max_end=0;
  prevchr[0] = '\0';
  while(get_sam_record(t_sam, IN_SAM)){
    if(LSEQ < t_sam->buffer_size/2){
      LSEQ=(t_sam->buffer_size+1)/2;
      pbuf = (char*)realloc_or_die(pbuf, LSEQ, "pbuf");
      sbuf = (char*)realloc_or_die(sbuf, LSEQ*2, "sbuf");
    }
    if(prevchr[0] == '\0'){
      strcpy(prevchr,t_sam->rname);
    }else if(strcmp(prevchr,t_sam->rname) != 0){
      if(ref_count>0){
        print_vertical(prevchr, &basic_column, pbuf,sbuf, &el_pool);
      }
      ref_count = 0;

      strcpy(prevchr, t_sam->rname);
      reset_i_list(i_list, NROWS);
      basic_column.end=-1;
      max_end = 0;
      flush_el_pool(&el_pool);
      flush_col_pool(&col_pool);
    }

    /* look for a free region and set */
    {
      int stt;
      int end;
      int stt2;
      int end2;
      int f_discard;
      int length_of_region;
      length_of_region = calculate_region_length(t_sam);
      if(!length_of_region){ continue; }
      stt = t_sam->pos;/* 0-originized */
      end = stt + length_of_region - 1;
      stt2 = stt;
      end2 = end;

      if(strcmp(t_sam->rname,t_sam->qname) == 0){/* this is a ref-ref hit */
        int length = strlen(t_sam->seq);
        if(length != t_sam->cigar[0].num){
          /*discard this record */
          continue;
        }
        if(basic_column.end < end){
          /* initialize basic_column[] */
          init_expand_b_col(&basic_column, end);
        }
        ++ref_count;
      }else{
        if(ref_count <= 0){
          continue;/* against chunks without ref-ref hits. */
        }
      }

      f_discard=1;
      {
        int s;
        for(s=stt; s<=end; ++s){
          basic_column.array[s].s_depth++;
          if(basic_column.array[s].s_depth <= valid_voters){
            f_discard=0;
          }
        }
      }
      if(f_discard){
        //continue;// commented out. // XXX
      }

      {
        int row = find_row_with_space(lists, i_list, LISTSIZE, NROWS, opt_separate, stt2, end2, INTERVAL);
        if(row >= NROWS) continue;
        /* column[row] is free */
        /* set */
        {
          int k;
          int pos=stt;
          int si = 0;
          gen_inserted_seq(sbuf, t_sam);
          for(k=0; k<t_sam->cigar_length; ++k){
            char sym = t_sam->cigar[k].sym;
            if(sym == 'M' || sym == 'D'){
              paste_match(&basic_column, &el_pool, &pos, &si, sbuf, t_sam->cigar[k].num, row, stt);
            }else if(sym == 'I'){
              paste_insert(&basic_column, &el_pool, &col_pool, pos, &si, sbuf, t_sam->cigar[k].num, row, stt);
            }else{
              fprintf(stderr,"strange sym %c\n",sym);
              exit(1);
            }
          }
        }
        lists[row][i_list[row]].stt = stt2; 
        lists[row][i_list[row]].end = end2; 
        ++i_list[row];
        if(max_end < end2){
          max_end = end2;
        }
      }
    }
  }

  if(ref_count>0){
    print_vertical(prevchr, &basic_column, pbuf,sbuf,&el_pool);
  }
  ref_count = 0;

  free_sam(t_sam);
  free(t_sam);
  if(IN_SAM)  fclose(IN_SAM);
  free(prevchr);
  free(basic_column.array);
  free_col_pool(&col_pool);
  free_el_pool(&el_pool);
  {
    int i;
    for(i=0;i<NROWS;++i){
      free(lists[i]);
    }
    free(lists);
  }
  free(i_list);
  free(pbuf);
  free(sbuf);

  return EXIT_SUCCESS;
}

int print_check(el_t * tmp, char * chr, int stt, int i_basic_column){
  if(tmp == NULL){
    fprintf(stderr, "NULL column 1 %s %d %d\n",chr,stt,i_basic_column);
    exit(1);
  }
  while(tmp && (tmp->el.base == ' ' || tmp->el.base == '-')){
    tmp = tmp->down;
  }
  if(tmp){
    return 1;/* will be printed */
  }else{
    return 0;
  }
}

void print_column(el_t * tmp, int d, char * pbuf){
  int index=0;
  if(tmp == NULL){
    fprintf(stderr, "NULL column 2\n");
    exit(1);
  }
  pbuf[index++]=tmp->el.base;
  while(tmp->down){
    tmp = tmp->down;
    pbuf[index++]=tmp->el.base;
  }
  pbuf[index] = '\0';
  if(d>=0){
    if(d==0){
      printf("%s %c %3d\n",pbuf,(char)(0+33),d);
    }else{
      if(d>=920){
        d = 920;
      }
      printf("%s %c %3d\n",pbuf,(char)(d/10+1+33),d);
    }
  }else{
    printf("%s\n",pbuf);
  }
}

void print_vertical(char * chr, b_col_t *basic_column,char * pbuf,char * sbuf, el_pool_t*el_pool){
  /* done valid voters */
  /* done terminal '-' */
  /* done chop short frags */

  int stt=basic_column->end+1;
  int end=0;

  int printed=0;
  int i;
  int si=0;
  for(i=0; i<=basic_column->end; ++i){
    el_t * tmp = basic_column->array[i].to_el;
    if(tmp){
      stt = i;
      break;
    }
  }
  for(i=basic_column->end; i>=stt; --i){
    el_t * tmp = basic_column->array[i].to_el;
    if(tmp){
      end = i;
      break;
    }
  }

  for(i=stt; i<=end; ++i){
    int d;
    el_t * tmp = basic_column->array[i].to_el;
    if(tmp == NULL){
      basic_column->array[i].to_el = get_el(el_pool);
      basic_column->array[i].to_el->el.base = 'N';
      basic_column->array[i].next = NULL;
      basic_column->array[i].score = 1;
      basic_column->array[i].s_depth = 1;
    }
    d = basic_column->array[i].s_depth;
    if(d==0){
      sbuf[si] = (char)(0+33);
    }else{
      if(d>=920){
        d = 920;
      }
      sbuf[si] = (char)(d/10+1+33);
    }
    ++si;
  }
  sbuf[si] = '\0';

  for(i=stt; i<=end; ++i){
    el_t* tmp = basic_column->array[i].to_el;
    column_t* tc;
    if(print_check(tmp,chr,stt,basic_column->end)){
      if(!printed){
        printf("%%%s\n",chr);
        //printf("#%s\n",sbuf); //TODO
      }
      print_column(tmp,-1,pbuf);
      printed=1;
    }
    tc = basic_column->array[i].next;
    while(tc){
      el_t * tmp = tc->to_el;
      if(print_check(tmp,chr,stt,basic_column->end)){
        if(!printed){
          printf("%%%s\n",chr);
        }
        print_column(tmp,-1,pbuf);
        printed=1;
      }
      tc = tc->next;
    }
  }
  if(printed)
    printf("\n");
}

