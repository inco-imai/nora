#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <limits.h>
#include <omp.h>

#define LNAME 256
int N_PLOIDY = 1;

int LBUF = 65536;

// ~depth
#define NROWS 256

#define BUFSIZE (4*NROWS)

#define LISTSIZE (1024)

#define N_NODES (10*1024*1024)
#define N_EDGES_PER_NODE (64)
unsigned long long N_EDGES;

unsigned long long max_read_length;
int NODEBUFMAX = 600000;


int term_delta = 42; // FIXME

int MIN_ALI_SCORE = 300;
int D_LARGE = 10000;

int MAX_IO_DEGREE = 100; //XXX

int debug=0;
int debug2=0;

typedef struct bl_t{
//  int Qid;
//  int Sid;
  char str_Qid[LNAME];
  char str_Sid[LNAME];
  int Qstt;
  int Qend;
  int Sstt;
  int Send;
  double aliscore;
  double evalue;
  double pident;
  int Qlen;
  int Slen;
}bl_t;

typedef struct edge_t{
  int valid;// 1: valid, 0: invalid
  int Qid;
  int Sid;
//  char str_Qid[LNAME];
//  char str_Sid[LNAME];
  char strand;// +:'+', -:'-'
  char r_or_l_or_a; // r:'r', l:'l', ambiguous: 'a'
  //     r           l      
  // q ----          ---- q 
  // s   ----      ----   s 
  int Qstt;
  int Qend;
  int Sstt;
  int Send;
  int u_sol;
  int ol_score;
  int will_be_deleted;
}edge_t;

typedef struct node_t{
  int valid;// 1: valid, 0: invalid, 2:term_flag (valid && term)
  int len;// read length
  unsigned long long si_edges[3];// the starting index of edges. i_edges[0] is for 'r', i_edges[1] is for 'l', i_edges[2] is for 'a'
  unsigned long long ei_edges[3];// the ending index of edges.
  int f_term;
}node_t;

typedef struct name_t{
  char name[256];
}name_t;

unsigned long long el_index=0ull;

unsigned long long col_index=0ull;

int weighted_mode=0;
int vertical_mode=1;
int depth_comment_mode=0;
int depth_filter_mode=0;
int evalue_filter_mode=0;
int bit_score_filter_mode=0;
int gz_mode=0;
int maxdepth = 50;
int valid_read_length_is_given=0;
int voter_limit_is_given=0;
int opt_chimeric_or_concatenated_filter=0;
int opt_debug=0;

void sort_edges(unsigned long long csi, char* pbuf,char * sbuf,edge_t * b_edges_1, edge_t * b_edges_2, edge_t * b_edges_3);
int ref_count;
void invalidate_reverse_edge(edge_t * e);
void node_checker();

node_t * nodes;
int l_nodes;
edge_t * edges;
unsigned long long l_edges;
unsigned long long csi_edges;
edge_t ** s_edges;
unsigned long long l_s_edges;

node_t * nodes_bak;
int l_nodes_bak;
edge_t * edges_bak;
unsigned long long l_edges_bak;

int unique_base_threshold=1;// XXX
int f_greedy = 0;
int f_pgs=0;

int * paths;

void f_forward(int qidx, unsigned long long currD, unsigned long long currE, int * paths, unsigned long long Emax2);
void f_backward(int qidx, unsigned long long currE, int * paths, unsigned long long Emax1);

void f_term_or_not(unsigned long long j, int * f_term);
void f_delete_node(unsigned long long j);

void print_graph_simple();
void print_graph_simple2();
int opt_bubblepopping = 0;
unsigned long long Dmax;

void check_edges(edge_t * e);
void invalidate_ambiguous_edges(unsigned long long j);

void clean_edges();
void swap_edges(unsigned long long a, unsigned long long b);

int get_max_depth(unsigned long long idx, int s, int limit);
int f_merge(unsigned long long idx, int s);
void will_invalidate_this_branch(unsigned long long idx, int s, int limit);
void delete_edges();

void copy_nodes(node_t * to, node_t * from, int len);
void copy_edges(edge_t * to, edge_t * from, unsigned long long len);

char opt_gfa_file_prefix[256];
FILE * OUT_GFA_ORIG;
FILE * OUT_GFA_INIT;
FILE * OUT_GFA_BP;
FILE * OUT_GFA_PRUNED;
char tmpbuf[500];

int * blacklist;

int get_degree(int idx);
int get_not_small_io_degree(int idx);

void invalidate_node(int idx);

char opt_blist[256];
FILE * IN_BLIST;
int f_blist=0;

void clean_graph();

void print_edge_population();
int opt_freq=0;

int get_number_of_io_edges(int idx, int direction);
void reduce_repeat(int idx, int direction);

int repeat_depth = 100;
int min_ali_score=0;

void reduce_edges(unsigned long long idx, int direction, int topx);

void loop_remover(unsigned long long idx, int Depth, unsigned long long root);

int main(int argc, char ** argv)
{
  int hitnum=0;
  N_EDGES = (N_EDGES_PER_NODE*(unsigned long long)N_NODES);
  opt_gfa_file_prefix[0] = '\0';
  opt_blist[0] = '\0';
  {
    int result;
    int tmp;
    while((result=getopt(argc,argv,"u:pd:gs:t:P:B:fr:a:")) != -1){
      switch(result){
        case 'd':
          tmp = atoi(optarg);
          if(tmp < 0){
            fprintf(stderr, "d(delta) must be >= 0: %s\n", optarg);
            exit(1);
          }
          term_delta = tmp;
          hitnum += 2;
          break;
        case 'g':
          f_greedy=1;
          ++hitnum;
          break;
        case 's':
          f_pgs=1;
          strcpy(opt_gfa_file_prefix,optarg);
          hitnum += 2;
          break;
        case 't':
          tmp = atoi(optarg);
          if(tmp < 0){
            fprintf(stderr, "t(transitive reduction depth) must be >= 0: %s\n", optarg);
            exit(1);
          }
          opt_bubblepopping = tmp;
          hitnum += 2;
          break;
        case 'P':
          N_PLOIDY = atoi(optarg);
          if(N_PLOIDY < 1){
            fprintf(stderr, "P (N_PLOIDY) must be >= 1: %s\n", optarg);
            exit(1);
          }
          hitnum += 2;
          break;
        case 'B':
          f_blist=1;
          strcpy(opt_blist,optarg);
          hitnum += 2;
          break;
        case 'f':
          opt_freq=1;
          ++hitnum;
          break;
        case 'r':
          tmp = atoi(optarg);
          if(tmp < 0){
            fprintf(stderr, "repeat_depth must be >= 0: %s\n", optarg);
            exit(1);
          }
          repeat_depth = tmp;
          hitnum += 2;
          break;
        case 'a':
          tmp = atoi(optarg);
          if(tmp < 0){
            fprintf(stderr, "a(min_ali_score) must be >= 0: %s\n", optarg);
            exit(1);
          }
          min_ali_score = tmp;
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
    if(evalue_filter_mode+depth_filter_mode+bit_score_filter_mode > 1){
      fprintf(stderr, "cannot set -e, -d and -b simultaneously\n");
      exit(1);
    }
  }

  if(opt_gfa_file_prefix[0] != '\0'){
    sprintf(tmpbuf,"%s.init.gfal",opt_gfa_file_prefix);
    OUT_GFA_INIT = fopen(tmpbuf,"w");
    if(OUT_GFA_INIT == NULL){
      fprintf(stderr, "cannot open the file %s\n", tmpbuf);
      exit(1);
    }
    sprintf(tmpbuf,"%s.bp.gfal",opt_gfa_file_prefix);
    OUT_GFA_BP = fopen(tmpbuf,"w");
    if(OUT_GFA_BP == NULL){
      fprintf(stderr, "cannot open the file %s\n", tmpbuf);
      exit(1);
    }
    sprintf(tmpbuf,"%s.pruned.gfal",opt_gfa_file_prefix);
    OUT_GFA_PRUNED = fopen(tmpbuf,"w");
    if(OUT_GFA_PRUNED == NULL){
      fprintf(stderr, "cannot open the file %s\n", tmpbuf);
      exit(1);
    }
    sprintf(tmpbuf,"%s.orig.gfal",opt_gfa_file_prefix);
    OUT_GFA_ORIG = fopen(tmpbuf,"w");
    if(OUT_GFA_ORIG == NULL){
      fprintf(stderr, "cannot open the file %s\n", tmpbuf);
      exit(1);
    }
  }

  if(opt_blist[0] != '\0'){
    IN_BLIST = fopen(opt_blist,"r");
    if(IN_BLIST == NULL){
      fprintf(stderr, "cannot open the file %s\n",opt_blist);
      exit(1);
    }
  }


    if(argc != 2+hitnum)
    {
      //fprintf(stderr, "USAGE: <this> <in.bl> <all.fasta> [-c olcutoff]\n");
      fprintf(stderr, "USAGE: <this> <in.bl>\n");
      fprintf(stderr, "\tin.bl MUST be sorted by 1st chromosome name (must), 2nd evalue and 3rd bitScore\n");
      //fprintf(stderr, "\t-g: use each left and right highest score only.\n");
      //fprintf(stderr, "\t-s: print graph\n");
      return 1;
    }

  char in_bfmy7[LNAME];
  strcpy(in_bfmy7, argv[1+hitnum]);

  FILE * IN_BFMT7;
  if(strcmp("-",in_bfmy7)==0){
    IN_BFMT7 = stdin;
  }
  else{
    if(gz_mode!=1){
      IN_BFMT7 = fopen(in_bfmy7,"r");
      if(IN_BFMT7 == NULL){
        fprintf(stderr, "cannot open the file %s\n", in_bfmy7);
        exit(1);
      }
    }
    else{
      fprintf(stderr, "cannot open the file %s\n", in_bfmy7);
      exit(1);
    }
  }

  char * buf = (char*)malloc(sizeof(char)*LBUF);
  if(buf == NULL){
    fprintf(stderr, "cannot allocate memory: buf\n");
    exit(1);
  }

  edges = (edge_t*)malloc(sizeof(edge_t)*N_EDGES);
  if(edges == NULL){
    fprintf(stderr, "cannot allocate memory: edges\n");
    exit(1);
  }
  l_edges = 0;
  csi_edges = 0;

  s_edges = (edge_t**)malloc(sizeof(edge_t*)*N_EDGES);
  if(s_edges == NULL){
    fprintf(stderr, "cannot allocate memory: s_edges\n");
    exit(1);
  }
  l_s_edges = 0;

  edge_t * b_edges_1 = (edge_t*)malloc(sizeof(edge_t)*LISTSIZE);
  if(b_edges_1 == NULL){
    fprintf(stderr, "cannot allocate memory: b_edges_1\n");
    exit(1);
  }
  edge_t * b_edges_2 = (edge_t*)malloc(sizeof(edge_t)*LISTSIZE);
  if(b_edges_2 == NULL){
    fprintf(stderr, "cannot allocate memory: b_edges_2\n");
    exit(1);
  }
  edge_t * b_edges_3 = (edge_t*)malloc(sizeof(edge_t)*LISTSIZE);
  if(b_edges_3 == NULL){
    fprintf(stderr, "cannot allocate memory: b_edges_3\n");
    exit(1);
  }
  edges_bak = (edge_t*)malloc(sizeof(edge_t)*N_EDGES);
  if(edges_bak == NULL){
    fprintf(stderr, "cannot allocate memory: edges_bak\n");
    exit(1);
  }
  l_edges_bak = 0;

  nodes = (node_t*)malloc(sizeof(node_t)*(unsigned long long)N_NODES);
  if(nodes == NULL){
    fprintf(stderr, "cannot allocate memory: nodes\n");
    exit(1);
  }
  l_nodes = 0;
  nodes_bak = (node_t*)malloc(sizeof(node_t)*(unsigned long long)N_NODES);
  if(nodes_bak == NULL){
    fprintf(stderr, "cannot allocate memory: nodes_bak\n");
    exit(1);
  }
  l_nodes_bak = 0;
  blacklist = (int*)malloc(sizeof(int)*(unsigned long long)N_NODES);
  if(blacklist == NULL){
    fprintf(stderr, "cannot allocate memory: blacklist\n");
    exit(1);
  }

  char * pbuf = (char*)malloc(sizeof(char)*LBUF);
  if(pbuf == NULL){
    fprintf(stderr, "cannot allocate memory: pbuf\n");
    exit(1);
  }

  char * sbuf = (char*)malloc(sizeof(char)*LBUF);
  if(sbuf == NULL){
    fprintf(stderr, "cannot allocate memory: sbuf\n");
    exit(1);
  }

  bl_t * t_bl;
  t_bl = (bl_t*)malloc(sizeof(bl_t));
  if(t_bl == NULL){
    fprintf(stderr, "cannot allocate memory: t_bl\n");
    exit(1);
  }

  int NUM_THREADS = omp_get_max_threads();
  paths = (int*)malloc(sizeof(int)*NODEBUFMAX*NUM_THREADS);
  if(paths == NULL){
    fprintf(stderr, "cannot allocate memory: paths\n");
    exit(1);
  }

  //name_t * names = (name_t*)malloc(sizeof(name_t)*N_NODES);
  //if(names == NULL){
   // fprintf(stderr, "cannot allocate memory: names\n");
    //exit(1);
  //}

  // read fasta db
  //char in_fa[LNAME];
  //strcpy(in_fa, argv[2+hitnum]);
  /*
  FILE * IN_FA;
  IN_FA = fopen(in_fa,"r");
  if(IN_FA == NULL){
    fprintf(stderr,"cannot open the file %s\n",in_fa);
    exit(1);
  }
  while(fgets(buf,LBUF,IN_FA) != NULL){
    if(buf[0] != '>'){
      fprintf(stderr,"strange fasta format (1)\n");
      exit(1);
    }
    int id;
    int pos;// TODO
    sscanf(&buf[1],"%d%n",&id,&pos);
    if(pos != strlen(buf)-2){
      fprintf(stderr,"strange fasta format\nname must be '>num\\n'\n");
      fprintf(stderr,"%s\n",buf);
      fprintf(stderr,"pos: %d len: %d\n",pos,(int)strlen(buf));
      exit(1);
    }
    if(id != l_nodes+1){
      fprintf(stderr,"input name must be >1, >2, >3, ...\n");
      exit(1);
    }
    if(fgets(buf,LBUF,IN_FA) == NULL){
      fprintf(stderr,"strange fasta format (2)\n");
      exit(1);
    }
    int len = strlen(buf);
    // chomp
    if(buf[len-1] == '\n'){
      buf[len-1] = '\0';
      --len;
    }

    // nodes[id-1] is for the data of the id
    nodes[l_nodes].valid = 0;
    nodes[l_nodes].len = len;
    nodes[l_nodes].si_edges[0] = ULLONG_MAX;
    nodes[l_nodes].si_edges[1] = ULLONG_MAX;
    nodes[l_nodes].ei_edges[0] = ULLONG_MAX;
    nodes[l_nodes].ei_edges[1] = ULLONG_MAX;
    ++l_nodes;
  }
  fclose(IN_FA);
  */

  {
    int i;
    for(i=0; i<N_NODES; ++i){
      nodes[i].valid = 0;
      nodes[i].len = 0;
      nodes[i].si_edges[0] = ULLONG_MAX;
      nodes[i].si_edges[1] = ULLONG_MAX;
      nodes[i].si_edges[2] = ULLONG_MAX;
      nodes[i].ei_edges[0] = ULLONG_MAX;
      nodes[i].ei_edges[1] = ULLONG_MAX;
      nodes[i].ei_edges[2] = ULLONG_MAX;
    }
  }

  {
    int i;
    for(i=0; i<N_NODES; ++i){
      blacklist[i] = 0;
    }
  }
  if(f_blist){
    while(fgets(buf,LBUF,IN_BLIST) != NULL){
      if(buf[0] == '#'){
        continue;
      }
      int id;
      int retv;
      retv = sscanf(buf,"%d", &id);
      if(retv != 1){
        fprintf(stderr, "strange blacklist\n");
        exit(1);
      }
      if(id <= N_NODES){
        blacklist[id-1] = 1;
      }
      else{
        fprintf(stderr, "too large id: %d (must be < %d)\n",id,N_NODES);
        exit(1);
      }
    }
    fclose(IN_BLIST);
  }

  char prev_chr[LNAME];
  prev_chr[0] = '\0';
  //unsigned long long max_end=0ull;
  while(gz_mode != 1 && fgets(buf,LBUF,IN_BFMT7) != NULL){
    int line_length;
    line_length = strlen(buf);
    if(line_length >= LBUF-1){
      while(line_length >= LBUF-1 && buf[line_length-1] != '\n'){
        char * newbuf;
        newbuf = (char*)realloc(buf,LBUF*2);
        if(newbuf == NULL){
          fprintf(stderr,"cannot realloc buf\n");
          exit(1);
        }
        buf = newbuf;
        fgets(buf+line_length,LBUF+1,IN_BFMT7);
        line_length = strlen(buf);

        char * newpbuf;
        newpbuf = (char*)realloc(pbuf,LBUF*2);
        if(newpbuf == NULL){
          fprintf(stderr,"cannot realloc pbuf\n");
          exit(1);
        }
        pbuf = newpbuf;

        char * newsbuf;
        newsbuf = (char*)realloc(sbuf,LBUF*2);
        if(newsbuf == NULL){
          fprintf(stderr,"cannot realloc sbuf\n");
          exit(1);
        }
        sbuf = newsbuf;

        LBUF *= 2;
      }
    }
    if(buf[0] == '#'){
      continue;
    }
    int sretval;
    // -outfmt '7 qseSid qstart qend sacc sstart send aliscore evalue pident qseq sseq'
    //
    sretval = sscanf(buf,"%s %d %d %s %d %d %lf %lf %lf %s %s %d %d", t_bl->str_Qid, &t_bl->Qstt, &t_bl->Qend, t_bl->str_Sid, &t_bl->Sstt, &t_bl->Send, &t_bl->aliscore, &t_bl->evalue, &t_bl->pident, sbuf, pbuf, &t_bl->Qlen, &t_bl->Slen);
    //sretval = sscanf(buf,"%s %d %d %s %d %d %lf %lf %lf", t_bl->str_Qid, &t_bl->Sstt, &t_bl->Send, t_bl->str_Sid, &t_bl->Qstt, &t_bl->Qend, &t_bl->aliscore, &t_bl->evalue, &t_bl->pident);
//if(atoi(t_bl->str_Qid) == 362){
//  fprintf(stderr, "362 in0 %d\n",atoi(t_bl->str_Sid));
//}
    if(sretval == 9){
    }
    else if(sretval != 13){
      fprintf(stderr, "sth strange: sretval %d\n",sretval);
      exit(1);
    }
    // TODO
    if(t_bl->aliscore < min_ali_score){
      continue;
    }

    if(prev_chr[0] == '\0'){
      strcpy(prev_chr,t_bl->str_Qid);
    }
    else if(strcmp(prev_chr,t_bl->str_Qid) != 0){
      if(ref_count>0){
        sort_edges(csi_edges, pbuf,sbuf,b_edges_1,b_edges_2,b_edges_3);
      }
      ref_count = 0;

      strcpy(prev_chr, t_bl->str_Qid);
//      l_edges=0;
      csi_edges = l_edges;
      //max_end = 0ull;

      el_index=0ull;
      col_index=0ull;

    }

//if(atoi(t_bl->str_Qid) == 362){
//  fprintf(stderr, "362 in %d\n",atoi(t_bl->str_Sid));
//}

    {
      if(strcmp(t_bl->str_Qid,t_bl->str_Sid) == 0){// this is a ref-ref hit
        ++ref_count;
      }
      else{
        if(ref_count <= 0){
          // mhseol always outputs ref-ref hit
          //continue;// against chunks without ref-ref hits.
        }
      }

      if(1){
        int Sstt,Send,Qstt,Qend;// s: ref, q: query
        Sstt = t_bl->Sstt;
        Send = t_bl->Send;
        Qstt = t_bl->Qstt;
        Qend = t_bl->Qend;
        if(Sstt<1 || Send<1 || Qstt<1 || Qend<1){
          fprintf(stderr,"strange stt/end: %d,%d,%d,%d\n",Sstt,Send,Qstt,Qend);
          exit(1);
        }
        if(Qstt>Qend){
          fprintf(stderr,"strange strand: %d,%d\n",Qstt,Qend);
          exit(1);
        }
        int Qlen,Slen;

        // check if this overlap is of natural type.
        int delta = term_delta; // FIXME

        //strcpy(edges[l_edges].str_Sid,t_bl->str_Sid);
        //strcpy(edges[l_edges].str_Qid,t_bl->str_Qid);
        {
          int pos;
          sscanf(t_bl->str_Sid,"%d%n",&edges[l_edges].Sid,&pos);
          if(pos != strlen(t_bl->str_Sid)){
            fprintf(stderr,"strange Sid: %s\n",t_bl->str_Sid);
            fprintf(stderr,"must be integer > 0\n");
            exit(1);
          }
          sscanf(t_bl->str_Qid,"%d%n",&edges[l_edges].Qid,&pos);
          if(pos != strlen(t_bl->str_Qid)){
            fprintf(stderr,"strange Qid: %s\n",t_bl->str_Qid);
            fprintf(stderr,"must be integer > 0\n");
            exit(1);
          }
        }
        //Slen = nodes[edges[l_edges].Sid-1].len;
        //Qlen = nodes[edges[l_edges].Qid-1].len;
        Qlen = t_bl->Qlen;
        Slen = t_bl->Slen;
        nodes[atoi(t_bl->str_Qid)-1].len = Qlen;
        nodes[atoi(t_bl->str_Sid)-1].len = Slen;
        //if(t_bl->Qlen != Qlen){
          //fprintf(stderr,"strange Qlen: %d %d\n",Qlen,t_bl->Qlen);
        //}
        //if(t_bl->Slen != Slen){
          //fprintf(stderr,"strange Slen: %d %d\n",Slen,t_bl->Slen);
        //}


        edges[l_edges].Sstt=Sstt;
        edges[l_edges].Send=Send;
        edges[l_edges].Qstt=Qstt;
        edges[l_edges].Qend=Qend;
        edges[l_edges].ol_score=t_bl->aliscore;

        //if(t_bl->aliscore < MIN_ALI_SCORE){// XXX
        //  continue;
        //}

        int valid = 0;

        // self hit
        if(strcmp(t_bl->str_Qid,t_bl->str_Sid) == 0){
          continue;
        }

        // contained
        if((Qstt <= delta && Qend > Qlen-delta) || (Qend <= delta && Qstt > Qlen-delta)){
          edges[l_edges].r_or_l_or_a='a';// ambiguous
          valid = 1;
        }
        else if(Sstt < Send){
          if(
            // s: +++>
            // q:   +++>
            Qstt <= delta &&
            Qend < Qlen-delta+1 &&
            delta < Sstt &&
            Slen-delta+1 <= Send
          )
          {
            edges[l_edges].strand='+';
            edges[l_edges].r_or_l_or_a='l';// left
            valid=1;
          }
          else if(
            // q: +++>
            // s:   +++>
            delta < Qstt &&
            Qlen-delta+1 <= Qend &&
            Sstt <= delta &&
            Send < Slen-delta+1
          )
          {
            edges[l_edges].strand='+';
            edges[l_edges].r_or_l_or_a='r';// right
            valid=1;
          }
          else{
            edges[l_edges].r_or_l_or_a='a';// amb
            valid=1;
          }
        }
        else if(Sstt > Send){
          fprintf(stderr,"not come here\n");
          exit(1);
          if(
            // s: <---
            // q:   +++>
            Qstt <= delta &&
            delta >= Send &&
            Qend < Qlen-delta+1 &&
            Slen-delta+1 > Sstt
          )
          {
            edges[l_edges].strand='-';
            edges[l_edges].r_or_l_or_a='l';// left
            valid = 1;
          }
          else if(
            // q: +++>
            // s:   <---
            Sstt >= Slen-delta+1 &&
            Qend >= Qlen-delta+1 &&
            Send > delta &&
            Qstt > delta
          )
          {
            edges[l_edges].strand='-';
            edges[l_edges].r_or_l_or_a='r';// right
              valid=1;
          }
          else{
            edges[l_edges].r_or_l_or_a='a';// amb
            valid=1;
          }
        }
        else{
          fprintf(stderr,"strange Sstt == Send: %d %d\n",Sstt,Send);
          exit(1);
        }

          //fprintf(stderr,"tmp valid: %d\n",valid);
        if(valid){
          if(edges[l_edges].r_or_l_or_a=='a'){
            //fprintf(stderr, "amb: %s -> %s\n",t_bl->str_Qid,t_bl->str_Sid);
          }
        
          nodes[edges[l_edges].Qid-1].valid = 1;
          //nodes[edges[l_edges].Sid-1].valid = 1;
          edges[l_edges].valid = 1;
          edges[l_edges].will_be_deleted = 0;
   //       if(edges[l_edges].Qid == 362){
    //        fprintf(stderr, "362 validated %d\n",edges[l_edges].Sid);
     //     }
    //printf("v/inv %d\n",edges[l_edges].valid);
    //printf("%llu\n",l_edges);
          ++l_edges;
          if(l_nodes < atoi(t_bl->str_Qid)){
            l_nodes = atoi(t_bl->str_Qid);
          }
          if(l_nodes < atoi(t_bl->str_Sid)){
            l_nodes = atoi(t_bl->str_Sid);
          }
        }
        else{
          // do not ++l_edges
        }
      }
    }
  }

  if(ref_count>0){
    sort_edges(csi_edges, pbuf,sbuf,b_edges_1,b_edges_2,b_edges_3);
  }

  ref_count = 0;

  if((gz_mode != 1 && strcmp("-",in_bfmy7) != 0 && fgets(buf,LBUF,IN_BFMT7) != NULL)){
    fprintf(stderr, "not reach the end of %s\n",in_bfmy7);
    exit(1);
  }

  if(f_pgs){
    //print_graph_simple2(OUT_GFA_INIT); //XXX
    //return 0;
  }

  if(opt_freq){
    print_edge_population();
    return 0;
  }

  if(f_pgs){
    print_graph_simple2(OUT_GFA_ORIG);// TODO
    fclose(OUT_GFA_ORIG);
    //return 0;
  }


  if(f_blist)
  {
    int i;
    for(i=0; i<l_nodes; ++i){
      if(blacklist[i]){
        invalidate_node(i);
      }
    }
    delete_edges();
    clean_edges();
  }


  {
    unsigned long long j;
    #pragma omp parallel for
    for(j=0; j<l_nodes; ++j){
      invalidate_ambiguous_edges(j);
    }
  }

  {
    // use a greedy algorithm for repeats
    int i;
    for(i=0; i<l_nodes; ++i){
      int r = get_number_of_io_edges(i,0);
      int l = get_number_of_io_edges(i,1);
      if(r >= repeat_depth && l >= repeat_depth){
        invalidate_node(i);
      }
      else if(r>=repeat_depth){
        reduce_repeat(i,0);
      }
      else if(l>= repeat_depth){
        reduce_repeat(i,1);
      }
      /*
      int d;
      for(d=0; d<=1; ++d){
        if(get_number_of_io_edges(i,d)>=repeat_depth){
          reduce_repeat(i,d);
        }
      }
      */
    }
    delete_edges();
    clean_edges();
  }

  {
    // keep the 5 highest ol_score edges and discard the rest
    int i;
    for(i=0; i<l_nodes; ++i){
      reduce_edges(i,0,MAX_IO_DEGREE);
      reduce_edges(i,1,MAX_IO_DEGREE);
    }
    delete_edges();
    clean_edges();
  }

  {
    // check edges
    unsigned long long i;
    for(i=0; i<l_edges; ++i){
      check_edges(&edges[i]);
    }
    delete_edges();
    clean_edges();
    //return 0;
    fprintf(stderr,"check_edges ok\n");
  }

  if(f_pgs){
    print_graph_simple2(OUT_GFA_INIT);// TODO
    fclose(OUT_GFA_INIT);
    //return 0;
  }

  // backup the graph
  l_nodes_bak = l_nodes;
  l_edges_bak = l_edges;
  copy_nodes(nodes_bak, nodes, l_nodes);
  copy_edges(edges_bak, edges, l_edges);
  {
    int i;
    for(i=0; i<l_nodes; ++i){
      blacklist[i] = 0;
    }
  }
  int prev_l_blacklist = 0;
  int max_iter = 10;
  int iter;

for(iter=0; iter<max_iter; ++iter){

//fprintf(stderr,"iter %d: 7722: %d\n",iter,get_not_small_io_degree(7722));

//fprintf(stderr,"l_nodes: %d\n",l_nodes);

  if(opt_bubblepopping>0)
  {
    fprintf(stderr,"start bubble popping\n");
    // transitive reduction
    unsigned long long j,D;
    //unsigned long long D = 5;
    Dmax = opt_bubblepopping;
    int thread_num;
    int idxstt;

    int k;
    #pragma omp parallel for private(k)
    for(j=0; j<l_nodes; ++j){
      for(k=1; k<=Dmax; ++k){
        loop_remover(j,k,j);
      }
    }
    clean_edges();

    // TODO
    #pragma omp parallel for private(D,thread_num,idxstt)
    for(j=0; j<l_nodes; ++j){
      thread_num = omp_get_thread_num();
      idxstt = NODEBUFMAX*thread_num;
      for(D=2; D<=Dmax; ++D){
        paths[idxstt] = (int)j;
        f_forward(j,0,0,paths+idxstt,D);
      }
    }
    fprintf(stderr,"bp done\n");
    clean_edges();

    #pragma omp parallel for private(k)
    for(j=0; j<l_nodes; ++j){
      for(k=1; k<=6*Dmax; ++k){
        loop_remover(j,k,j);
      }
    }
    clean_edges();

    // TODO
    #pragma omp parallel for private(D,thread_num,idxstt)
    for(j=0; j<l_nodes; ++j){
      thread_num = omp_get_thread_num();
      idxstt = NODEBUFMAX*thread_num;
      for(D=2; D<=6*Dmax; ++D){
        paths[idxstt] = (int)j;
        f_forward(j,0,0,paths+idxstt,D);
      }
    }
    fprintf(stderr,"bp2 done\n");

    clean_edges();

  if(f_pgs){
    print_graph_simple2(OUT_GFA_BP);// TODO
//    fprintf(stderr,"bp done 7722: %d\n",get_not_small_io_degree(7722));
    //return 0;
  }


    // pruning
    int limit = 30;
    int L_IDX_MAX = 100;
    unsigned long long stt,end,foo;
    int s;
    for(s=0; s<=1; ++s){
      #pragma omp parallel for private(stt,end,foo)
      for(j=0; j<l_nodes; ++j){
        stt = nodes[j].si_edges[s];
        end = nodes[j].ei_edges[s];
        if(stt == ULLONG_MAX){// no edge
          continue;
        }
        if(stt == end){// straight
          continue;
        }
        // split
        int max_depth = 0;
        unsigned long long idx_max[100];
        int l_idx = 0;
        int tmp;
        unsigned long long target_idx=ULLONG_MAX;
        for(foo=stt; foo<=end; ++foo){
          if(edges[foo].valid){
            tmp = get_max_depth(edges[foo].Sid-1,s,limit);
            if(tmp < D_LARGE){
              //fprintf(stderr,"d kita %d\n",tmp);
            }
            if(tmp >= D_LARGE){
              if(l_idx < L_IDX_MAX){
                idx_max[l_idx] = foo;
                ++l_idx;
              }
              else{
                fprintf(stderr, "too many edges\n");
                exit(1);
              }
            }
            else if(l_idx == 0 && tmp > max_depth){
              target_idx = foo;
              max_depth = tmp;
            }
          }
        }
        if(l_idx == 0){
          for(foo=stt; foo<=end; ++foo){
            if(edges[foo].valid){
              if(target_idx != ULLONG_MAX && foo != target_idx){
//fprintf(stderr,"kita l_idx=%d, max_depth=%d\n",l_idx,max_depth);
                will_invalidate_this_branch(edges[foo].Sid-1,s,limit);
                edges[foo].will_be_deleted = 1;
              }
            }
          }
        }
        else{
          for(foo=stt; foo<=end; ++foo){
            int flag=0;
            int bar;
            for(bar=0; bar<l_idx; ++bar){
              if(foo == idx_max[bar]){
                flag=1;
                break;
              }
            }
            if(flag == 0){
//fprintf(stderr,"kita itb\n");
              will_invalidate_this_branch(edges[foo].Sid-1,s,limit);
              edges[foo].will_be_deleted = 1;
            }
          }
        }
      }
    }
    delete_edges();
    clean_edges();
    fprintf(stderr,"pruning done\n");
  }

  {
    int i;
    int lbl = 0;
    for(i=0; i<l_nodes; ++i){
      if(get_not_small_io_degree(i)>N_PLOIDY){
        blacklist[i] = 1;
      }
      //if(i+1==7722){
        //fprintf(stderr,"kita not_small_io_degree: %d\n",get_not_small_io_degree(i));
      //}
    }
    for(i=0; i<l_nodes; ++i){
      if(blacklist[i]){
        //invalidate_node(i);
        ++lbl;
      }
    }
    //delete_edges();
    //clean_edges();
    if(lbl == prev_l_blacklist){
fprintf(stderr,"iteration done: %d th\n",iter+1);
      break;
    }
    else{
      prev_l_blacklist = lbl;
    }
    copy_nodes(nodes, nodes_bak, l_nodes_bak);
    copy_edges(edges, edges_bak, l_edges_bak);
    for(i=0; i<l_nodes; ++i){
      if(blacklist[i]){
        invalidate_node(i);
      }
    }
    delete_edges();
    clean_edges();
  }
}
  if(f_pgs){
    //print_graph_simple2(OUT_GFA_PRUNED);
    //return 0;
  }



  {
    int i;
    int f_fork=0;
    for(i=0; i<l_nodes; ++i){
      if(get_not_small_io_degree(i)>N_PLOIDY){
        //nodes[i].valid = 0;
        invalidate_node(i);
        f_fork=1;
        fprintf(stderr,"fork id: %d\n",i+1);
      }
    }
    if(f_fork){
      fprintf(stderr, "WARNING: forking graph. forks were deleted\n");
    }
    delete_edges();
    clean_edges();
  }

  clean_graph();
  if(f_pgs){
    print_graph_simple2(OUT_GFA_PRUNED);// TODO
  }

edge_t * get_highest_score_edge(unsigned long long stt_index, unsigned long long end_index);
void print_overhanging(edge_t * e, int strand);

  // TODO make unitigs from the graph (& retry ol & tr)
  // unitig <=>(def.) left outdegree == 1 && right outdegree == 1.
  // the rest is not unitig
  {
    int i;
    int j;
    // term or not
    for(i=0; i<l_nodes; ++i){
      if(!nodes[i].valid){
        continue;
      }
      if(nodes[i].si_edges[0] == ULLONG_MAX && nodes[i].si_edges[1] != ULLONG_MAX){ // right edge is empty
        nodes[i].valid = 2;
        //fprintf(stderr,"right edge is empty: id=%d\n",edges[nodes[i].si_edges[1]].Qid);
      }
      else if(nodes[i].si_edges[1] == ULLONG_MAX && nodes[i].si_edges[0] != ULLONG_MAX){ // left edge is empty
        nodes[i].valid = 3;
        //fprintf(stderr,"left edge is empty: id=%d\n",edges[nodes[i].si_edges[0]].Qid);
      }
      else if(nodes[i].si_edges[1] != ULLONG_MAX && nodes[i].si_edges[0] != ULLONG_MAX){ // has right and left edges
        nodes[i].valid = 1;
      }
      else if(nodes[i].si_edges[1] == ULLONG_MAX && nodes[i].si_edges[0] == ULLONG_MAX){
        nodes[i].valid = 0;
      }
    }

    // from a term to a term


    int u_id=1;
    for(i=0; i<l_nodes; ++i){
      if(nodes[i].valid == 3){// term. left edge is empty.
        // travers right edges
        for(j=0; j<1; ++j){
          //unsigned long long k;
          unsigned long long stt = nodes[i].si_edges[j];// i_edges[0] is for 'r'
          unsigned long long end = nodes[i].ei_edges[j];
          if(stt == ULLONG_MAX){
            continue;
          }

          int f_stopped = 0;

          edge_t * e = get_highest_score_edge(stt,end);
          if(e == NULL){
            fprintf(stderr,"sth strange. no edge?\n");
            return 1;
          }
          printf(">%d,%c\n",u_id++,'l');// delimiter. linear
          int strand = 1;// ---> +
          //int prev_dir = j;
          //int next_dir = prev_dir;
          int next_dir = j;
          while(nodes[e->Sid-1].valid == 1){// Given nodes without edge were discarded, this is a unitig defined above
            print_overhanging(e, strand);
            e->valid = 0; // used
            invalidate_reverse_edge(e);
            nodes[e->Sid-1].valid = 8+u_id-1;// visited. u_id <- 1, 2, 3...// XXX

            //next_dir = prev_dir;
            if(e->strand == '-'){
              fprintf(stderr,"never come here 100\n");
              exit(1);
              next_dir = 1-next_dir;
              strand = 1-strand;
            }
            //prev_dir = next_dir;

            //e = &edges[nodes[e->Sid-1].si_edges[next_dir]];
            e = get_highest_score_edge(nodes[e->Sid-1].si_edges[next_dir],nodes[e->Sid-1].ei_edges[next_dir]);
            //e = get_highest_score_edge(nodes[e->Sid-1].si_edges[j],nodes[e->Sid-1].ei_edges[j]);
            if(e == NULL){
              //fprintf(stderr,"sth strange 002. no edge?\n");
              fprintf(stderr,"kita1\n");
              printf("#footer\tlinear_end\trepeat_may_be\n");
              f_stopped = 1;
              break;
              //return 1;
            }
//fprintf(stderr,"kita while loop\n");
          }
          if(!f_stopped){
            if(nodes[e->Sid-1].valid == 2){// right edge is empty
              if(e->Sid-1 == i){// circular
                print_overhanging(e, strand);
                printf("#footer\tlinear_end\trepeat_may_be\n");
              }
              else{// linear
                char s;
                if(strand == 1){
                  s = '+';
                }
                else if(strand == 0){
                  s = '-';
                }
                else{
                  fprintf(stderr,"bug 1\n");
                  exit(1);
                }
                // print all the read, not the subsequence. the other terminal
                //printf("%d\t%d\t%d\t%c\n",e->Qid,1,nodes[e->Qid-1].len,s);
                print_overhanging(e, strand);
                printf("%d\t%d\t%d\t%c\n",e->Sid,1,nodes[e->Sid-1].len,s);
                printf("#footer\tlinear_end\n");
              }


              e->valid = 0;
              invalidate_reverse_edge(e);
            }
            else{
              fprintf(stderr,"node's valid: %d\n",nodes[e->Sid-1].valid);
              fprintf(stderr,"kita2\n");
              fprintf(stderr, "bug 2\n");
              exit(1);
            }
          }
        }
      }
    }
//node_checker();
    // the rest edges are in circular chromosome(s)
    for(i=0; i<l_edges; ++i){
      if(!edges[i].valid){
        continue;
      }
      else if(edges[i].r_or_l_or_a != 'r'){
        continue;
      }
      //else if(!nodes[edges[i].Qid-1].valid || nodes[edges[i].Qid-1].valid >= 8){// the node is invalid or already visited
      //  continue;
      //}
      //else if(!nodes[edges[i].Sid-1].valid || nodes[edges[i].Sid-1].valid >= 8){
      //  continue;
      //}
      else if(!nodes[edges[i].Qid-1].valid || !nodes[edges[i].Sid-1].valid || nodes[edges[i].Qid-1].valid >= 8 || nodes[edges[i].Sid-1].valid >= 8){// the node is invalid or already visited
        continue;
      }
      
      edge_t * e = &edges[i];
      int strand=1;
      printf(">%d,%c\n",u_id++,'c');// circular
      //int stt = edges[i].Qid;
      nodes[edges[i].Qid-1].valid = 8+u_id-1;// visited. u_id <- 1, 2, 3...
      int f_stopped = 0;
      //int prev_dir = 0; // right
      //int next_dir = prev_dir;
      int next_dir = 0;// right
      while(nodes[e->Sid-1].valid && nodes[e->Sid-1].valid < 8){ // valid && not visited
        print_overhanging(e, strand);
        e->valid = 0; // used
        invalidate_reverse_edge(e);
        nodes[e->Qid-1].valid = 8+u_id-1;// visited
        if(nodes[e->Sid-1].valid >= 8){
          printf("#footer\tcircular_end\trepeat_may_be\n");
          f_stopped = 1;
          break;
        }

        //next_dir = prev_dir;
        if(e->strand == '-'){
          fprintf(stderr, "never come here\n");
          exit(1);
          next_dir = 1-next_dir;
          strand = 1-strand;
        }
        //prev_dir = next_dir;

        /*
        if(e->r_or_l_or_a == 'r'){
          next_dir = 0;// 'r'
        }
        else if(e->r_or_l_or_a == 'l'){
          next_dir = 1;// 'l'
          //fprintf(stderr, "never come here\n");
          exit(1);
        }
        else{
          fprintf(stderr,"strange r_or_l_or_a: %c\n",e->r_or_l_or_a);
          exit(1);
        }
        if(e->strand == '-'){
          next_dir = 1-next_dir;
          strand = 1-strand;
          fprintf(stderr, "never come here2\n");
          exit(1);
        }
        */
        //fprintf(stderr,"printed Qid %d\n",e->Qid);
        //fprintf(stderr,"printed Sid %d\n",e->Sid);
        //fprintf(stderr,"printed node is valid %d\n",nodes[e->Sid-1].valid);
        //fprintf(stderr,"%llu %llu\n",nodes[e->Sid-1].si_edges[next_dir],nodes[e->Sid-1].ei_edges[next_dir]);
        //node_checker();
        e = get_highest_score_edge(nodes[e->Sid-1].si_edges[next_dir],nodes[e->Sid-1].ei_edges[next_dir]);
        if(e == NULL){
          //fprintf(stderr,"sth strange 003. no edge?| stt %d\n",stt);
          printf("#footer\tcircular_end\trepeat_may_be\n");
          f_stopped = 1;
          break;
          //return 1;
        }
        //e = &edges[nodes[e->Sid-1].si_edges[next_dir]];
        //if(!e->valid){
        //  fprintf(stderr,"strange. does not reached the stt\n");
        //  fprintf(stderr,"# %d %d\n",stt,e->Sid);
        //  exit(1);
        //}
      }
      if(!f_stopped){
        print_overhanging(e, strand);
        printf("#footer\tcircular_end\n");
        e->valid = 0; // used
        invalidate_reverse_edge(e);
      }
      //printf("# %d %d %d\n",stt,e->Qid,e->Sid);
    }
  }

  free(t_bl);
  free(buf);
  free(edges);
  free(s_edges);
  free(b_edges_1);
  free(b_edges_2);
  free(b_edges_3);
  free(edges_bak);
  free(nodes);
  free(nodes_bak);
  free(pbuf);
  free(sbuf);
  free(paths);
  free(blacklist);
  //free(names);

  fclose(OUT_GFA_BP);
  fclose(OUT_GFA_PRUNED);

  return 0;
}

void sort_edges(unsigned long long csi, char * pbuf,char * sbuf,edge_t * b_edges_1, edge_t * b_edges_2, edge_t * b_edges_3){
  // sort by l_or_r stably
  {
    unsigned long long i;
    unsigned long long i_b1=0;
    unsigned long long i_b2=0;
    unsigned long long i_b3=0;
    for(i=csi; i<l_edges; ++i){
      if(edges[i].r_or_l_or_a == 'r'){
        b_edges_1[i_b1++] = edges[i];
      }
      else if(edges[i].r_or_l_or_a == 'l'){
        b_edges_2[i_b2++] = edges[i];
      }
      else if(edges[i].r_or_l_or_a == 'a'){
        b_edges_3[i_b3++] = edges[i];
        // do nothing
      }
      else{
        fprintf(stderr,"unexpected r_or_l_or_a: %c\n",edges[i].r_or_l_or_a);
        exit(1);
      }
    }
    for(i=0; i<i_b1; ++i){
      edges[csi+i] = b_edges_1[i];
    }
    for(i=0; i<i_b2; ++i){
      edges[csi+i_b1+i] = b_edges_2[i];
    }
    for(i=0; i<i_b3; ++i){
      edges[csi+i_b1+i_b2+i] = b_edges_3[i];
    }
    int ti_n = edges[csi].Qid-1;
    nodes[ti_n].si_edges[0] = csi;// 'r'
    nodes[ti_n].si_edges[1] = csi+i_b1;// 'l'
    nodes[ti_n].si_edges[2] = csi+i_b1+i_b2;// 'a'
    nodes[ti_n].ei_edges[0] = csi+i_b1-1;
    nodes[ti_n].ei_edges[1] = csi+i_b1+i_b2-1;
    nodes[ti_n].ei_edges[2] = csi+i_b1+i_b2+i_b3-1;
    if(i_b1 == 0){
      nodes[ti_n].si_edges[0] = ULLONG_MAX;
      nodes[ti_n].ei_edges[0] = ULLONG_MAX;
    }
    if(i_b2 == 0){
      nodes[ti_n].si_edges[1] = ULLONG_MAX;
      nodes[ti_n].ei_edges[1] = ULLONG_MAX;
    }
    if(i_b3 == 0){
      nodes[ti_n].si_edges[2] = ULLONG_MAX;
      nodes[ti_n].ei_edges[2] = ULLONG_MAX;
    }
  }
}

void print_graph_simple(FILE * FH){
  unsigned long long i;
  for(i=0; i<l_edges; ++i){
    if(edges[i].valid){
      if(edges[i].r_or_l_or_a == 'l' && edges[i].strand == '+'){
        fprintf(FH, "L\t%d\t%c\t%d\t%c\t%dM\n",edges[i].Sid,'+',edges[i].Qid,'+',50); // 50, 100, 25 are dummy
        //printf("L\t%d\t%c\t%d\t%c\t%dM\n",edges[i].Qid,'+',edges[i].Sid,'+',50); // 50, 100, 25 are dummy
      }
      else if(edges[i].r_or_l_or_a == 'r' && edges[i].strand == '+'){
        fprintf(FH, "L\t%d\t%c\t%d\t%c\t%dM\n",edges[i].Qid,'+',edges[i].Sid,'+',100);
      }
      // TODO
      /*
      else if(edges[i].r_or_l_or_a == 'l' && edges[i].strand == '-'){
        printf("L\t%d\t%c\t%d\t%c\t%dM\n",edges[i].Qid,'-',edges[i].Sid,'+',50);
      }
      else if(edges[i].r_or_l_or_a == 'r' && edges[i].strand == '-'){
        printf("L\t%d\t%c\t%d\t%c\t%dM\n",edges[i].Qid,'+',edges[i].Sid,'-',100);
      }
      */
      else if(edges[i].r_or_l_or_a == 'a'){
        fprintf(FH, "L\t%d\t%c\t%d\t%c\t%dM\n",edges[i].Qid,'+',edges[i].Sid,'+',25);
      }
      else{
        fprintf(stderr,"strange data 001\n");
      }
    }
  }
}

int edge_sol_comparator(const void * a, const void * b){
  edge_t * ep_a = *(edge_t**)a;
  edge_t * ep_b = *(edge_t**)b;
  //fprintf(stderr,"%d %d\n",ep_a->u_sol,ep_b->u_sol);
  return ep_a->u_sol - ep_b->u_sol;
  //return 1;
}

void print_overhanging(edge_t * e, int strand){
  char s;
  if(strand == 1){
    s = '+';
  }
  else if(strand == 0){
    s = '-';
  }
  else{
    fprintf(stderr,"bug 3\n");
    exit(1);
  }
  int fl_stt;
  int fl_end;
  if(e->Qend > nodes[e->Qid-1].len-term_delta && e->Qstt > term_delta){
    fl_stt = 1;
    fl_end = e->Qstt-1;
  }
  else if(e->Qstt <= term_delta && e->Qend <= nodes[e->Qid-1].len-term_delta){
    fl_stt = e->Qend+1;
    fl_end = nodes[e->Qid-1].len;
  }
  else{
    fprintf(stderr,"bug 4\n");
    fprintf(stderr,"e->Qstt %d\n",e->Qstt);
    fprintf(stderr,"e->Qend %d\n",e->Qend);
    fprintf(stderr,"len-t %d\n",nodes[e->Qid-1].len-term_delta);
    exit(1);
  }
  printf("%d\t%d\t%d\t%c\n",e->Qid,fl_stt,fl_end,s);
}

void invalidate_reverse_edge(edge_t * e){
  unsigned long long tsi = nodes[e->Sid-1].si_edges[0];
  unsigned long long tei = nodes[e->Sid-1].ei_edges[0];
  unsigned long long tsi2 = nodes[e->Sid-1].si_edges[1];
  unsigned long long tei2 = nodes[e->Sid-1].ei_edges[1];
  unsigned long long tsi3 = nodes[e->Sid-1].si_edges[2];
  unsigned long long tei3 = nodes[e->Sid-1].ei_edges[2];
  unsigned long long l;
  int count = 0;
  if(tsi != ULLONG_MAX && tei != ULLONG_MAX){
    for(l = tsi; l<=tei; ++l){
      if(edges[l].Sid == e->Qid){
        edges[l].valid = 0;
        ++count;
      }
    }
  }
  if(tsi2 != ULLONG_MAX && tei2 != ULLONG_MAX){
    for(l = tsi2; l<=tei2; ++l){
      if(edges[l].Sid == e->Qid){
        edges[l].valid = 0;
        ++count;
      }
    }
  }
  if(tsi3 != ULLONG_MAX && tei3 != ULLONG_MAX){
    for(l = tsi3; l<=tei3; ++l){
      if(edges[l].Sid == e->Qid){
        edges[l].valid = 0;
        ++count;
      }
    }
  }
  if(count > 1){
    fprintf(stderr, "unexpected count=%d\n",count);
    exit(1);
  }
}

edge_t * get_highest_score_edge(unsigned long long stt_index, unsigned long long end_index){
  int maxscore = -1;
  unsigned long long e_index = ULLONG_MAX;
  unsigned long long k;
  if(stt_index == ULLONG_MAX){
    return NULL;
  }
  for(k=stt_index; k<=end_index; ++k){
    edge_t * e = &edges[k];
    if(!e->valid){
      continue;
    }
    else if(maxscore < e->ol_score){
      //if(nodes[e->Qid-1].valid && nodes[e->Sid-1].valid && nodes[e->Qid-1].valid < 8 && nodes[e->Sid-1].valid < 8) // valid and not visited
      if(nodes[e->Qid-1].valid && nodes[e->Sid-1].valid){
        maxscore = e->ol_score;
        e_index = k;
      }
    }
  }
  if(e_index != ULLONG_MAX){
    return &edges[e_index];
  }
  else{
    fprintf(stderr, "no edge found\n");
    return NULL;
  }
}

void node_checker(){
  // check edges
  unsigned long long i;
  for(i=0; i<l_nodes; ++i){
    if(!nodes[i].valid){
      continue;
    }
    if(nodes[i].si_edges[0] > nodes[i].ei_edges[0]){
      fprintf(stderr,"strange 1: i= %llu\n",i);
      return;
    }
    if(nodes[i].si_edges[1] > nodes[i].ei_edges[1]){
      fprintf(stderr,"strange 2: i= %llu\n",i);
      return;
    }
  }
}
void f_backward(int qidx, unsigned long long currE, int * paths, unsigned long long Emax1){
  unsigned long long tstt2, tend2;
  tstt2 = nodes[qidx].si_edges[1];
  tend2 = nodes[qidx].ei_edges[1];
  if(tstt2 == ULLONG_MAX){
    return;
  }
  unsigned long long foo,i;
  int next;
  for(foo=tstt2; foo<=tend2; ++foo){
    if(edges[foo].valid){
      next = edges[foo].Sid-1;
    }
    else{
      continue;
    }
    int f_visited = 0;
    for(i=0; i<=Dmax+currE; ++i){
      if(currE < Emax1-1){
        if(next == paths[i]){
          f_visited = 1;
          break;
        }
      }
      else{
        if(next == paths[i] && next != paths[0]){
          f_visited = 1;
          break;
        }
      }
    }
    if(f_visited){
      continue;
    }
    paths[Dmax+currE+1] = next;
    if(currE == Emax1-1){
      if(next == paths[0]){
        edges[foo].valid = 0;
        invalidate_reverse_edge(&edges[foo]);
      }
      if(foo != tend2){
        continue;
      }
      else{
        return;
      }
    }
    f_backward(next,currE+1,paths,Emax1);
  }
}

void f_forward(int qidx, unsigned long long currD, unsigned long long currE, int * paths, unsigned long long Emax2){
  if(Emax2>currD){
    unsigned long long tstt, tend;
    tstt = nodes[qidx].si_edges[0];
    tend = nodes[qidx].ei_edges[0];
    if(tstt == ULLONG_MAX){
      return;
    }
    unsigned long long foo,i;
    int next;
    for(foo=tstt; foo<=tend; ++foo){
      if(edges[foo].valid){
        next = edges[foo].Sid-1;
      }
      else{
        continue;
      }
      int f_visited = 0;
      for(i=0; i<=currD; ++i){
        if(next == paths[i]){
          f_visited = 1;
          break;
        }
      }
      if(f_visited){
        continue;
      }
      paths[currD+1] = next;
      f_forward(next,currD+1,currE,paths,Emax2);
    }
  }
  else{
    unsigned long long E;
    for(E=1; E<=Emax2; ++E){
      f_backward(qidx,0,paths,E);
    }
  }
}

void f_term_or_not(unsigned long long j, int * f_term){
  int s;
  for(s=0; s<=1; ++s){
    unsigned long long stt,end,foo;
    stt = nodes[j].si_edges[s];
    end = nodes[j].ei_edges[s];
    if(stt == ULLONG_MAX){
      *f_term = 1;
      return;
    }
    int tmp=0;
    for(foo=stt; foo<=end; ++foo){
      tmp += edges[foo].valid;
    }
    if(tmp == 0){
      //fprintf(stderr, "kita edge 0\n");
      *f_term = 1;
      return;
    }
    else{
      *f_term = 0;
    }
  }
}

void f_delete_node(unsigned long long j){
  int s;
  for(s=0; s<=1; ++s){
    unsigned long long stt,end,foo;
    stt = nodes[j].si_edges[s];
    end = nodes[j].ei_edges[s];
    if(stt == ULLONG_MAX){
      continue;
    }
    for(foo=stt; foo<=end; ++foo){
      edges[foo].valid = 0;
      invalidate_reverse_edge(&edges[foo]);
    }
  }
  nodes[j].valid = 0;
}

void check_edges(edge_t * e){
  unsigned long long tsi = nodes[e->Sid-1].si_edges[0];
  unsigned long long tei = nodes[e->Sid-1].ei_edges[0];
  unsigned long long tsi2 = nodes[e->Sid-1].si_edges[1];
  unsigned long long tei2 = nodes[e->Sid-1].ei_edges[1];
  unsigned long long tsi3 = nodes[e->Sid-1].si_edges[2];
  unsigned long long tei3 = nodes[e->Sid-1].ei_edges[2];
  unsigned long long l;
  int flag=0;
  if(tsi != ULLONG_MAX && tei != ULLONG_MAX){
    for(l = tsi; l<=tei; ++l){
      if(edges[l].Sid == e->Qid){
        flag=1;
      }
    }
  }
  if(tsi2 != ULLONG_MAX && tei2 != ULLONG_MAX){
    for(l = tsi2; l<=tei2; ++l){
      if(edges[l].Sid == e->Qid){
        flag=1;
      }
    }
  }
  if(tsi3 != ULLONG_MAX && tei3 != ULLONG_MAX){
    for(l = tsi3; l<=tei3; ++l){
      if(edges[l].Sid == e->Qid){
        flag=1;
      }
    }
  }
  if(!flag){
    //fprintf(stderr, "no reverse edge %d - %d (will be removed)\n",e->Qid, e->Sid);
    e->will_be_deleted=1;
    //exit(1);
  }
}

void invalidate_ambiguous_edges(unsigned long long j){
  unsigned long long tsi3 = nodes[j].si_edges[2];
  unsigned long long tei3 = nodes[j].ei_edges[2];
  unsigned long long l;
  if(tsi3 != ULLONG_MAX && tei3 != ULLONG_MAX){
    for(l = tsi3; l<=tei3; ++l){
      edges[l].valid = 0;
      invalidate_reverse_edge(&edges[l]);
    }
  }
}

void swap_edges(unsigned long long a, unsigned long long b){
  if(a == b){
    return;
  }
  edge_t tmp;
  tmp = edges[a];
  edges[a] = edges[b];
  edges[b] = tmp;
}

void clean_edges(){
  unsigned long long j,stt,end,foo,n_valid;
  int s;
  for(j=0; j<l_nodes; ++j){
    for(s=0; s<=1; ++s){
      stt = nodes[j].si_edges[s];
      end = nodes[j].ei_edges[s];
      if(stt == ULLONG_MAX){
        continue;
      }
      //if(nodes[j].valid)
      n_valid = 0;
      for(foo=stt; foo<=end; ++foo){
        if(edges[foo].valid){
          if(foo != stt+n_valid){
            swap_edges(foo,stt+n_valid);
          }
          ++n_valid;
        }
      }
      if(n_valid == 0){
        nodes[j].si_edges[s] = ULLONG_MAX;
        nodes[j].ei_edges[s] = ULLONG_MAX;
      }
      else{
        nodes[j].si_edges[s] = stt;
        nodes[j].ei_edges[s] = stt+n_valid-1;
      }
    }
  }
}


int get_max_depth(unsigned long long idx, int s, int limit){
  if(limit <= 0){
    return D_LARGE;
  }
  unsigned long long stt,end,foo;
  stt = nodes[idx].si_edges[s];
  end = nodes[idx].ei_edges[s];
  if(stt == ULLONG_MAX){
    return 1;
  }
  int max, tmp;
  max = 0;
  for(foo=stt; foo<=end; ++foo){
    if(!edges[foo].valid){
      continue;
    }
    tmp = get_max_depth(edges[foo].Sid-1,s,limit-1);
    if(tmp > max){
      max = tmp;
    }
  }
  return (1+max);
}

int f_merge(unsigned long long idx, int s){
  if(!(s==0 || s==1)){
    fprintf(stderr,"strange direction s %d. s must be 0 or 1\n",s);
    exit(1);
  }
  unsigned long long stt, end;
  stt = nodes[idx].si_edges[1-s];
  end = nodes[idx].ei_edges[1-s];
  if(stt == ULLONG_MAX){
    return 0;// never come here
  }
  else if(stt == end){
    return 0;// not merging
  }
  else{
    return 1; // merging
  }
}

void will_invalidate_this_branch(unsigned long long idx, int s, int limit){
  if(limit <= 0){
    return;
  }
  unsigned long long stt,end,foo;
  stt = nodes[idx].si_edges[s];
  end = nodes[idx].ei_edges[s];
  if(stt == ULLONG_MAX){
    return;
  }
//fprintf(stderr,"limit %d idx %llu\n",limit,idx);
  for(foo=stt; foo<=end; ++foo){
    if(edges[foo].valid){
      will_invalidate_this_branch(edges[foo].Sid-1,s,limit-1);
      edges[foo].will_be_deleted = 1;
    }
  }
}

void delete_edges(){
  unsigned long long j;
  for(j=0; j<l_edges; ++j){
    if(edges[j].will_be_deleted){
      edges[j].valid = 0;
      invalidate_reverse_edge(&edges[j]);
    }
  }
}

int get_degree(int idx){
  int ret=0;
  int s;
  unsigned long long stt,end;
  for(s=0; s<=1; ++s){
    stt = nodes[idx].si_edges[s];
    end = nodes[idx].ei_edges[s];
    if(stt != ULLONG_MAX){
      ret += end-stt+1;
    }
  }
  return ret;
}

void invalidate_node(int idx){
  unsigned long long stt,end,foo;
  int s;
  for(s=0; s<=1; ++s){
    stt = nodes[idx].si_edges[s];
    end = nodes[idx].ei_edges[s];
    if(stt == ULLONG_MAX || end == ULLONG_MAX){
      continue;
    }
    for(foo=stt; foo<=end; ++foo){
      edges[foo].will_be_deleted = 1;
      //edges[foo].valid = 0;
      //invalidate_reverse_edge(&edges[foo]);
    }
  }
  nodes[idx].valid = 0;
}

void copy_nodes(node_t * to, node_t * from, int len){
  int i,j;
  for(i=0; i<len; ++i){
    to[i].valid = from[i].valid;
    to[i].len = from[i].len;
    for(j=0; j<3; ++j){
      to[i].si_edges[j] = from[i].si_edges[j];
      to[i].ei_edges[j] = from[i].ei_edges[j];
    }
    to[i].f_term = from[i].f_term;
  }
}

void copy_edges(edge_t * to, edge_t * from, unsigned long long len){
  unsigned long long i;
  for(i=0; i<len; ++i){
    to[i] = from[i];
  }
}

int get_not_small_io_degree(int idx){
  int iodegrees[2];
  iodegrees[0] = iodegrees[1] = 0;
  int s;
  unsigned long long stt,end;
  for(s=0; s<=1; ++s){
    stt = nodes[idx].si_edges[s];
    end = nodes[idx].ei_edges[s];
    if(stt != ULLONG_MAX && end != ULLONG_MAX){
      iodegrees[s] = end-stt+1;
    }
  }
  //if(idx+1 == 7722){
    //fprintf(stderr,"r,l: %d, %d\n",iodegrees[0],iodegrees[1]);
  //}
  if(iodegrees[0] >= iodegrees[1]){
    return iodegrees[0];
  }
  else{
    return iodegrees[1];
  }
}

void print_graph_simple2(FILE * FH){
  int i;
  int s;
  for(i=0; i<l_nodes; ++i){
    if(!nodes[i].valid){
      continue;
    }
    for(s=0; s<=1; ++s){
      unsigned long long stt,end;
      stt = nodes[i].si_edges[s];
      end = nodes[i].ei_edges[s];
      if(stt == ULLONG_MAX || end == ULLONG_MAX){
        continue;
      }
      unsigned long long foo;
      for(foo=stt; foo<=end; ++foo){
        if(edges[foo].valid){
          if(edges[foo].r_or_l_or_a == 'l' && edges[foo].strand == '+'){
            fprintf(FH, "L\t%d\t%c\t%d\t%c\t%dM\n",edges[foo].Sid,'+',edges[foo].Qid,'+',50); // 50, 100, 25 are dummy
          }
          else if(edges[foo].r_or_l_or_a == 'r' && edges[foo].strand == '+'){
            fprintf(FH, "L\t%d\t%c\t%d\t%c\t%dM\n",edges[foo].Qid,'+',edges[foo].Sid,'+',100);
          }
          else if(edges[foo].r_or_l_or_a == 'a'){
            fprintf(FH, "L\t%d\t%c\t%d\t%c\t%dM\n",edges[foo].Qid,'+',edges[foo].Sid,'+',25);
          }
          else{
            fprintf(stderr,"strange data 001\n");
          }
        }
      }
    }
  }
}

void clean_graph(){
  unsigned long long j,stt,end,foo,n_valid;
  int s;
  for(j=0; j<l_nodes; ++j){
    for(s=0; s<=1; ++s){
      stt = nodes[j].si_edges[s];
      end = nodes[j].ei_edges[s];
      if(stt == ULLONG_MAX){
        continue;
      }
      if(nodes[j].valid){
        n_valid = 0;
        for(foo=stt; foo<=end; ++foo){
          if(edges[foo].valid){
            if(foo != stt+n_valid){
              swap_edges(foo,stt+n_valid);
            }
            ++n_valid;
          }
        }
        if(n_valid == 0){
          nodes[j].si_edges[s] = ULLONG_MAX;
          nodes[j].ei_edges[s] = ULLONG_MAX;
        }
        else{
          nodes[j].si_edges[s] = stt;
          nodes[j].ei_edges[s] = stt+n_valid-1;
        }
      }
      else{
        for(foo=stt; foo<=end; ++foo){
          invalidate_reverse_edge(&edges[foo]);
          edges[foo].valid = 0;
        }
      }
    }
  }
}

int get_number_of_io_edges(int idx, int direction){
  unsigned long long stt,end;
  if(direction < 0 || direction > 3){
    fprintf(stderr, "direction must be 0 or 1 or 2\n");
    return 0;
  }
  stt = nodes[idx].si_edges[direction];
  end = nodes[idx].ei_edges[direction];
  if(stt != ULLONG_MAX && end != ULLONG_MAX){
    return (int)(end-stt+1);
  }
  else{
    return 0;
  }
}
  
void print_edge_population(){
  int i;
  int d;
  for(i=0; i<l_nodes; ++i){
    for(d=0; d<=1; ++d){
      printf("%d\n",get_number_of_io_edges(i,d));
    }
  }
}

void reduce_repeat(int idx, int direction){
  unsigned long long stt,end,foo,MI;
  if(direction < 0 || direction > 3){
    fprintf(stderr, "direction must be 0 or 1 or 2\n");
    return;
  }
  stt = nodes[idx].si_edges[direction];
  end = nodes[idx].ei_edges[direction];
  if(stt == ULLONG_MAX || end == ULLONG_MAX){
    return;
  }
  int max_score = -1;
  MI = ULLONG_MAX;
  for(foo=stt; foo<=end; ++foo){
    if(edges[foo].valid && max_score < edges[foo].ol_score){
      max_score = edges[foo].ol_score;
      MI = foo;
    }
  }
  if(MI == ULLONG_MAX){
    fprintf(stderr, "strange edge scores\n");
    exit(1);
  }
  for(foo=stt; foo<=end; ++foo){
    if(foo != MI){
      invalidate_reverse_edge(&edges[foo]);
      edges[foo].valid = 0;
    }
  }
}

void reduce_edges(unsigned long long idx, int direction, int topx){
  unsigned long long stt,end,foo,i;
  if(direction < 0 || direction > 3){
    fprintf(stderr, "direction must be 0 or 1 or 2\n");
    return;
  }
  stt = nodes[idx].si_edges[direction];
  end = nodes[idx].ei_edges[direction];
  if(stt == ULLONG_MAX || end == ULLONG_MAX){
    return;
  }
  for(i=0; i<topx; ++i){
    if(end>i){
      for(foo=end-i; foo>stt; --foo){
        if(edges[foo].ol_score > edges[foo-1].ol_score){
          swap_edges(foo,foo-1);
          //edge_t tmp = edges[foo-1];
          //edges[foo-1] = edges[foo];
          //edges[foo] = tmp;
        }
      }
    }
  }
  for(foo=stt+topx; foo<=end; ++foo){
    invalidate_reverse_edge(&edges[foo]);
    edges[foo].valid = 0;
  }
}

void loop_remover(unsigned long long idx, int Depth, unsigned long long root){
  unsigned long long stt,end,foo;
  stt = nodes[idx].si_edges[0];
  end = nodes[idx].ei_edges[0];
  if(stt == ULLONG_MAX || end == ULLONG_MAX){
    return;
  }
  for(foo=stt; foo<=end; ++foo){
    if(Depth>0){
      loop_remover(edges[foo].Sid-1,Depth-1,root);
    }
    else{
      if(end != stt){// if the node has two or more right edges
        if(edges[foo].Sid-1==root){
          invalidate_reverse_edge(&edges[foo]);
          edges[foo].valid = 0;
        }
      }
    }
  }
}

