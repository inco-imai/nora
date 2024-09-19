// this program uses Heikki Hyyro's algorithm
// "A Bit-Vector Algorithm for Computing Levenshtein and Damerau Edit Distances"
// The Prague Stringology Conference 2002
// http://www.stringology.org/event/2002/p6.html
//
// Heikki 2005 "Bit-parallel approximate string matching algorithms with transposition"
//
// implemented by Imai

int BW=63;
int DEBUG=0;
int DEBUG2=1;
int DEBUG3=1;

int AW=0;// XXX obsolete
int max_read_length = 60000;
double error_rate = 0.15; // TODO
int l_check = 500; // TODO

//unsigned long long * D0s;
//unsigned long long * HPs;
//unsigned long long * VPs;
//unsigned long long * HNs;
//unsigned long long * VNs;

//int opt_fastq=0;
int opt_print_center=0;
int opt_verbose=0;
int opt_l4=0;
int opt_l3=0;
int opt_l2=0;
int opt_norm=0;
int center=0;

typedef struct v128_t{
  unsigned long long v[2];
}v128_t;

typedef struct v256_t{
  unsigned long long v[4];
}v256_t;

typedef struct v512_t{
  unsigned long long v[8];
}v512_t;

v128_t * D0s_v128_t;
v128_t * HPs_v128_t;
v128_t * VPs_v128_t;
v128_t * HNs_v128_t;
v128_t * VNs_v128_t;

void lshift_v128_t(v128_t * r, v128_t * a, int w){
  // r = a<<w
  if(w<0){
    fprintf(stderr,"WARNING: negative lshift %d\n",w);
    r->v[0] = r->v[1] = 0ull;
    return;
  }
  if(w>=128){
    r->v[0] = r->v[1] = 0ull;
    return;
  }
  if(w>=64){
    r->v[1] = a->v[0] << (w-64);
    r->v[0] = 0ull;
    return;
  }
  if(w==0){
    r->v[0] = a->v[0];
    r->v[1] = a->v[1];
    return;
  }
  // 0<w<=63
  unsigned long long tmp = a->v[0] >> (64-w);
  r->v[0] = a->v[0] << w;
  r->v[1] = (a->v[1] << w) | tmp;
  return;
}

void rshift_v128_t(v128_t * r, v128_t * a, int w){
  // r = a>>w
  if(w<0){
    fprintf(stderr,"WARNING: negative rshift %d\n",w);
    r->v[0] = r->v[1] = 0ull;
    return;
  }
  if(w>=128){
    r->v[1] = r->v[0] = 0ull;
    return;
  }
  if(w>=64){
    r->v[0] = a->v[1] >> (w-64);
    r->v[1] = 0ull;
    return;
  }
  if(w==0){
    r->v[0] = a->v[0];
    r->v[1] = a->v[1];
    return;
  }
  // 0<w<=63
  unsigned long long tmp = a->v[1] << (64-w);
  r->v[1] = a->v[1] >> w;
  r->v[0] = tmp | (a->v[0] >> w);
  return;
}

void add_v128_t(v128_t * r, v128_t * a, v128_t * b){
  // r = a+b
  unsigned long long c = 0ull;
  unsigned long long tmp;
  tmp = a->v[0]+b->v[0];
  if(tmp < a->v[0]){
    //fprintf(stderr,"carry kita\n");
    c = 1ull;
  }
  r->v[0] = tmp;
  r->v[1] = a->v[1]+b->v[1]+c; // ignore a new carry (if exists)
  return;
}

void sub_v128_t(v128_t * r, v128_t * a, v128_t * b){
  // r = a-b
  unsigned long long B = 0ull;
  if(a->v[0] >= b->v[0]){
    r->v[0] = a->v[0]-b->v[0];
  }
  else{
    r->v[0] = a->v[0]+~b->v[0]+1ull;
    B=1ull;
  }
  if(a->v[1] >= b->v[1]){
    r->v[1] = a->v[1]-b->v[1]-B;// ignore a new borrow (if exists)
  }
  else{
    r->v[1] = a->v[1]+~b->v[1]+1ull-B;
    //B=1ull;
  }
  return;
}

void and_v128_t(v128_t * r, v128_t * a, v128_t * b){
  // r = a & b;
  r->v[0] = a->v[0] & b->v[0];
  r->v[1] = a->v[1] & b->v[1];
  return;
}

void xor_v128_t(v128_t * r, v128_t * a, v128_t * b){
  // r = a ^ b;
  r->v[0] = a->v[0] ^ b->v[0];
  r->v[1] = a->v[1] ^ b->v[1];
  return;
}

void or_v128_t(v128_t * r, v128_t * a, v128_t * b){
  // r = a | b;
  r->v[0] = a->v[0] | b->v[0];
  r->v[1] = a->v[1] | b->v[1];
  return;
}

void not_v128_t(v128_t * r, v128_t * a){
  // r = ~a;
  r->v[0] = ~a->v[0];
  r->v[1] = ~a->v[1];
  return;
}

int bool_v128_t(v128_t * k){
  if(k->v[0] | k->v[1]){
    return 1;// true
  }
  return 0;// false
}

void print_str_stderr(char * s, int len){
  int i;
  for(i=0; i<len; ++i){
    fprintf(stderr,"%c",s[i]);
  }
  fprintf(stderr,"\n");
}

void print_rle_stderr(char * buf){
  if(buf == NULL){
    fprintf(stderr," ");
    return;
  }
  int l_r=1;
  int l_b=1;
  if(buf[0] != '\0'){
    while(buf[l_b] != '\0'){
      if(buf[l_b-1] == buf[l_b]){
        ++l_r;
      }
      else{
        fprintf(stderr,"%c%d",buf[l_b-1],l_r);
        l_r = 1;
      }
      ++l_b;
    }
    fprintf(stderr,"%c%d",buf[l_b-1],l_r);
  }
  else{
    fprintf(stderr,"%c%d",'0',0);
  }
  fprintf(stderr," ");
}

void print_3_cont(char * buf){
  int l_b = strlen(buf);
  if(l_b>=256){
    fprintf(stderr,"cannot handle string longer than 256\n");
    exit(1);
  }
  char buf2[256];
  int pop[3];
  char * stts[3];
  stts[0] = &buf2[0];
  {
    int i;
    int i_pop=0;
    int c_pop=1;
    for(i=0; i<l_b-1; ++i){
      if(buf[i+1] != buf[i]){
        pop[i_pop] = c_pop;
        break;
      }
      else{
        ++c_pop;
      }
    }
    //fprintf(stderr,"%d\n",pop[0]);
    //fprintf(stderr,"%d\n",l_b);
    if(i < l_b-1){
      i_pop=2;
      c_pop=1;
      for(i=l_b-1; i>0; --i){
        if(buf[i-1] != buf[i]){
          pop[i_pop] = c_pop;
          break;
        }
        else{
          ++c_pop;
        }
      }
      if(i==0){
        fprintf(stderr,"bug\n");
        exit(1);
      }
      if(pop[2]>=256){
        fprintf(stderr,"bug pop[2] %d\n",pop[2]);
        exit(1);
      }
      //fprintf(stderr,"i=0 %d\n",pop[0]);
      //fprintf(stderr,"i=2 %d\n",pop[2]);
      //fprintf(stderr,"l_b %d\n",l_b);
      pop[1] = l_b-pop[0]-pop[2];
      //fprintf(stderr,"i=1 %d\n",pop[1]);
      int j,k;
      for(j=0; j<3; ++j){
        //fprintf(stderr,"i %d\n",pop[j]);
      }
      for(j=1; j<3; ++j){
        pop[j] += pop[j-1];
      }

      int l_b2=0;
      for(k=0; k<pop[0] && k<256; ++k){
        buf2[l_b2++] = buf[k];
      }
      buf2[l_b2++] = '\0';
      for(j=1; j<3; ++j){
        stts[j] = &buf2[l_b2];
        for(k=pop[j-1]; k<pop[j] && k<256; ++k){
          buf2[l_b2++] = buf[k];
        }
        buf2[l_b2++] = '\0';
      }
      buf2[l_b2-1]='\0';
    }
    else{
      buf[i]='\0';
      strcpy(buf2,buf);
      stts[0] = &buf2[0];
      stts[1] = stts[2] = NULL;
    }
  }
  //fprintf(stderr,"%s\n",stts[1]);
  if(opt_print_center){
    print_rle_stderr(stts[1]);
  }
  else
  {
    int i;
    for(i=0; i<3; ++i){
      //fprintf(stderr,"%s ",stts[i]);
      print_rle_stderr(stts[i]);
    }
  }
  fprintf(stderr,"\n");
  //fprintf(stderr,"%s\n",buf2);
}

void print_bits_stderr_v128_t(v128_t k);

void print_bits_stderr_v128_t_summary(v128_t k){
  if(opt_verbose){
    print_bits_stderr_v128_t(k);
  }
  int i;
  int dig = 2;
  char buf[256];
  //char buf2[256];
  int l_b = 0;
  for(i=dig; i>0; --i)
  {
    unsigned long long mask = 1ull<<63;
    int j;
    unsigned long long num = k.v[i-1];
    for(j=64; j>0; --j){
      if(num & mask){
        buf[l_b] = '+';
        ++l_b;
        //fprintf(stderr,"1");
        //printf("0");
      }
      else{
        buf[l_b] = '-';
        ++l_b;
        //fprintf(stderr,"0");
        //printf("0");
      }
      mask >>= 1;
      if(j%4 == 1){
        //fprintf(stderr," ");
        //printf(" ");
      }
    }
  }
  buf[l_b]='\0';
  print_3_cont(buf);
//  fprintf(stderr,"%s\n",buf);
  //fprintf(stderr,"\n");
  //printf("\n");
  return;
}

void print_bits_stderr_v128_t(v128_t k){
  int i;
  int dig = 2;
  for(i=dig; i>0; --i)
  {
    unsigned long long mask = 1ull<<63;
    int j;
    unsigned long long num = k.v[i-1];
    for(j=64; j>0; --j){
      if(num & mask){
        fprintf(stderr,"1");
        //printf("0");
      }
      else{
        fprintf(stderr,"0");
        //printf("0");
      }
      mask >>= 1;
      if(j%4 == 1){
        fprintf(stderr," ");
        //printf(" ");
      }
    }
  }
  fprintf(stderr,"\n");
  //printf("\n");
  return;
}

void print_bits_stderr_wl_v128_t(v128_t k,char * label){
  fprintf(stderr,"%s: ",label);
  int i;
  int dig = 2;
  for(i=dig; i>0; --i)
  {
    unsigned long long mask = 1ull<<63;
    int j;
    unsigned long long num = k.v[i-1];
    for(j=64; j>0; --j){
      if(num & mask){
        fprintf(stderr,"1");
      }
      else{
        fprintf(stderr,"0");
      }
      mask >>= 1;
      if(j%8 == 1){
        fprintf(stderr," ");
      }
    }
  }
  fprintf(stderr,"\n");
  return;
}

void print_bits(unsigned long long num){
  unsigned long long mask = 1ull<<63;
  {
    int j;
    for(j=64; j>0; --j){
      if(num & mask){
        printf("1");
      }
      else{
        printf("0");
      }
      mask >>= 1;
      if(j%4 == 1){
        printf(" ");
      }
    }
  }
  printf("\n");
  return;
}

void print_bits_stderr(unsigned long long num);

void print_bits_stderr_summary(unsigned long long num){
  if(opt_verbose){
    print_bits_stderr(num);
  }
  char buf[256];
  int l_b=0;
  unsigned long long mask = 1ull<<63;
  {
    int j;
    for(j=64; j>0; --j){
      if(num & mask){
        buf[l_b++] = '+';
        //fprintf(stderr,"1");
      }
      else{
        buf[l_b++] = '-';
        //fprintf(stderr,"0");
      }
      mask >>= 1;
      if(j%4 == 1){
        //fprintf(stderr," ");
      }
    }
  }
  buf[l_b] = '\0';
  //fprintf(stderr,"\n");
  print_3_cont(buf);
  return;
}

void print_bits_stderr(unsigned long long num){
  unsigned long long mask = 1ull<<63;
  {
    int j;
    for(j=64; j>0; --j){
      if(num & mask){
        fprintf(stderr,"1");
      }
      else{
        fprintf(stderr,"0");
      }
      mask >>= 1;
      if(j%4 == 1){
        fprintf(stderr," ");
      }
    }
  }
  fprintf(stderr,"\n");
  return;
}

void print_bits_stderr_wl(unsigned long long num,char * label){
  fprintf(stderr,"%s: ",label);
  unsigned long long mask = 1ull<<63;
  {
    int j;
    for(j=64; j>0; --j){
      if(num & mask){
        fprintf(stderr,"1");
      }
      else{
        fprintf(stderr,"0");
      }
      mask >>= 1;
      if(j%4 == 1){
        fprintf(stderr," ");
      }
    }
  }
  fprintf(stderr,"\n");
  return;
}

int min_of_the_3(int diag, int vert, int hori, int* bt){
  if(diag <= vert && diag<=hori){
    *bt=1;
    return diag;
  }
  else if(vert <= diag && vert<=hori){
    *bt=2;
    return vert;
  }
  else{
    *bt=4;
    return hori;
  }
  return -1;// error
}
void swap(int* a, int* b){
  *a=*a+*b;
  *b=*a-*b;
  *a=*a-*b;
  return;
}

void chomp(char * s){
  int len = strlen(s);
  if(s[len-1] == '\n'){
    s[len-1] = '\0';
  }
}

void normalize_str(char * str){
  int l_str = strlen(str);
  int i;
  for(i=0; i<l_str; ++i){
    if(str[i] == 'A'){
    }
    else if(str[i] == 'C'){
    }
    else if(str[i] == 'G'){
    }
    else if(str[i] == 'T'){
    }
    else if(str[i] == 'a'){
      str[i] = 'A';
    }
    else if(str[i] == 'c'){
      str[i] = 'C';
    }
    else if(str[i] == 'g'){
      str[i] = 'G';
    }
    else if(str[i] == 't'){
      str[i] = 'T';
    }
    else{
      str[i] = 'N';
    }
  }
}

void normalize_str2(char * str, int len){
  int l_str = len;
  int i;
  for(i=0; i<l_str; ++i){
    if(str[i] == 'A'){
    }
    else if(str[i] == 'C'){
    }
    else if(str[i] == 'G'){
    }
    else if(str[i] == 'T'){
    }
    else if(str[i] == 'a'){
      str[i] = 'A';
    }
    else if(str[i] == 'c'){
      str[i] = 'C';
    }
    else if(str[i] == 'g'){
      str[i] = 'G';
    }
    else if(str[i] == 't'){
      str[i] = 'T';
    }
    else{
      str[i] = 'N';
    }
  }
}

void b_heikki_128_re_center(char* txt, char* pat, char * txt_buf, char * pat_buf, char * ali_buf, char * sname, char * qname, int l_txt, int l_pat, int center){

  // w = (1+2*bw)+1;
  int WW = (1+2*BW)+1;
  if(WW > 128){
    //fprintf(stderr, "WARNING: W must be <= 64 in the current implementation.\nWe used W=63\n");
    BW=63;
    WW=(1+2*BW)+1;
    //exit(1);
  }

  //int f_tra=0;

  int orig_center = center;
  if(center < 0){
    //f_tra=1;
    center = -center;
    char * tmp;
    tmp = pat_buf;
    pat_buf = txt_buf;
    txt_buf = tmp;

    tmp = pat;
    pat = txt;
    txt = tmp;

    int tmpi;
    tmpi = l_txt;
    l_txt = l_pat;
    l_pat = tmpi;

    tmp = sname;
    sname = qname;
    qname = tmp;

    /*
    tmpi = t_offset;
    t_offset = p_offset;
    p_offset = tmpi;
    */
    /*
    int tmp = m5buf[l_ah].qid;
    m5buf[l_ah].qid = m5buf[l_ah].sid;
    m5buf[l_ah].sid = tmp;
    tmp = m5buf[l_ah].qlen;
    m5buf[l_ah].qlen = m5buf[l_ah].slen;
    m5buf[l_ah].slen = tmp;
    */
  }
  int i_stt = (center-BW>=0) ? center-BW : 0;

  {
    if(DEBUG2){
      fprintf(stderr,"%s\n",sname);
      print_str_stderr(txt,l_txt);
      //fprintf(stderr,"%s\n",txt);
      fprintf(stderr,"%s\n",qname);
      print_str_stderr(pat,l_pat);
      //fprintf(stderr,"%s\n",pat);
      fprintf(stderr,"l_txt %d\n",l_txt);
      fprintf(stderr,"l_pat %d\n",l_pat);
      fprintf(stderr,"center %d\n",center);
      fprintf(stderr,"pat_sbstr\n    ");
      int i;
      for(i=0; i<WW; ++i){
        fprintf(stderr,"%c",pat[i_stt+i]);
        if(i%8==7){
          fprintf(stderr," ");
        }
      }
      fprintf(stderr,"\n");

      fprintf(stderr,"rev\n    ");
      for(i=0; i<WW; ++i){
        fprintf(stderr,"%c",pat[i_stt+WW-1-i]);
        if(i%8==7){
          fprintf(stderr," ");
        }
      }
      fprintf(stderr,"\n");
    }

    v128_t peq[256];
    {
      int i;
      for(i=0; i<256; i++){
        xor_v128_t(&peq[i],&peq[i],&peq[i]);
        //peq[i] = 0;
      }
    }

    v128_t m_ww;
    v128_t one;
    one.v[0] = 1ull;
    one.v[1] = 0ull;
    lshift_v128_t(&m_ww,&one,WW);
    sub_v128_t(&m_ww,&m_ww,&one);
    //print_bits_stderr_wl_v128_t(m_ww,"t_ww");
    //m_ww = (1ull<<WW)-1ull;

    char ALPHABET[]="ACGT";
    int l_alph = strlen(ALPHABET);
    int ss = (BW-center >= 0) ? BW-center : 0;
    {
      int i;
      v128_t f_mat;
      lshift_v128_t(&f_mat,&one,ss);
      //lshift_v128_t(&f_mat,&one,BW-center);
      //lshift_v128_t(&f_mat,&one,BW);
      //unsigned long long f_mat = 1ull << BW;
      if(DEBUG){
        print_bits_stderr_wl_v128_t(one,"1 ");
        print_bits_stderr_wl_v128_t(f_mat,"fm");
        fprintf(stderr,"ss %d\n",ss);
      }
      for(i=0; i<2*BW+2-ss && i_stt+i<l_pat; ++i){
        or_v128_t(&peq[(int)pat[i_stt+i]],&peq[(int)pat[i_stt+i]],&f_mat);
        //peq[(int)pat[i]] |= f_mat;
        lshift_v128_t(&f_mat,&f_mat,1);
        //f_mat <<= 1;
        /*
        {
          int i;
          if(DEBUG){
            v128_t tmp1,tmp2;
            xor_v128_t(&tmp1,&peq[(int)ALPHABET[0]],&peq[(int)ALPHABET[1]]);
            xor_v128_t(&tmp2,&peq[(int)ALPHABET[2]],&peq[(int)ALPHABET[3]]);
            xor_v128_t(&tmp1,&tmp1,&tmp2);
            for(i=0; i<l_alph; ++i){
              print_bits_stderr_wl_v128_t(peq[(int)ALPHABET[i]],"pq");
            }
            print_bits_stderr_wl_v128_t(tmp1,"xr");
          }
        }
        */
      }
      for(; i<2*BW+2-ss; ++i){
        int j;
        for(j=0; j<l_alph; ++j){
          or_v128_t(&peq[(int)ALPHABET[j]],&peq[(int)ALPHABET[j]],&f_mat);
        }
        lshift_v128_t(&f_mat,&f_mat,1);
      }
    }
    {
      int i;
      if(DEBUG){
        v128_t tmp1,tmp2;
        xor_v128_t(&tmp1,&peq[(int)ALPHABET[0]],&peq[(int)ALPHABET[1]]);
        xor_v128_t(&tmp2,&peq[(int)ALPHABET[2]],&peq[(int)ALPHABET[3]]);
        xor_v128_t(&tmp1,&tmp1,&tmp2);
        for(i=0; i<l_alph; ++i){
          print_bits_stderr_wl_v128_t(peq[(int)ALPHABET[i]],"pq");
        }
        print_bits_stderr_wl_v128_t(tmp1,"xr");
      }
    }

    v128_t vn;
    v128_t vp;
    xor_v128_t(&vn,&vn,&vn);
    xor_v128_t(&vp,&vp,&vp);
    //print_bits_stderr_wl_v128_t(vp,"def. vp");
    {
      // vp = 1^(2*BW+1-ss)0^(1+ss)
      int i;
      v128_t bit;
      lshift_v128_t(&bit,&one,1+ss);
      for(i=0; i<2*BW+1-ss; ++i){
        or_v128_t(&vp,&vp,&bit);
        lshift_v128_t(&bit,&bit,1);
      }
      /*
      int i;
      for(i=0; i<BW+1+center; ++i){
        or_v128_t(&vp,&vp,&one);
        lshift_v128_t(&vp,&vp,1);
        //print_bits_stderr_wl_v128_t(vp,"lsh. vp");
        //vp |= 1ull;
        //vp <<= 1;
      }
      rshift_v128_t(&vp,&vp,1);
      //print_bits_stderr_wl_v128_t(vp,"rsh. vp");
      //vp >>= 1;

      lshift_v128_t(&vp,&vp,WW-(BW+1+center));
      //print_bits_stderr_wl_v128_t(vp,"lsh. vp");
      //fprintf(stderr,"%d\n",WW-(BW+1));
      //vp <<= WW-(BW+1);
      */
    }
    //print_bits_stderr_wl_v128_t(one,"one");
    int ed = 0;
    {
      int j;
      v128_t eq;
      v128_t hp;
      v128_t hn;
      v128_t d0;
      //int l_check = 100;
      int l_mx_stt=0;
      //int l_mx_stt=l_vbuf*omp_get_thread_num();
      int l_mx=l_mx_stt;
      //int l_mx=0;
      //fprintf(stderr,"j-BW %d\n",0-BW);
      //fprintf(stderr,"l_pat %d\n",l_pat);
      for(j=0; j<l_txt && center-BW+j<l_pat; ++j){
        eq = peq[(int)txt[j]];
        //print_bits_stderr_wl_v128_t(vn,"vn");
        //print_bits_stderr_wl_v128_t(vp,"vp");
        if(DEBUG){
          print_bits_stderr_wl_v128_t(eq,"eq");
          print_bits_stderr_wl_v128_t(vp,"vp");
          print_bits_stderr_wl_v128_t(vn,"vn");
        }
        and_v128_t(&d0,&eq,&vp);
        add_v128_t(&d0,&d0,&vp);
        xor_v128_t(&d0,&d0,&vp);
        or_v128_t(&d0,&d0,&eq);
        or_v128_t(&d0,&d0,&vn);
        //and_v128_t(&d0,&d0,&m_ww);
        //d0 = (((eq&vp)+vp)^vp)|eq|vn;
        if(DEBUG){
          print_bits_stderr_wl_v128_t(d0,"d0");
        }

        D0s_v128_t[l_mx] = d0;

        v128_t f_ed;
        lshift_v128_t(&f_ed,&one,BW);
        and_v128_t(&hp,&d0,&f_ed);// hp is just a buffer in this line
        if(!bool_v128_t(&hp))
        //if(!(d0 & (1ull<<BW)))
        {
          ++ed;
        }

        or_v128_t(&hp,&d0,&vp);
        not_v128_t(&hp,&hp);
        or_v128_t(&hp,&hp,&vn);
        //hp = vn | ~( d0 | vp);
        and_v128_t(&hn,&vp,&d0);
        //hn = vp & d0;
        if(DEBUG){
          print_bits_stderr_wl_v128_t(hp,"hp");
          print_bits_stderr_wl_v128_t(hn,"hn");
        }

        HPs_v128_t[l_mx] = hp;

        rshift_v128_t(&d0,&d0,1);
        or_v128_t(&vp,&d0,&hp);
        not_v128_t(&vp,&vp);
        or_v128_t(&vp,&vp,&hn);
        //vp = hn | ~((d0>>1) | hp);
        and_v128_t(&vn,&hp,&d0);
        //vn = hp & (d0>>1);
        if(DEBUG){
          print_bits_stderr_wl_v128_t(vn,"vn");
          print_bits_stderr_wl_v128_t(vp,"vp");
        }

        VPs_v128_t[l_mx] = vp;
        ++l_mx;

        /*
        if((j+1)%l_check == 0){
          if(ed > (j+1)*2.0*error_rate){
            // TODO
            //break;
          }
          else{
            ed = 0;
          }
        }
        */

        {
          int i;
          for(i=0; i<l_alph; ++i){
            rshift_v128_t(&peq[(int)ALPHABET[i]],&peq[(int)ALPHABET[i]],1);
            //peq[(int)ALPHABET[i]] >>= 1;
          }
          int i_new = (i_stt+2*BW+2-ss-1)+(j+1);
          //int i_new = (center+BW+1)+(j+1);
          v128_t lm_bit;
          lshift_v128_t(&lm_bit,&one,WW-1);
          if(i_new < l_pat){
            or_v128_t(&peq[(int)pat[i_new]],&peq[(int)pat[i_new]],&lm_bit);
            //peq[(int)pat[i_new]] |= 1ull<<(WW-1);
          }
          else{
            for(i=0; i<l_alph; ++i){
              or_v128_t(&peq[(int)ALPHABET[i]],&peq[(int)ALPHABET[i]],&lm_bit);
            }
          }
        }
        if(DEBUG){
          fprintf(stderr,"loop %d done\n",j);
        }
      }
      --j;
      // TB
      {
        //fprintf(stderr,"#l_mx %d\n",l_mx);
        int W_in=WW-1;
        int scr=0;
        int rs=0;
        int s_min=max_read_length;
        int l_ali=0;
        int i_ctop=center-BW+j;
        int i_cbottom=center+BW+j;
        int i_min = i_cbottom;
        int j_min = j;
        int d_min = W_in-1;//0-origin
        if(DEBUG2){
          fprintf(stderr,"i_ctop %d\n",i_ctop);
          fprintf(stderr,"i_cbottom %d\n",i_cbottom);
          fprintf(stderr,"i_min %d\n",i_min);
          fprintf(stderr,"j_min %d\n",j_min);
        }
        int vl;
        if(i_ctop>=0){
          vl = W_in-1;// -1: exclude the top of a vector
        }
        else{
          vl = i_cbottom+1-1;
        }
        int k;
        if(DEBUG2){
          fprintf(stderr,"vl %d\n",vl);
          fprintf(stderr,"DEF. scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
        }
        for(k=0; k<vl; ++k){
          v128_t bit,t1,t2;
          lshift_v128_t(&bit,&one,(W_in-1)-k);
          and_v128_t(&t1,&vn,&bit);
          and_v128_t(&t2,&vp,&bit);
          if(bool_v128_t(&t1))
          //if(vn & (1ull << ((W_in-1)-k)))
          {
            ++rs;
          }
          else
          if(bool_v128_t(&t2))
          //if(vp & (1ull << ((W_in-1)-k)))
          {
            --rs;
          }
          scr = rs;
          if(scr <= s_min){
            s_min = scr;
            if(center+BW+j-k > l_pat-1){
              i_min = l_pat-1;
              j_min = j-((center+BW+j-k)-(l_pat-1));
              //fprintf(stderr,"m scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
              //fprintf(stderr,"morning\n");
              //fprintf(stderr,"i_min j_min %d %d\n",i_min,j_min);
            }
            else{
              i_min = center+BW+j-k;
              j_min = j;
              //fprintf(stderr,"e scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
              //fprintf(stderr,"night\n");
              //fprintf(stderr,"i_min j_min %d %d\n",i_min,j_min);
            }
            d_min = (W_in-1)-k;
          }
          //fprintf(stderr,"scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
        }
        l_mx = l_mx_stt+j_min+1;
        if(DEBUG2){
          fprintf(stderr,"i_min %d\n",i_min);
          fprintf(stderr,"j_min %d\n",j_min);
        }

        int i_pat = i_min;
        int i_txt = j_min;
        if(DEBUG2){
          fprintf(stderr,"i_txt %d\n",i_txt);
          fprintf(stderr,"l_txt %d\n",l_txt);
          fprintf(stderr,"i_pat %d\n",i_pat);
          fprintf(stderr,"l_pat %d\n",l_pat);
          fprintf(stderr,"l_mx %d\n",l_mx);
        }
        /*
        pat_buf[l_ali] = pat[i_pat];
        --i_pat;
        txt_buf[l_ali] = txt[i_txt];
        --i_txt;
        if(txt_buf[l_ali] == pat_buf[l_ali]){
          ali_buf[l_ali] = '|';
        }
        else{
          ali_buf[l_ali] = 'x';
        }
        ++l_ali;
        --l_mx;
        */
        /*
        if(f_tra){
          c_m5->qstt = 0;
          c_m5->sstt = i_stt;
          c_m5->qend = 0+j_min+1;
          c_m5->send = i_stt+i_min+1;
        }
        else{
          c_m5->qstt = i_stt;
          c_m5->sstt = 0;
          c_m5->qend = i_stt+i_min+1;
          c_m5->send = 0+j_min+1;
        }
        */

        int score=0;
        int l_match=0;
        int l_mismatch=0;
        int l_ins=0;
        int l_del=0;

        int i_mxs=l_mx-1;
        if(DEBUG2){
          fprintf(stderr,"###\n");
          fprintf(stderr,"i_txt %d\n",i_txt);
          fprintf(stderr,"i_pat %d\n",i_pat);
          fprintf(stderr,"i_mxs %d\n",i_mxs);
          fprintf(stderr,"d_min %d\n",d_min);
        }
        v128_t t1,t2,t3,b1,b2;
        int p_i_mxs = -1;
        int p_i_txt = -1;
        int p_i_pat = -1;
        while(i_mxs>=l_mx_stt && i_txt>=0 && i_pat>=i_stt){
          if(DEBUG2){
            fprintf(stderr,"d_min i_mxs(l_mx_stt) i_txt i_pat %d %d(%d) %d %d\n",d_min,i_mxs,l_mx_stt,i_txt,i_pat);
            if(p_i_mxs == i_mxs && p_i_txt == i_txt && p_i_pat == i_pat){
              fprintf(stderr,"same\n");
              exit(1);
            }
            p_i_mxs = i_mxs;
            p_i_txt = i_txt;
            p_i_pat = i_pat;
          }
          if(d_min <= 0 || d_min >= WW-2){
            fprintf(stderr,"WARNING: out of the band: d_min %d\n",d_min);
            fprintf(stderr,"c\n");
            return;
          }
          if(d_min<0){
            fprintf(stderr,"nega d_min %d\n",d_min);
          }
          if(d_min-1<0){
            fprintf(stderr,"nega d_min-1 %d\n",d_min-1);
          }
          lshift_v128_t(&b1,&one,d_min);
          lshift_v128_t(&b2,&one,d_min-1);
          and_v128_t(&t1,&D0s_v128_t[i_mxs],&b1);
          and_v128_t(&t2,&HPs_v128_t[i_mxs],&b1);
          and_v128_t(&t3,&VPs_v128_t[i_mxs],&b2);
          if(pat[i_pat] == txt[i_txt]){
            // match
            --i_mxs;
            pat_buf[l_ali] = pat[i_pat];
            --i_pat;
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            ali_buf[l_ali] = '|';
            ++l_match;
          }
          else{
            int f_e=0;
            if(i_pat >= i_txt){
              if(bool_v128_t(&t3))
              //else if(VPs[i_mxs] & (1ull << (d_min-1)))
              {
                // upward
                --d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
                ++l_ins;
                ++score;
              }
              else if(bool_v128_t(&t2))
              //else if(HPs[i_mxs] & (1ull << d_min))
              {
                // left
                --i_mxs;
                ++d_min;
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
                ++l_del;
                ++score;
              }
              else if(!bool_v128_t(&t1))
              //else if(!(D0s[i_mxs] & (1ull << d_min)))
              {
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
                ++l_mismatch;
                ++score;
              }
              else{
                f_e = 1;
              }
            }
            else{
              if(bool_v128_t(&t2))
              //else if(HPs[i_mxs] & (1ull << d_min))
              {
                // left
                --i_mxs;
                ++d_min;
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '!';
                ++l_del;
                ++score;
              }
              else if(bool_v128_t(&t3))
              //else if(VPs[i_mxs] & (1ull << (d_min-1)))
              {
                // upward
                --d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
                ++l_ins;
                ++score;
              }
              else if(!bool_v128_t(&t1))
              //else if(!(D0s[i_mxs] & (1ull << d_min)))
              {
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
                ++l_mismatch;
                ++score;
              }
              else{
                f_e = 1;
              }
            }
            if(f_e){
              fprintf(stderr,"not come here 6\n");
              txt_buf[l_ali] = '\0';
              ali_buf[l_ali] = '\0';
              pat_buf[l_ali] = '\0';
              {
                int i;
                int loop = l_ali/2;
                for(i=0; i<loop; ++i){
                  char tmp;
                  tmp = txt_buf[l_ali-1-i];
                  txt_buf[l_ali-1-i] = txt_buf[i];
                  txt_buf[i] = tmp;
                  tmp = ali_buf[l_ali-1-i];
                  ali_buf[l_ali-1-i] = ali_buf[i];
                  ali_buf[i] = tmp;
                  tmp = pat_buf[l_ali-1-i];
                  pat_buf[l_ali-1-i] = pat_buf[i];
                  pat_buf[i] = tmp;
                }
              }
              fprintf(stderr,"%s\n",sname);
              fprintf(stderr,"%s\n",qname);
              fprintf(stderr,"i_mxs %d\n",i_mxs);
              fprintf(stderr,"i_txt %d\n",i_txt);
              fprintf(stderr,"i_pat %d\n",i_pat);
              fprintf(stderr,"%s#\n",txt_buf);
              fprintf(stderr,"%s#\n",ali_buf);
              fprintf(stderr,"%s#\n",pat_buf);
              fprintf(stderr,"##\n");
              v128_t bit;
              v128_t bit2;
              v128_t tmp;
              lshift_v128_t(&bit,&one,d_min);
              lshift_v128_t(&bit2,&one,d_min-1);
              fprintf(stderr,"d_min %d\n",d_min);
              and_v128_t(&tmp,&D0s_v128_t[i_mxs],&bit);
              fprintf(stderr,"D0: %d\n",bool_v128_t(&tmp));
              and_v128_t(&tmp,&HPs_v128_t[i_mxs],&bit);
              fprintf(stderr,"HP: %d\n",bool_v128_t(&tmp));
              //fprintf(stderr,"HP+1: %llu\n",(HPs[i_mxs] & (1ull << (d_min+1))));
              //fprintf(stderr,"HP-1: %llu\n",(HPs[i_mxs] & (1ull << (d_min-1))));
              and_v128_t(&tmp,&VPs_v128_t[i_mxs],&bit2);
              fprintf(stderr,"VP: %d\n",bool_v128_t(&tmp));
              //fprintf(stderr,"VP: %d\n",(VPs[i_mxs] & (1ull << (d_min-1))));
              print_bits_stderr_wl_v128_t(HPs_v128_t[i_mxs],"HP");
              print_bits_stderr_wl_v128_t(VPs_v128_t[i_mxs],"VP");
              print_bits_stderr_wl_v128_t(D0s_v128_t[i_mxs],"D0");
              fprintf(stderr,"##\n");
              fprintf(stderr,">%s\n",sname);
              print_str_stderr(txt,l_txt);
              fprintf(stderr,">%s\n",qname);
              print_str_stderr(pat,l_pat);
              fprintf(stderr,"center %d\n",center);
              fprintf(stderr,"orig_center %d\n",orig_center);
              fprintf(stderr,"txt %c\n",txt[i_txt]);
              fprintf(stderr,"pat %c\n",pat[i_pat]);
              exit(1);
            }
          }
          ++l_ali;
        }
        {
          int i;
          int match = 0;
          for(i=0; i<=l_ali; ++i){
            if(ali_buf[i] == '|'){
              ++match;
            }
          }
          double r = ((double)match/(double)l_ali);
          if(DEBUG2){
            fprintf(stderr,"mrate %f\n",r);
          }
        }
        while(i_pat>=i_stt){
          txt_buf[l_ali] = '-';
          pat_buf[l_ali] = pat[i_pat--];
          ali_buf[l_ali] = '@';
          ++l_ali;
          ++l_ins;
          ++score;
        }
        while(i_txt>=0){
          txt_buf[l_ali] = txt[i_txt--];
          pat_buf[l_ali] = '-';
          ali_buf[l_ali] = '*';
          ++l_ali;
          ++l_del;
          ++score;
        }
        /*
        c_m5->score = score;
        c_m5->l_match = l_match;
        c_m5->l_mismatch = l_mismatch;
        if(f_tra){
          c_m5->l_del = l_ins;
          c_m5->l_ins = l_del;
        }
        else{
          c_m5->l_ins = l_ins;
          c_m5->l_del = l_del;
        }
        */

        txt_buf[l_ali] = '\0';
        ali_buf[l_ali] = '\0';
        pat_buf[l_ali] = '\0';
        if(DEBUG2){
          fprintf(stderr,"i_mxs %d\n",i_mxs);
          fprintf(stderr,"i_txt %d\n",i_txt);
          fprintf(stderr,"i_pat %d\n",i_pat);
          fprintf(stderr,"%s#\n",txt_buf);
          fprintf(stderr,"%s#\n",ali_buf);
          fprintf(stderr,"%s#\n",pat_buf);
          fprintf(stderr,"##\n");
        }
        {
          int i;
          int loop = l_ali/2;
          for(i=0; i<loop; ++i){
            char tmp;
            tmp = txt_buf[l_ali-1-i];
            txt_buf[l_ali-1-i] = txt_buf[i];
            txt_buf[i] = tmp;
            tmp = ali_buf[l_ali-1-i];
            ali_buf[l_ali-1-i] = ali_buf[i];
            ali_buf[i] = tmp;
            tmp = pat_buf[l_ali-1-i];
            pat_buf[l_ali-1-i] = pat_buf[i];
            pat_buf[i] = tmp;
          }
        }
        //fprintf(stderr,"%s\n",txt_buf);
        //fprintf(stderr,"%s\n",ali_buf);
        //fprintf(stderr,"%s\n",pat_buf);
        //fprintf(stderr,"##\n");

      }
    }
    //fprintf(stderr,"@%s\n%s\n+\t%s\n%s\n",&nls[i][1],reads[i],eds[i],eds[i]);
    return;
  }
}
void b_heikki_128_re(char* txt, char* pat, char * txt_buf, char * pat_buf, char * ali_buf, char * name1, char * name2){
  {

    int l_pat = strlen(pat);
    int l_txt = strlen(txt);
    /*
    if(l_txt < l_pat){
      //fprintf(stderr,"your longer string must be longer than your short string.\n");
      char * tmp = txt;
      txt = pat;
      pat = tmp;
      l_pat = strlen(pat);
      l_txt = strlen(txt);
      tmp = name1;
      name1 = name2;
      name2 = name1;
    }
    */
    //if(l_pat > 64){
    //  fprintf(stderr,"pat length must be <= 64\n");
    //  exit(1);
    //}
    normalize_str(pat);
    normalize_str(txt);

    printf("%s\n",name1);
    printf("%s\n",txt);
    printf("%s\n",name2);
    printf("%s\n",pat);
    printf("l_txt %d\n",l_txt);
    printf("l_pat %d\n",l_pat);

    v128_t peq[256];
    {
      int j;
      for(j=0; j<256; j++){
        xor_v128_t(&peq[j],&peq[j],&peq[j]);
        //peq[j] = 0;
      }
    }

    // w = (1+2*bw)+1;
    int WW = (1+2*BW)+1;
    if(WW > 128){
      //fprintf(stderr, "WARNING: W must be <= 64 in the current implementation.\nWe used W=63\n");
      BW=63;
      WW=(1+2*BW)+1;
      //exit(1);
    }
    v128_t m_ww;
    v128_t one;
    one.v[0] = 1ull;
    one.v[1] = 0ull;
    lshift_v128_t(&m_ww,&one,WW);
    sub_v128_t(&m_ww,&m_ww,&one);
    //print_bits_stderr_wl_v128_t(m_ww,"t_ww");
    //m_ww = (1ull<<WW)-1ull;

    char ALPHABET[]="ACGT";
    int l_alph = strlen(ALPHABET);
    {
      int j;
      v128_t f_mat;
      lshift_v128_t(&f_mat,&one,BW);
      //unsigned long long f_mat = 1ull << BW;
      for(j=0; j<BW+2 && j<l_pat; j++){
        or_v128_t(&peq[(int)pat[j]],&peq[(int)pat[j]],&f_mat);
        //peq[(int)pat[j]] |= f_mat;
        lshift_v128_t(&f_mat,&f_mat,1);
        //f_mat <<= 1;
      }
      for(j=l_pat; j<BW+2 ; j++){
        int k;
        for(k=0; k<l_alph; ++k){
          or_v128_t(&peq[(int)ALPHABET[k]],&peq[(int)ALPHABET[k]],&f_mat);
        }
        lshift_v128_t(&f_mat,&f_mat,1);
      }
    }
    {
      int i;
      if(DEBUG){
        for(i=0; i<l_alph; ++i){
          print_bits_stderr_wl_v128_t(peq[(int)ALPHABET[i]],"pq");
        }
      }
    }

    v128_t vn;
    v128_t vp;
    xor_v128_t(&vn,&vn,&vn);
    xor_v128_t(&vp,&vp,&vp);
    //print_bits_stderr_wl_v128_t(vp,"def. vp");
    {
      int i;
      for(i=0; i<BW+1; ++i){
        or_v128_t(&vp,&vp,&one);
        lshift_v128_t(&vp,&vp,1);
        //print_bits_stderr_wl_v128_t(vp,"lsh. vp");
        //vp |= 1ull;
        //vp <<= 1;
      }
      rshift_v128_t(&vp,&vp,1);
      //print_bits_stderr_wl_v128_t(vp,"rsh. vp");
      //vp >>= 1;

      lshift_v128_t(&vp,&vp,WW-(BW+1));
      //print_bits_stderr_wl_v128_t(vp,"lsh. vp");
      //printf("%d\n",WW-(BW+1));
      //vp <<= WW-(BW+1);
    }
    //print_bits_stderr_wl_v128_t(one,"one");
    int ed = 0;
    {
      int j;
      v128_t eq;
      v128_t hp;
      v128_t hn;
      v128_t d0;
      //int l_check = 100;
      int l_mx=0;
      //printf("j-BW %d\n",0-BW);
      //printf("l_pat %d\n",l_pat);
      for(j=0; j<l_txt && j-BW<l_pat; ++j){
        eq = peq[(int)txt[j]];
        //print_bits_stderr_wl_v128_t(vn,"vn");
        //print_bits_stderr_wl_v128_t(vp,"vp");
        if(DEBUG){
          print_bits_stderr_wl_v128_t(eq,"eq");
          print_bits_stderr_wl_v128_t(vp,"vp");
          print_bits_stderr_wl_v128_t(vn,"vn");
        }
        and_v128_t(&d0,&eq,&vp);
        add_v128_t(&d0,&d0,&vp);
        xor_v128_t(&d0,&d0,&vp);
        or_v128_t(&d0,&d0,&eq);
        or_v128_t(&d0,&d0,&vn);
        //and_v128_t(&d0,&d0,&m_ww);
        //d0 = (((eq&vp)+vp)^vp)|eq|vn;
        if(DEBUG){
          print_bits_stderr_wl_v128_t(d0,"d0");
        }

        D0s_v128_t[l_mx] = d0;

        v128_t f_ed;
        lshift_v128_t(&f_ed,&one,BW);
        and_v128_t(&hp,&d0,&f_ed);// hp is just a buffer in this line
        if(!bool_v128_t(&hp))
        //if(!(d0 & (1ull<<BW)))
        {
          ++ed;
        }

        or_v128_t(&hp,&d0,&vp);
        not_v128_t(&hp,&hp);
        or_v128_t(&hp,&hp,&vn);
        //hp = vn | ~( d0 | vp);
        and_v128_t(&hn,&vp,&d0);
        //hn = vp & d0;
        if(DEBUG){
          print_bits_stderr_wl_v128_t(hp,"hp");
          print_bits_stderr_wl_v128_t(hn,"hn");
        }

        HPs_v128_t[l_mx] = hp;

        rshift_v128_t(&d0,&d0,1);
        or_v128_t(&vp,&d0,&hp);
        not_v128_t(&vp,&vp);
        or_v128_t(&vp,&vp,&hn);
        //vp = hn | ~((d0>>1) | hp);
        and_v128_t(&vn,&hp,&d0);
        //vn = hp & (d0>>1);
        if(DEBUG){
          print_bits_stderr_wl_v128_t(vn,"vn");
          print_bits_stderr_wl_v128_t(vp,"vp");
        }

        VPs_v128_t[l_mx] = vp;
        ++l_mx;

        /*
        if((j+1)%l_check == 0){
          if(ed > (j+1)*2.0*error_rate){
            // TODO
            //break;
          }
          else{
            ed = 0;
          }
        }
        */

        {
          int i;
          for(i=0; i<l_alph; ++i){
            rshift_v128_t(&peq[(int)ALPHABET[i]],&peq[(int)ALPHABET[i]],1);
            //peq[(int)ALPHABET[i]] >>= 1;
          }
          int i_new = (BW+1)+(j+1);
          v128_t lm_bit;
          lshift_v128_t(&lm_bit,&one,WW-1);
          if(i_new < l_pat){
            or_v128_t(&peq[(int)pat[i_new]],&peq[(int)pat[i_new]],&lm_bit);
            //peq[(int)pat[i_new]] |= 1ull<<(WW-1);
          }
          else{
            for(i=0; i<l_alph; ++i){
              or_v128_t(&peq[(int)ALPHABET[i]],&peq[(int)ALPHABET[i]],&lm_bit);
            }
          }
        }
      }
      --j;
      // TB
      {
        //printf("#l_mx %d\n",l_mx);
        int W_in=WW-1;
        int scr=0;
        int rs=0;
        int s_min=max_read_length;
        int l_ali=0;
        int i_ctop=BW+j-(W_in-1);
        int i_cbottom=BW+j;
        int i_min = i_cbottom;
        int j_min = j;
        int d_min = W_in-1;//0-origin
        printf("i_ctop %d\n",i_ctop);
        printf("i_cbottom %d\n",i_cbottom);
        printf("i_min %d\n",i_min);
        printf("j_min %d\n",j_min);
        int vl;
        if(i_ctop>=0){
          vl = W_in-1;
        }
        else{
          vl = BW+j+1;
        }
        printf("vl %d\n",vl);
        int k;
        printf("DEF. scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
        for(k=0; k<vl; ++k){
          v128_t bit,t1,t2;
          lshift_v128_t(&bit,&one,(W_in-1)-k);
          and_v128_t(&t1,&vn,&bit);
          and_v128_t(&t2,&vp,&bit);
          if(bool_v128_t(&t1))
          //if(vn & (1ull << ((W_in-1)-k)))
          {
            ++rs;
          }
          else
          if(bool_v128_t(&t2))
          //if(vp & (1ull << ((W_in-1)-k)))
          {
            --rs;
          }

          if(BW+j-k > l_pat-1){
            scr = rs;
            //scr = rs-(BW+j-k-(l_pat-1));// all mismatch
          }
          else{
            scr = rs;
          }
          if(scr <= s_min){
            s_min = scr;
            if(BW+j-k > l_pat-1){
              // TODO
              i_min = l_pat-1;
              j_min = j-((BW+j-k)-(l_pat-1));
              //printf("m scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
              //printf("morning\n");
              //printf("i_min j_min %d %d\n",i_min,j_min);
            }
            else{
              i_min = BW+j-k;
              j_min = j;
              //printf("e scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
              //printf("night\n");
              //printf("i_min j_min %d %d\n",i_min,j_min);
            }
            d_min = (W_in-1)-k;
          }
          //printf("scr rs i_min j_min %d %d %d %d\n",scr,rs,i_min,j_min);
        }
        l_mx = j_min+1;
        printf("i_min %d\n",i_min);
        printf("j_min %d\n",j_min);

        int i_pat = i_min;
        int i_txt = j_min;
        printf("i_txt %d\n",i_txt);
        printf("l_txt %d\n",l_txt);
        printf("i_pat %d\n",i_pat);
        printf("l_pat %d\n",l_pat);
        printf("l_mx %d\n",l_mx);
        /*
        pat_buf[l_ali] = pat[i_pat];
        --i_pat;
        txt_buf[l_ali] = txt[i_txt];
        --i_txt;
        if(txt_buf[l_ali] == pat_buf[l_ali]){
          ali_buf[l_ali] = '|';
        }
        else{
          ali_buf[l_ali] = 'x';
        }
        ++l_ali;
        --l_mx;
        */

        int i_mxs=l_mx-1;
        printf("###\n");
        printf("i_txt %d\n",i_txt);
        printf("i_pat %d\n",i_pat);
        printf("i_mxs %d\n",i_mxs);
        printf("d_min %d\n",d_min);
        v128_t t1,t2,t3,b1,b2;
        while(i_mxs>=0 && i_txt>=0 && i_pat>=0){
          if(d_min <= 0 || d_min >= WW-2){
            fprintf(stdout,"WARNING: out of the band: d_min %d\n",d_min);
            return;
          }
          if(d_min<0){
            fprintf(stderr,"nega d_min %d\n",d_min);
          }
          if(d_min-1<0){
            fprintf(stderr,"nega d_min-1 %d\n",d_min-1);
          }
          lshift_v128_t(&b1,&one,d_min);
          lshift_v128_t(&b2,&one,d_min-1);
          and_v128_t(&t1,&D0s_v128_t[i_mxs],&b1);
          and_v128_t(&t2,&HPs_v128_t[i_mxs],&b1);
          and_v128_t(&t3,&VPs_v128_t[i_mxs],&b2);
          if(pat[i_pat] == txt[i_txt]){
            // match
            --i_mxs;
            pat_buf[l_ali] = pat[i_pat];
            --i_pat;
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            ali_buf[l_ali] = '|';
          }
          else{
            int f_e=0;
            if(i_pat >= i_txt){
              if(bool_v128_t(&t3))
              //else if(VPs[i_mxs] & (1ull << (d_min-1)))
              {
                // upward
                --d_min;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
              }
              else if(bool_v128_t(&t2))
              //else if(HPs[i_mxs] & (1ull << d_min))
              {
                // left
                --i_mxs;
                ++d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
              }
              else if(!bool_v128_t(&t1))
              //else if(!(D0s[i_mxs] & (1ull << d_min)))
              {
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
              }
              else{
                f_e = 1;
              }
            }
            else{
              if(bool_v128_t(&t2))
              //else if(HPs[i_mxs] & (1ull << d_min))
              {
                // left
                --i_mxs;
                ++d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
              }
              else if(bool_v128_t(&t3))
              //else if(VPs[i_mxs] & (1ull << (d_min-1)))
              {
                // upward
                --d_min;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
              }
              else if(!bool_v128_t(&t1))
              //else if(!(D0s[i_mxs] & (1ull << d_min)))
              {
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
              }
              else{
                f_e = 1;
              }
            }
            if(f_e){
              fprintf(stderr,"not come here 6\n");
              txt_buf[l_ali] = '\0';
              ali_buf[l_ali] = '\0';
              pat_buf[l_ali] = '\0';
              {
                int i;
                int loop = l_ali/2;
                for(i=0; i<loop; ++i){
                  char tmp;
                  tmp = txt_buf[l_ali-1-i];
                  txt_buf[l_ali-1-i] = txt_buf[i];
                  txt_buf[i] = tmp;
                  tmp = ali_buf[l_ali-1-i];
                  ali_buf[l_ali-1-i] = ali_buf[i];
                  ali_buf[i] = tmp;
                  tmp = pat_buf[l_ali-1-i];
                  pat_buf[l_ali-1-i] = pat_buf[i];
                  pat_buf[i] = tmp;
                }
              }
              fprintf(stderr,"%s\n",name1);
              fprintf(stderr,"%s\n",name2);
              fprintf(stderr,"i_mxs %d\n",i_mxs);
              fprintf(stderr,"i_txt %d\n",i_txt);
              fprintf(stderr,"i_pat %d\n",i_pat);
              fprintf(stderr,"%s#\n",txt_buf);
              fprintf(stderr,"%s#\n",ali_buf);
              fprintf(stderr,"%s#\n",pat_buf);
              fprintf(stderr,"##\n");
              v128_t bit;
              v128_t bit2;
              v128_t tmp;
              lshift_v128_t(&bit,&one,d_min);
              lshift_v128_t(&bit2,&one,d_min-1);
              and_v128_t(&tmp,&D0s_v128_t[i_mxs],&bit);
              fprintf(stderr,"D0: %d\n",bool_v128_t(&tmp));
              and_v128_t(&tmp,&HPs_v128_t[i_mxs],&bit);
              fprintf(stderr,"HP: %d\n",bool_v128_t(&tmp));
              //fprintf(stderr,"HP+1: %llu\n",(HPs[i_mxs] & (1ull << (d_min+1))));
              //fprintf(stderr,"HP-1: %llu\n",(HPs[i_mxs] & (1ull << (d_min-1))));
              and_v128_t(&tmp,&VPs_v128_t[i_mxs],&bit2);
              fprintf(stderr,"VP: %d\n",bool_v128_t(&tmp));
              //fprintf(stderr,"VP: %d\n",(VPs[i_mxs] & (1ull << (d_min-1))));
              print_bits_stderr_wl_v128_t(HPs_v128_t[i_mxs],"HP");
              print_bits_stderr_wl_v128_t(VPs_v128_t[i_mxs],"VP");
              print_bits_stderr_wl_v128_t(D0s_v128_t[i_mxs],"D0");
              fprintf(stderr,"##\n");
              exit(1);
            }
          }
          ++l_ali;
        }
        {
          int i;
          int match = 0;
          for(i=0; i<=l_ali; ++i){
            if(ali_buf[i] == '|'){
              ++match;
            }
          }
          double r = ((double)match/(double)l_ali);
          fprintf(stdout,"mrate %f\n",r);
        }
        while(i_pat>=0){
          txt_buf[l_ali] = '-';
          pat_buf[l_ali] = pat[i_pat--];
          ali_buf[l_ali] = '@';
          ++l_ali;
        }
        while(i_txt>=0){
          txt_buf[l_ali] = txt[i_txt--];
          pat_buf[l_ali] = '-';
          ali_buf[l_ali] = '*';
          ++l_ali;
        }
        txt_buf[l_ali] = '\0';
        ali_buf[l_ali] = '\0';
        pat_buf[l_ali] = '\0';
        printf("i_mxs %d\n",i_mxs);
        printf("i_txt %d\n",i_txt);
        printf("i_pat %d\n",i_pat);
        printf("%s#\n",txt_buf);
        printf("%s#\n",ali_buf);
        printf("%s#\n",pat_buf);
        printf("##\n");
        {
          int i;
          int loop = l_ali/2;
          for(i=0; i<loop; ++i){
            char tmp;
            tmp = txt_buf[l_ali-1-i];
            txt_buf[l_ali-1-i] = txt_buf[i];
            txt_buf[i] = tmp;
            tmp = ali_buf[l_ali-1-i];
            ali_buf[l_ali-1-i] = ali_buf[i];
            ali_buf[i] = tmp;
            tmp = pat_buf[l_ali-1-i];
            pat_buf[l_ali-1-i] = pat_buf[i];
            pat_buf[i] = tmp;
          }
        }
        //printf("%s\n",txt_buf);
        //printf("%s\n",ali_buf);
        //printf("%s\n",pat_buf);
        //printf("##\n");

      }
    }
    //printf("@%s\n%s\n+\t%s\n%s\n",&nls[i][1],reads[i],eds[i],eds[i]);
    return;
  }
}

void b_heikki_re_core(char* txt, int l_txt, char* pat, int l_pat, char * txt_buf, char * pat_buf, char * ali_buf,
  unsigned long long * D0s,unsigned long long * HPs,unsigned long long * VPs,unsigned long long * HNs,unsigned long long *VNs, int max_read_length, double error_rate, int l_check);
// TODO USED
void b_heikki_re(char* txt, char* pat, char * txt_buf, char * pat_buf, char * ali_buf,
  unsigned long long * D0s,unsigned long long * HPs,unsigned long long * VPs,unsigned long long * HNs,unsigned long long *VNs, int max_read_length, double error_rate, int l_check, int l_chunk){
  {
    //int l_chunk = 5;
    int l_txt = strlen(txt);
    int l_pat = strlen(pat);
    normalize_str2(txt,l_txt);
    normalize_str2(pat,l_pat);

//printf("txt: %s\n",txt);
//printf("pat: %s\n",pat);
    int f_swap = 0;
    if(l_txt < l_pat){
      //fprintf(stderr,"your longer string must be longer than your short string.\n");
      char * tmp = txt;
      txt = pat;
      pat = tmp;
      tmp = txt_buf;
      txt_buf = pat_buf;
      pat_buf = tmp;
      int l_tmp = l_pat;
      l_pat = l_txt;
      l_txt = l_tmp;
      f_swap = 1;
    }

//fprintf(stderr,"b_heikki swap done\n");
//fprintf(stderr,"txt: %s\n",txt);
//fprintf(stderr,"pat: %s\n",pat);
    int loop;
    loop = l_pat/l_chunk;
    int r_txt = l_txt-l_chunk*loop;
    int r_pat = l_pat-l_chunk*loop;

    int off_txt=0;
    int off_pat=0;
    int off_btxt=0;
    int off_bpat=0;
    int off_bali=0;

    int i;
    for(i=0; i<loop; ++i){
      b_heikki_re_core(txt+off_txt, l_chunk, pat+off_pat, l_chunk, txt_buf+off_btxt, pat_buf+off_bpat, ali_buf+off_bali, D0s,HPs,VPs,HNs,VNs, max_read_length, error_rate, l_check);
      if(txt_buf[off_btxt] == '\0'){
        txt_buf[0] = '\0';
        ali_buf[0] = '\0';
        pat_buf[0] = '\0';
        return;
      }
      off_txt += l_chunk;
      off_pat += l_chunk;
      int l_diff=l_chunk;
      while(txt_buf[off_btxt+l_diff] != '\0'){
        ++l_diff;
      }
      /*
      if(l_diff != (int)strlen(txt_buf+off_btxt)){
        fprintf(stderr,"unexpected l_diff:%d, strlen: %d\n",l_diff,(int)strlen(txt_buf+off_btxt));
        exit(1);
      }
      if(l_diff != (int)strlen(ali_buf+off_bali)){
        fprintf(stderr,"unexpected l_diff:%d, strlen: %d\n",l_diff,(int)strlen(ali_buf+off_bali));
        exit(1);
      }
      if(l_diff != (int)strlen(pat_buf+off_bpat)){
        fprintf(stderr,"unexpected l_diff:%d, strlen: %d\n",l_diff,(int)strlen(pat_buf+off_bpat));
        exit(1);
      }
      */
      off_btxt += l_diff;
      off_bali += l_diff;
      off_bpat += l_diff;
      //off_btxt += strlen(txt_buf+off_btxt);
      //off_bali += strlen(ali_buf+off_bali);
      //off_bpat += strlen(pat_buf+off_bpat);
      int f_done=0;
      while(txt_buf[off_btxt-1] == '-'){
        --off_btxt;
        --off_bali;
        --off_bpat;
        --off_pat;
        ++r_pat;
        f_done=1;
      }
      while(!f_done && pat_buf[off_bpat-1] == '-'){
        --off_btxt;
        --off_bali;
        --off_bpat;
        --off_txt;
        ++r_txt;
      }
    }
    if(r_pat <= r_txt && r_pat > 0){
      //b_heikki_re_core(txt+off_txt, r_txt, pat+off_pat, r_pat, txt_buf+off_btxt, pat_buf+off_bpat, ali_buf+off_bali, D0s,HPs,VPs,HNs,VNs, max_read_length, error_rate, l_check);
      b_heikki_re_core(txt+off_txt, r_pat, pat+off_pat, r_pat, txt_buf+off_btxt, pat_buf+off_bpat, ali_buf+off_bali, D0s,HPs,VPs,HNs,VNs, max_read_length, error_rate, l_chunk);
    }
    else if(r_txt < r_pat && r_txt > 0){
      b_heikki_re_core(txt+off_txt, r_txt, pat+off_pat, r_txt, txt_buf+off_btxt, pat_buf+off_bpat, ali_buf+off_bali, D0s,HPs,VPs,HNs,VNs, max_read_length, error_rate, l_chunk);
    }
    if(txt_buf[off_btxt] == '\0'){
      txt_buf[0] = '\0';
      ali_buf[0] = '\0';
      pat_buf[0] = '\0';
      return;
    }

    if(f_swap){
      char * tmp = pat_buf;
      pat_buf = txt_buf;
      txt_buf = tmp;
      tmp = txt;
      txt = pat;
      pat = tmp;
    }
  }
}
void b_heikki_re_core(char* txt, int l_txt, char* pat, int l_pat, char * txt_buf, char * pat_buf, char * ali_buf,
  unsigned long long * D0s,unsigned long long * HPs,unsigned long long * VPs,unsigned long long * HNs,unsigned long long *VNs, int max_read_length, double error_rate, int l_chunk){
  {

    //if(l_pat > 64){
    //  fprintf(stderr,"pat length must be <= 64\n");
    //  exit(1);
    //}

    //normalize_str(pat);
    //normalize_str(txt);

    //printf("%s\n",txt);
    //printf("%s\n",pat);
    //printf("l_txt %d\n",l_txt);
    //printf("l_pat %d\n",l_pat);

    if(l_pat > l_txt){
      fprintf(stderr,"l_pat > l_txt: %d, %d\n",l_pat,l_txt);
      exit(1);
    }
    unsigned long long peq[256];
    {
      int j;
      for(j=0; j<256; j++){
        peq[j] = 0;
      }
    }

    // w = (1+2*bw)+1;
    int WW = (1+2*BW)+1;
    if(WW > 64){
      //fprintf(stderr, "WARNING: W must be <= 64 in the current implementation.\nWe used W=63\n");
      BW=31;
      WW=(1+2*BW)+1;
      //exit(1);
    }
    //unsigned long long m_ww;
    //m_ww = (1ull<<WW)-1ull;

    char ALPHABET[]="ACGT";
    int l_alph = strlen(ALPHABET);
    {
      int j;
      unsigned long long f_mat = 1ull << BW;
      for(j=0; j<BW+2 && j<l_pat; j++){
        peq[(int)pat[j]] |= f_mat;
        f_mat <<= 1;
      }
      for(j=l_pat; j<BW+2; ++j){
        int k;
        for(k=0; k<l_alph; ++k){
          peq[(int)ALPHABET[k]] |= f_mat;
        }
        f_mat <<= 1;
      }
    }
    {
      int i;
      if(DEBUG){
        for(i=0; i<l_alph; ++i){
          print_bits_stderr_wl(peq[(int)ALPHABET[i]],"pq");
        }
      }
    }

    unsigned long long vn= 0ull;
    unsigned long long vp= 0ull;
    {
      int i;
      for(i=0; i<BW+1; ++i){
        vp |= 1ull;
        vp <<= 1;
      }
      vp >>= 1;

      vp <<= WW-(BW+1);
    }
    int ed = 0;
    {
      int j;
      unsigned long long eq;
      unsigned long long hp;
      unsigned long long hn;
      unsigned long long d0;
      int l_mx=0;
      //printf("j-BW %d\n",0-BW);
      //printf("l_pat %d\n",l_pat);
      for(j=0; j<l_txt && j-BW<l_pat; ++j){
        eq = peq[(int)txt[j]];
        if(DEBUG){
          print_bits_stderr_wl(eq,"eq");
          print_bits_stderr_wl(vp,"vp");
          print_bits_stderr_wl(vn,"vn");
        }
        d0 = (((eq&vp)+vp)^vp)|eq|vn;
        if(DEBUG){
          print_bits_stderr_wl(d0,"d0");
        }

        D0s[l_mx] = d0;

        hp = vn | ~( d0 | vp);
        if(hp & (1ull << (BW-1-j))){
          //printf("ok\n");
        }
        else{
          //printf("bad\n");
        }
        hn = vp & d0;
        if(DEBUG){
          print_bits_stderr_wl(hp,"hp");
          print_bits_stderr_wl(hn,"hn");
        }

        HPs[l_mx] = hp;

        vp = hn | ~((d0>>1) | hp);
        vn = hp & (d0>>1);
        if(DEBUG){
          print_bits_stderr_wl(vn,"vn");
          print_bits_stderr_wl(vp,"vp");
        }

        VPs[l_mx] = vp;
        ++l_mx;

        // buggy ? // TODO
        if(!(d0 & (1ull<<BW))){
          ++ed;
        }
        // turned off 2024/08/25 for speed and quantity of output
        // discontinue
        /*
        if((j+1)%l_check == 0){
          //if(ed > (j+1)*1.2*error_rate)
          if((double)ed > ((double)l_check)*error_rate){
            txt_buf[0] = '\0';
            ali_buf[0] = '\0';
            pat_buf[0] = '\0';
            return;
          }
          else{
            ed = 0;
          }
        }
        */

        {
          int i;
          for(i=0; i<l_alph; ++i){
            peq[(int)ALPHABET[i]] >>= 1;
          }
          int i_new = (BW+1)+(j+1);
          if(i_new < l_pat){
            peq[(int)pat[i_new]] |= 1ull<<(WW-1);
          }
          else{
            for(i=0; i<l_alph; ++i){
              peq[(int)ALPHABET[i]] |= 1ull<<(WW-1);
            }
          }
        }
      }
      // check 20240904
      /*
      if(ed > (int)(((double)l_chunk)*error_rate)){
        txt_buf[0] = '\0';
        ali_buf[0] = '\0';
        pat_buf[0] = '\0';
        return;
      }
      */
//fprintf(stderr,"j: %d, l_txt: %d, l_pat: %d\n",j,l_txt,l_pat);
      --j;
      // TB
      {
        int W_in=WW-1;
        int scr=0;
        int rs=0;
        int s_min=max_read_length; // TODO
        int l_ali=0;
        int i_ctop=BW+j-(W_in-1);
        int i_cbottom=BW+j;
        int i_min = i_cbottom;
        int j_min = j;
        int d_min = W_in-1;//0-origin
        //printf("i_ctop %d\n",i_ctop);
        //printf("i_cbottom %d\n",i_cbottom);
        //printf("i_min %d\n",i_min);
        //printf("j_min %d\n",j_min);
        int vl;
        int n_match=0;
        int n_overhanging=0;
        if(i_ctop>=0){
          vl = W_in-1;
        }
        else{
          vl = BW+j+1;
        }
        //printf("vl %d\n",vl);
        int k;
        for(k=0; k<vl; ++k){
          if(vn & (1ull << ((W_in-1)-k))){
            ++rs;
          }
          else if(vp & (1ull << ((W_in-1)-k))){
            --rs;
          }

          if(BW+j-k > l_pat-1){
            scr = rs;
            //scr = rs-(BW+j-k-(l_pat-1));// all mismatch
          }
          else{
            scr = rs;
          }
          if(scr <= s_min){
            s_min = scr;
            if(BW+j-k > l_pat-1){
              //printf("morning\n");
              i_min = l_pat-1;
              j_min = j-((BW+j-k)-(l_pat-1));
            }
            else{
              //printf("night\n");
              i_min = BW+j-k;
              j_min = j;
            }
            d_min = (W_in-1)-k;
          }
          //printf("scr rs %d %d\n",scr,rs);
        }
        l_mx = j_min+1;
        //printf("i_min %d\n",i_min);
        //printf("j_min %d\n",j_min);

        int i_pat = i_min;
        int i_txt = j_min;
        //printf("i_txt %d\n",i_txt);
        //printf("l_txt %d\n",l_txt);
        //printf("i_pat %d\n",i_pat);
        //printf("l_pat %d\n",l_pat);
        //printf("l_mx %d\n",l_mx);
        /*
        pat_buf[l_ali] = pat[i_pat];
        --i_pat;
        txt_buf[l_ali] = txt[i_txt];
        --i_txt;
        if(txt_buf[l_ali] == pat_buf[l_ali]){
          ali_buf[l_ali] = '|';
        }
        else{
          ali_buf[l_ali] = 'x';
        }
        ++l_ali;
        --l_mx;
        */

        int i_mxs=l_mx-1;
        //printf("###\n");
        //printf("i_txt %d\n",i_txt);
        //printf("i_pat %d\n",i_pat);
        //printf("i_mxs %d\n",i_mxs);
        if(l_pat-1-i_pat>0 && l_txt-1-i_txt>0){
          fprintf(stderr,"unexpected\n");
          exit(1);
        }
        else if(l_pat-1-i_pat>0){
          int foo;
          for(foo=0; foo<l_pat-1-i_pat; ++foo){
            pat_buf[l_ali] = pat[l_pat-1-foo];
            ali_buf[l_ali] = '@';
            txt_buf[l_ali] = '-';
            ++l_ali;
            ++n_overhanging;
          }
        }
        else if(l_txt-1-i_txt>0){
          int foo;
          for(foo=0; foo<l_txt-1-i_txt; ++foo){
            txt_buf[l_ali] = txt[l_txt-1-foo];
            ali_buf[l_ali] = '*';
            pat_buf[l_ali] = '-';
            ++l_ali;
            ++n_overhanging;
          }
        }
        while(i_mxs>=0 && i_txt>=0 && i_pat>=0){
          if(d_min <= 0 || d_min >= WW-2)
          //if(d_min < 0 || d_min-1<0 || d_min >= WW-1)
          {
            //fprintf(stderr,"WARNING: out of the band: d_min %d\n",d_min); // XXX
            // MODIFIED XXX
            txt_buf[0] = '\0';
            ali_buf[0] = '\0';
            pat_buf[0] = '\0';
            return;
          }
          if(d_min < 0){
            fprintf(stderr,"negative lshift : d_min %d\n",d_min);
          }
          else if(d_min-1 < 0){
            fprintf(stderr,"negative lshift : d_min-1 %d\n",d_min-1);
          }


            //if(i_pat >= i_txt)
            {
              if(VPs[i_mxs] & (1ull << (d_min-1))){
                // upward
                --d_min;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
              }
              else if(HPs[i_mxs] & (1ull << d_min)){
                // left
                --i_mxs;
                ++d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
              }
            else if(pat[i_pat] != txt[i_txt]){
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
              }
            else{
                // match
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '|';
                ++n_match;
              }
            }
            /*
            else
            {
              if(HPs[i_mxs] & (1ull << d_min)){
                // left
                --i_mxs;
                ++d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
              }
              else if(VPs[i_mxs] & (1ull << (d_min-1))){
                // upward
                --d_min;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
              }
            else if(pat[i_pat] != txt[i_txt]){
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
              }
              else if(!(HNs[i_mxs] & (1ull << d_min))){
                // left
                --i_mxs;
                ++d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
              }
            else if(!(VNs[i_mxs] & (1ull << (d_min-1)))){
                // upward
                --d_min;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
              }
            else{
                // match
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '|';
              }
            }
            */





          /*
          if(pat[i_pat] == txt[i_txt]){
            // match
            --i_mxs;
            pat_buf[l_ali] = pat[i_pat];
            --i_pat;
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            ali_buf[l_ali] = '|';
          }
          else{
            int f_e=0;
            if(i_pat >= i_txt)
            {
              if(VPs[i_mxs] & (1ull << (d_min-1))){
                // upward
                --d_min;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
              }
              else if(HPs[i_mxs] & (1ull << d_min)){
                // left
                --i_mxs;
                ++d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
              }
              else if(!(D0s[i_mxs] & (1ull << d_min)))
              {
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
              }
              else{
                f_e = 1;
              }
            }
            else{
              if(HPs[i_mxs] & (1ull << d_min)){
                // left
                --i_mxs;
                ++d_min;
                if(d_min < 0){
                  fprintf(stderr,"unexpected d_min: %d\n",d_min);
                  exit(1);
                }
                pat_buf[l_ali] = '-';
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = '*';
              }
              else if(VPs[i_mxs] & (1ull << (d_min-1))){
                // upward
                --d_min;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;;
                txt_buf[l_ali] = '-';
                ali_buf[l_ali] = '@';
              }
              else if(!(D0s[i_mxs] & (1ull << d_min)))
              {
                // mismatch
                --i_mxs;
                pat_buf[l_ali] = pat[i_pat];
                --i_pat;
                txt_buf[l_ali] = txt[i_txt];
                --i_txt;
                ali_buf[l_ali] = 'x';
              }
              else{
                f_e = 1;
              }
            }
            if(f_e){
              fprintf(stderr,"not come here 6\n");
              txt_buf[l_ali] = '\0';
              ali_buf[l_ali] = '\0';
              pat_buf[l_ali] = '\0';
              {
                int i;
                int loop = l_ali/2;
                for(i=0; i<loop; ++i){
                  char tmp;
                  tmp = txt_buf[l_ali-1-i];
                  txt_buf[l_ali-1-i] = txt_buf[i];
                  txt_buf[i] = tmp;
                  tmp = ali_buf[l_ali-1-i];
                  ali_buf[l_ali-1-i] = ali_buf[i];
                  ali_buf[i] = tmp;
                  tmp = pat_buf[l_ali-1-i];
                  pat_buf[l_ali-1-i] = pat_buf[i];
                  pat_buf[i] = tmp;
                }
              }
              fprintf(stderr,"i_mxs %d\n",i_mxs);
              fprintf(stderr,"i_txt %d\n",i_txt);
              fprintf(stderr,"i_pat %d\n",i_pat);
              fprintf(stderr,"%s#\n",txt_buf);
              fprintf(stderr,"%s#\n",ali_buf);
              fprintf(stderr,"%s#\n",pat_buf);
              fprintf(stderr,"##\n");
              fprintf(stderr,"D0: %llu\n",(D0s[i_mxs] & (1ull << d_min)));
              fprintf(stderr,"HP: %llu\n",(HPs[i_mxs] & (1ull << d_min)));
              //fprintf(stderr,"HP+1: %llu\n",(HPs[i_mxs] & (1ull << (d_min+1))));
              //fprintf(stderr,"HP-1: %llu\n",(HPs[i_mxs] & (1ull << (d_min-1))));
              fprintf(stderr,"VP: %llu\n",(VPs[i_mxs] & (1ull << (d_min-1))));
              print_bits_stderr_wl(HPs[i_mxs],"HP");
              print_bits_stderr_wl(VPs[i_mxs],"VP");
              print_bits_stderr_wl(D0s[i_mxs],"D0");
              fprintf(stderr,"##\n");
              exit(1);
            }
          }
          */
          ++l_ali;
        }
        while(i_pat>=0){
          txt_buf[l_ali] = '-';
          pat_buf[l_ali] = pat[i_pat--];
          ali_buf[l_ali] = '@';
          ++l_ali;
        }
        while(i_txt>=0){
          txt_buf[l_ali] = txt[i_txt--];
          pat_buf[l_ali] = '-';
          ali_buf[l_ali] = '*';
          ++l_ali;
        }
        txt_buf[l_ali] = '\0';
        ali_buf[l_ali] = '\0';
        pat_buf[l_ali] = '\0';
        if((float)n_match/(float)(l_ali-n_overhanging) < (1.0-error_rate)){
          txt_buf[0] = '\0';
          ali_buf[0] = '\0';
          pat_buf[0] = '\0';
          return;
        }
        //printf("i_mxs %d\n",i_mxs);
        //printf("i_txt %d\n",i_txt);
        //printf("i_pat %d\n",i_pat);
        if(DEBUG){
          fprintf(stderr,"#1#\n");
          fprintf(stderr,"%s#\n",txt_buf);
          fprintf(stderr,"%s#\n",ali_buf);
          fprintf(stderr,"%s#\n",pat_buf);
          fprintf(stderr,"#1#\n");
        }
        {
          int i;
          int loop = l_ali/2;
          for(i=0; i<loop; ++i){
            char tmp;
            tmp = txt_buf[l_ali-1-i];
            txt_buf[l_ali-1-i] = txt_buf[i];
            txt_buf[i] = tmp;
            tmp = ali_buf[l_ali-1-i];
            ali_buf[l_ali-1-i] = ali_buf[i];
            ali_buf[i] = tmp;
            tmp = pat_buf[l_ali-1-i];
            pat_buf[l_ali-1-i] = pat_buf[i];
            pat_buf[i] = tmp;
          }
        }
        if(DEBUG){
          fprintf(stderr,"#2#\n");
          fprintf(stderr,"%s#\n",txt_buf);
          fprintf(stderr,"%s#\n",ali_buf);
          fprintf(stderr,"%s#\n",pat_buf);
          fprintf(stderr,"#2#\n");
        }

      }
    }
    //printf("@%s\n%s\n+\t%s\n%s\n",&nls[i][1],reads[i],eds[i],eds[i]);
    return;
  }
}
void b_heikki_128(char* txt, char* pat, char * txt_buf, char * pat_buf, char * ali_buf){
  {

    int l_pat = strlen(pat);
    int l_txt = strlen(txt);
    if(l_txt < l_pat)
    {
      //fprintf(stderr,"your longer string must be longer than your short string.\n");
      char * tmp = txt;
      txt = pat;
      pat = tmp;
      l_pat = strlen(pat);
      l_txt = strlen(txt);
    }
    normalize_str(pat);
    normalize_str(txt);

    printf("%s\n",txt);
    printf("%s\n",pat);

    v128_t peq[256];
    {
      int j;
      for(j=0; j<256; j++){
        xor_v128_t(&peq[j],&peq[j],&peq[j]);
      }
    }

    // c = bw, w = 1+2*bw;
    // bw > 1;
    int WW = 1+2*BW;
    int CW = BW;
    if(WW > 128){
      fprintf(stderr, "W must be <= 128 in the current implementation.\n");
      exit(1);
    }
    v128_t m_ww;
    v128_t one;
    one.v[0] = 1ull;
    one.v[1] = 0ull;
    {
      lshift_v128_t(&m_ww,&one,WW);
      sub_v128_t(&m_ww,&m_ww,&one);
    }

    v128_t f_mat;
    {
      lshift_v128_t(&f_mat,&one,128-CW);
    }
    char ALPHABET[]="ACGT";
    {
      int j;
      for(j=0; j<CW && j<l_pat; j++){
        or_v128_t(&peq[(int)pat[j]],&peq[(int)pat[j]],&f_mat);
        //peq[(int)pat[j]] |= f_mat;
        lshift_v128_t(&f_mat,&f_mat,1);
        //f_mat <<= 1;
      }
      if(CW > l_pat){
        for(j=l_pat; j<CW; ++j){
          int l_alph = strlen(ALPHABET);
          int i;
          for(i=0; i<l_alph; ++i){
            or_v128_t(&peq[(int)ALPHABET[i]],&peq[(int)ALPHABET[i]],&f_mat);
            //peq[(int)ALPHABET[i]] |= f_mat;
          }
          lshift_v128_t(&f_mat,&f_mat,1);
          //f_mat <<= 1;
        }
      }
    }

    v128_t mv;
    xor_v128_t(&mv,&mv,&mv);
    //unsigned long long mv= 0ull;
    // pv = 1^(c-w')0^(w')1^(w')0^(w-(c+w')
    // // pv = 1^(c+1-w')0^(w')1^(w')0^(w-(c+1+w') <- wrong
    // w' <- AW
    v128_t pv;
    xor_v128_t(&pv,&pv,&pv);
    //unsigned long long pv= 0ull;
    {
      int i;
      for(i=0; i<CW-AW; ++i){
        or_v128_t(&pv,&pv,&one);
        lshift_v128_t(&pv,&pv,1);
        //pv |= 1ull;
        //pv <<= 1;
      }
      rshift_v128_t(&pv,&pv,1);
      //pv >>= 1;
      lshift_v128_t(&pv,&pv,AW);
      /*
      for(i=0; i<AW; ++i){
        lshift_v128_t(&pv,&pv,1);
        //pv <<= 1;
      }
      */
      for(i=0; i<AW; ++i){
        lshift_v128_t(&pv,&pv,1);
        or_v128_t(&pv,&pv,&one);
        //pv <<= 1;
        //pv |= 1ull;
      }
      lshift_v128_t(&pv,&pv,128-(CW+AW));
      /*
      for(i=0; i<WW-(CW+AW); ++i){
        lshift_v128_t(&pv,&pv,1);
        //pv <<= 1;
      }
      */
    }
    //int ed = 0;
    //v128_t f_ed;
    //lshift_v128_t(&f_ed, &one,WW-(CW+1));
    //unsigned long long f_ed = 1ull<<(WW-(CW+1));
    {
      int j;
      v128_t eq;
      v128_t tmp;
      v128_t ph;
      v128_t mh;
      v128_t d0;
      //int l_check = 128;
      int l_mx=0;
      v128_t lm_bit;
      lshift_v128_t(&lm_bit,&one,127);
      for(j=0; j<l_txt; j++){
        {
          int l_alph = strlen(ALPHABET);
          int i;
          for(i=0; i<l_alph; ++i){
            rshift_v128_t(&peq[(int)ALPHABET[i]],&peq[(int)ALPHABET[i]],1);
            //peq[(int)ALPHABET[i]] >>= 1;
          }
          if(j+CW <= l_pat-1){
            or_v128_t(&peq[(int)pat[j+CW]],&peq[(int)pat[j+CW]],&lm_bit);
            //peq[(int)pat[j+CW]] |= (1ull<<(WW-1));
          }
          else if(j+CW-(WW-1) <= l_pat-1){
            for(i=0; i<l_alph; ++i){
              or_v128_t(&peq[(int)ALPHABET[i]],&peq[(int)ALPHABET[i]],&lm_bit);
              //peq[(int)ALPHABET[i]] |= (1ull<<(WW-1));
            }
          }
          else{
            fprintf(stderr,"break\n");
            break;
          }
        }
        eq = peq[(int)txt[j]];
        and_v128_t(&tmp,&eq,&pv);
        add_v128_t(&tmp,&tmp,&pv);
        xor_v128_t(&tmp,&tmp,&pv);
        or_v128_t(&tmp,&tmp,&eq);
        or_v128_t(&d0,&tmp,&mv);
        //tmp = (eq & pv) + pv;
        //tmp = (tmp ^ pv) | eq;
        //d0 = tmp | mv;

        D0s_v128_t[l_mx] = d0;
        //D0s[l_mx] = d0;

        or_v128_t(&tmp,&d0,&pv);
        not_v128_t(&tmp,&tmp);
        or_v128_t(&ph,&mv,&tmp);
        //ph = mv | ~( d0 | pv);
        and_v128_t(&mh,&pv,&d0);
        //mh = pv & d0;

        rshift_v128_t(&d0,&d0,1);
        or_v128_t(&d0,&d0,&lm_bit);
        //d0 >>= 1;
        //d0 |= (1ull<<63);
        or_v128_t(&tmp,&d0,&ph);
        not_v128_t(&tmp,&tmp);
        or_v128_t(&pv,&mh,&tmp);
        //pv = mh | ~(d0 | ph);
        and_v128_t(&mv,&ph,&d0);
        //mv = ph & d0;

        //HPs_v128_t[l_mx] = ph;
        //VPs_v128_t[l_mx] = pv;
        HNs_v128_t[l_mx] = mh;
        VNs_v128_t[l_mx] = mv;
        ++l_mx;

        /*
        and_v128_t(&tmp,&d0,&f_ed);
        if(bool_v128_t(&tmp))
        //if(d0 & f_ed)
        {
        }
        else{
          ++ed;
        }
        //printf("%d\n", ed);
        if((j+1)%l_check == 0){
          if(ed > (j+1)*2.0*error_rate){
            // TODO
            //break;
          }
          else{
            ed = 0;
          }
        }
        */

        /*
        fprintf(stderr,"%d mv pv ph mh d0\n",j);
        print_bits_stderr_v128_t_summary(mv);
        print_bits_stderr_v128_t_summary(pv);
        print_bits_stderr_v128_t_summary(ph);
        print_bits_stderr_v128_t_summary(mh);
        print_bits_stderr_v128_t_summary(d0);
        */
      }
      --j;
      //int j_txt = j;
      //eds[i][j] = '\0';
      // TB
      //printf("\n");
      //print_bits_stderr(pv);
      //print_bits_stderr(mv);
      //print_bits_stderr_v128_t(ph);
      //print_bits_stderr_v128_t_summary(ph);
      //print_bits_stderr_v128_t(pv);
      //print_bits_stderr_v128_t_summary(pv);
      //print_bits_stderr_v128_t(mv);
      //print_bits_stderr_v128_t_summary(mv);
      {
        int i;
        int scr=0;
        int s_min=max_read_length;
        //int i_min=0;
        int i_offset=0;
        int l_ali=0;
        v128_t m_t = lm_bit;
        v128_t tmp;
        v128_t tmp2;
        int i_ctop = CW+j-(WW-1);
        int i_cbottom = CW+j;
        //unsigned long long m_t=1ull;
        printf("WW %d\n",WW);
        int oloop;
        if(i_ctop>=0){
          oloop = WW;
        }
        else{
          oloop = i_cbottom+1;
        }
        for(i=0; i<oloop; ++i){
          and_v128_t(&tmp,&pv,&m_t);
          and_v128_t(&tmp2,&mv,&m_t);
          if(bool_v128_t(&tmp))
          //if(pv & m_t)
          {
            --scr;
            if(scr < s_min){
              s_min = scr;
              i_offset = i;
            }
          }
          else if(bool_v128_t(&tmp2))
          //else if(mv & m_t)
          {
            ++scr;
          }
          rshift_v128_t(&m_t,&m_t,1);
          //m_t >>= 1;
        }
        printf("i_offset %d\n",i_offset);

        int i_txt = j;
        int i_pat = i_cbottom-i_offset;
        //int i_pat = (CW+1)+(j-1)-(WW-1)+i_min;
        printf("i_txt %d\n",i_txt);
        printf("i_pat %d\n",i_pat);
        printf("l_txt %d\n",l_txt);
        printf("l_pat %d\n",l_pat);
        printf("l_mx %d\n",l_mx);
        //int i_pat = j-1+i_min;
        if(i_pat > l_pat-1){
          int diff = i_pat-(l_pat-1);
          //printf("kita\n");
          printf("diff %d\n",diff);
          l_mx -= diff;
          i_txt -= diff;
          i_pat -= diff;
          //i_pat = l_pat-1;
        }
        else{
          //printf("kitenai\n");
        }
        printf("i_txt %d\n",i_txt);
        printf("i_pat %d\n",i_pat);
        printf("l_pat %d\n",l_pat);
        printf("l_mx %d\n",l_mx);
        pat_buf[l_ali] = pat[i_pat];
        --i_pat;
        txt_buf[l_ali] = txt[i_txt];
        --i_txt;
        if(txt_buf[l_ali] == pat_buf[l_ali]){
          ali_buf[l_ali] = '|';
        }
        else{
          ali_buf[l_ali] = 'x';
        }
        ++l_ali;

        /*
        if(i_pat > l_pat-1){
          pat_buf[l_ali] = ' ';
        }
        else{
          pat_buf[l_ali] = pat[i_pat];
        }
        --i_pat;

        txt_buf[l_ali] = txt[i_txt];
        --i_txt;
        if(txt_buf[l_ali] == pat_buf[l_ali]){
          ali_buf[l_ali] = '|';
        }
        else if(pat_buf[l_ali] == ' '){
          ali_buf[l_ali] = ' ';
        }
        else{
          ali_buf[l_ali] = 'x';
        }
        ++l_ali;
        */

        int i_mxs=l_mx-1;
        while(i_mxs>=0 && i_txt>=0 && i_pat>=0){
          //if(pat[i_pat] == txt[i_txt] || ((VPs[i_mxs] ^ HPs[i_mxs]) & ~(VNs[i_mxs] | HNs[i_mxs])) & (1ull << i_min))
          v128_t focus;
          v128_t vfocus;
          v128_t tmp;
          v128_t tmp2;
          v128_t tmp3;
          lshift_v128_t(&focus,&one,127-i_offset);
          lshift_v128_t(&vfocus,&one,127-i_offset-1);
          and_v128_t(&tmp,&D0s_v128_t[i_mxs],&focus);
          and_v128_t(&tmp2,&HPs_v128_t[i_mxs],&focus);
          and_v128_t(&tmp3,&VPs_v128_t[i_mxs],&vfocus);
          if(pat[i_pat] == txt[i_txt] || !bool_v128_t(&tmp))
          //if(pat[i_pat] == txt[i_txt] || !(D0s[i_mxs] & (1ull << i_min)))
          //if(D0s[i_mxs] & (1ull << i_min) && pat[i_pat] == txt[i_txt])
          {
            --i_mxs;
            if(i_pat > l_pat-1){
              pat_buf[l_ali] = ' ';
              fprintf(stderr,"not come here\n");
              exit(1);
            }
            else{
              pat_buf[l_ali] = pat[i_pat];
            }
            --i_pat;
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            if(txt_buf[l_ali] == pat_buf[l_ali]){
              ali_buf[l_ali] = '|';
            }
            else if(pat_buf[l_ali] == ' '){
              ali_buf[l_ali] = ' ';
            }
            else{
              ali_buf[l_ali] = 'x';
            }
          }
          else if(bool_v128_t(&tmp2))
          //else if(HPs[i_mxs] & (1ull << i_min))
          {
            --i_mxs;
            --i_offset;
            if(i_offset){
              fprintf(stderr,"unexpected i_offset: %d\n",i_offset);
              exit(1);
            }
            pat_buf[l_ali] = '-';
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            ali_buf[l_ali] = '*';
          }
          else if(bool_v128_t(&tmp3))
          //else if(VPs[i_mxs] & (1ull << i_min))
          {
            ++i_offset;
            if(i_pat > l_pat-1){
              pat_buf[l_ali] = ' ';
              fprintf(stderr,"not come here\n");
              exit(1);
            }
            else{
              pat_buf[l_ali] = pat[i_pat];
            }
            --i_pat;;
            txt_buf[l_ali] = '-';
            ali_buf[l_ali] = '@';
          }
          else{
            --i_mxs;
            if(i_pat > l_pat-1){
              pat_buf[l_ali] = ' ';
              fprintf(stderr,"not come here\n");
              exit(1);
            }
            else{
              pat_buf[l_ali] = pat[i_pat];
            }
            --i_pat;
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            if(txt_buf[l_ali] == pat_buf[l_ali]){
              ali_buf[l_ali] = '|';
            }
            else if(pat_buf[l_ali] == ' '){
              ali_buf[l_ali] = ' ';
            }
            else{
              ali_buf[l_ali] = 'x';
            }
          }
          ++l_ali;
        }
        txt_buf[l_ali] = '\0';
        ali_buf[l_ali] = '\0';
        pat_buf[l_ali] = '\0';
        printf("i_mxs %d\n",i_mxs);
        printf("i_txt %d\n",i_txt);
        printf("i_pat %d\n",i_pat);
        printf("%s#\n",txt_buf);
        printf("%s#\n",ali_buf);
        printf("%s#\n",pat_buf);
        printf("##\n");
        {
          int i;
          int loop = l_ali/2;
          for(i=0; i<loop; ++i){
            char tmp;
            tmp = txt_buf[l_ali-1-i];
            txt_buf[l_ali-1-i] = txt_buf[i];
            txt_buf[i] = tmp;
            tmp = ali_buf[l_ali-1-i];
            ali_buf[l_ali-1-i] = ali_buf[i];
            ali_buf[i] = tmp;
            tmp = pat_buf[l_ali-1-i];
            pat_buf[l_ali-1-i] = pat_buf[i];
            pat_buf[i] = tmp;
          }
        }
        //printf("%s\n",txt_buf);
        //printf("%s\n",ali_buf);
        //printf("%s\n",pat_buf);
        //printf("##\n");

      }
      /*
      int i_stt;
      int j_stt;
      int i_mxs;
      int i_hmin;
      {
        int i_top = CW+j_txt-(WW-1);
        fprintf(stderr,"i_top %d\n",i_top);
        int i_bottom = CW+j_txt;
        if(i_top >= l_pat-1)
        {
          // `  `  
          //  `__` 
          //
          i_stt = i_top;
          if(i_stt > l_pat-1){
            fprintf(stderr,"%d %d\n",i_stt,l_pat-1);
            fprintf(stderr,"unexpected i_stt %d (> %d)\n",i_stt,l_pat-1);
            exit(1);
          }
          int i_to = j_txt-(WW-1);
          if(i_to < 0){
            i_to = 0;
          }
          v128_t t1;
          v128_t t2;
          int scr=0;
          int m_sc=scr;
          int m_idx=j_txt;
          int i;
          fprintf(stderr,"i_to %d\n",i_to);
          fprintf(stderr,"j_txt %d\n",j_txt);
          for(i=j_txt; i>=i_to; --i){
            if(j_txt-i<0){
              fprintf(stderr,"unexpected j_txt-i: %d\n",j_txt-i);
              exit(1);
            }
            rshift_v128_t(&t1,&HPs_v128_t[i],j_txt-i);
            and_v128_t(&t1,&t1,&one);
            rshift_v128_t(&t2,&HNs_v128_t[i],j_txt-i);
            //print_bits_stderr_v128_t(t2);
            and_v128_t(&t2,&t2,&one);
            if(bool_v128_t(&t1)){
              --scr;
              if(scr < m_sc){
                m_idx = i-1;
                m_sc = scr;
                i_hmin = j_txt-i;
              }
            }
            else if(bool_v128_t(&t2)){
              ++scr;
            }
            else{
              // +-0
            }
          }
          j_stt = m_idx;
          i_mxs = m_idx;
          fprintf(stderr,"j_stt m_idx %d\n",j_stt);
        }
        else if(i_bottom <= l_pat-1){
          // `  `  
          //  `  ` 
          //   ` | 
          //    `| 
          //
          j_stt = j_txt;
          i_mxs = j_txt;
          int i;
          v128_t t1;
          v128_t t2;
          int scr=0;
          int m_sc=scr;
          int m_idx=0;
          for(i=0; i<WW; ++i){
            rshift_v128_t(&t1,&VNs_v128_t[j_txt],i);
            and_v128_t(&t1,&t1,&one);
            rshift_v128_t(&t2,&VPs_v128_t[j_txt],i);
            and_v128_t(&t2,&t2,&one);
            if(bool_v128_t(&t1)){
              --scr;
              if(scr < m_sc){
                m_idx=i+1;
                m_sc = scr;
                i_hmin = i;
              }
            }
            else if(bool_v128_t(&t2)){
              ++scr;
            }
            else{
            }
          }
          i_stt = i_top+m_idx;
          fprintf(stderr,"m_idx i_stt %d %d\n",m_idx,i_stt);
          //i_stt = CW+j_txt-(WW-1)+m_idx;
        }
        else{
          // `  `  
          //  `  ` 
          //   `_| 
          //
          int i_to = (l_pat-1)-CW;
          if(i_to<0){
            i_to=0;
          }
          int i;
          int offset;
          v128_t t1;
          v128_t t2;
          int o_scr=0;
          int scr=o_scr;
          int m_hsc=o_scr;
          int m_hidx=j_txt;
          int m_vsc=o_scr;
          int vcells = WW-(i_bottom-(l_pat-1));
          int m_vidx=vcells-1;
          if(vcells<1){
            fprintf(stderr,"unexpected vcells: %d\n",vcells);
            exit(1);
          }
          int t_ih1;
          int t_ih2;
          for(i=j_txt; i>=i_to; --i){
            offset = vcells-1;
            if(j_txt-i<0){
              fprintf(stderr,"uNexpected j_txt-i: %d\n",j_txt-i);
              exit(1);
            }
            rshift_v128_t(&t1,&HPs_v128_t[i],offset+j_txt-i);
            and_v128_t(&t1,&t1,&one);
            rshift_v128_t(&t2,&HNs_v128_t[i],offset+j_txt-i);
            and_v128_t(&t2,&t2,&one);
            if(bool_v128_t(&t1)){
              --scr;
              if(scr < m_hsc){
                m_hidx = i-1;
                m_hsc = scr;
                t_ih1 = offset+j_txt-i;
              }
            }
            else if(bool_v128_t(&t2)){
              ++scr;
            }
            else{
              // +-0
            }
          }
          scr=o_scr;
          for(i=vcells-2; i>=0; --i){
            rshift_v128_t(&t1,&VPs_v128_t[j_txt],i);
            and_v128_t(&t1,&t1,&one);
            rshift_v128_t(&t2,&VNs_v128_t[j_txt],i);
            and_v128_t(&t2,&t2,&one);
            if(bool_v128_t(&t1)){
              --scr;
              if(scr < m_vsc){
                m_vidx=i;
                m_vsc = scr;
                t_ih2 = i;
              }
            }
            else if(bool_v128_t(&t2)){
              ++scr;
            }
            else{
              // +-0
            }
          }
          if(m_hsc <= m_vsc){
            fprintf(stderr,"l_pat %d\n",l_pat);
            i_stt = l_pat-1;
            j_stt = m_hidx;
            fprintf(stderr,"j_stt m_hidx %d\n",j_stt);
            i_mxs = m_hidx;
            i_hmin = t_ih1;
          }
          else{
            fprintf(stderr,"m_vidx %d\n",m_vidx);
            i_stt = i_top+m_vidx;
            j_stt = j_txt;
            fprintf(stderr,"j_stt j_txt %d\n",j_stt);
            i_mxs = j_txt;
            i_hmin = t_ih2;
          }
          fprintf(stderr,"i_stt %d\n",i_stt);
        }

        int i_pat = i_stt;
        int i_txt = j_stt;
        fprintf(stderr,"i_txt %d\n",i_txt);
        fprintf(stderr,"i_mxs %d\n",i_mxs);// TODO bug test5.fa
        fprintf(stderr,"%d %d\n",i_pat,i_txt);
        int l_ali = 0;
        pat_buf[l_ali] = pat[i_pat];
        --i_pat;
        txt_buf[l_ali] = txt[i_txt];
        --i_txt;
        if(txt_buf[l_ali] == pat_buf[l_ali]){
          ali_buf[l_ali] = '|';
        }
        else{
          ali_buf[l_ali] = 'x';
        }
        ++l_ali;

        while(i_mxs>=0 && i_txt>=0 && i_pat>=0){
          //if(pat[i_pat] == txt[i_txt] || ((VPs[i_mxs] ^ HPs[i_mxs]) & ~(VNs[i_mxs] | HNs[i_mxs])) & (1ull << i_hmin))
          v128_t focus;
          v128_t vfocus;
          v128_t tmp;
          v128_t tmp2;
          v128_t tmp3;
          lshift_v128_t(&focus,&one,i_hmin);
          if(i_hmin == 0){
            xor_v128_t(&vfocus,&vfocus,&vfocus);
          }
          else{
            lshift_v128_t(&vfocus,&one,i_hmin-1);
          }
          and_v128_t(&tmp,&D0s_v128_t[i_mxs],&focus);
          and_v128_t(&tmp2,&HPs_v128_t[i_mxs],&focus);
          and_v128_t(&tmp3,&VPs_v128_t[i_mxs],&vfocus);
          if(pat[i_pat] == txt[i_txt] || !bool_v128_t(&tmp))
          //if(pat[i_pat] == txt[i_txt] || !(D0s[i_mxs] & (1ull << i_min)))
          //if(D0s[i_mxs] & (1ull << i_min) && pat[i_pat] == txt[i_txt])
          {
            --i_mxs;
            if(i_pat > l_pat-1){
              pat_buf[l_ali] = ' ';
              fprintf(stderr,"not come here 1 %d %d\n",i_pat,l_pat);
              exit(1);
            }
            else{
              pat_buf[l_ali] = pat[i_pat];
            }
            --i_pat;
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            if(txt_buf[l_ali] == pat_buf[l_ali]){
              ali_buf[l_ali] = '|';
            }
            else if(pat_buf[l_ali] == ' '){
              ali_buf[l_ali] = ' ';
            }
            else{
              ali_buf[l_ali] = 'x';
            }
          }
          else if(bool_v128_t(&tmp2))
          //else if(HPs[i_mxs] & (1ull << i_hmin))
          {
            --i_mxs;
            ++i_hmin;
            pat_buf[l_ali] = '-';
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            ali_buf[l_ali] = '*';
          }
          else if(bool_v128_t(&tmp3))
          //else if(VPs[i_mxs] & (1ull << i_hmin))
          {
            --i_hmin;
            if(i_pat > l_pat-1){
              pat_buf[l_ali] = ' ';
              fprintf(stderr,"not come here 2\n");
              exit(1);
            }
            else{
              pat_buf[l_ali] = pat[i_pat];
            }
            --i_pat;;
            txt_buf[l_ali] = '-';
            ali_buf[l_ali] = '@';
          }
          else{
            --i_mxs;
            if(i_pat > l_pat-1){
              pat_buf[l_ali] = ' ';
              fprintf(stderr,"not come here 3\n");
              exit(1);
            }
            else{
              pat_buf[l_ali] = pat[i_pat];
            }
            --i_pat;
            txt_buf[l_ali] = txt[i_txt];
            --i_txt;
            if(txt_buf[l_ali] == pat_buf[l_ali]){
              ali_buf[l_ali] = '|';
            }
            else if(pat_buf[l_ali] == ' '){
              ali_buf[l_ali] = ' ';
            }
            else{
              ali_buf[l_ali] = 'x';
            }
          }
          ++l_ali;
        }
        txt_buf[l_ali] = '\0';
        ali_buf[l_ali] = '\0';
        pat_buf[l_ali] = '\0';
        printf("i_mxs %d\n",i_mxs);
        printf("i_txt %d\n",i_txt);
        printf("i_pat %d\n",i_pat);
        printf("%s#\n",txt_buf);
        printf("%s#\n",ali_buf);
        printf("%s#\n",pat_buf);
        printf("##\n");
        {
          int i;
          int loop = l_ali/2;
          for(i=0; i<loop; ++i){
            char tmp;
            tmp = txt_buf[l_ali-1-i];
            txt_buf[l_ali-1-i] = txt_buf[i];
            txt_buf[i] = tmp;
            tmp = ali_buf[l_ali-1-i];
            ali_buf[l_ali-1-i] = ali_buf[i];
            ali_buf[i] = tmp;
            tmp = pat_buf[l_ali-1-i];
            pat_buf[l_ali-1-i] = pat_buf[i];
            pat_buf[i] = tmp;
          }
        }
        //printf("%s\n",txt_buf);
        //printf("%s\n",ali_buf);
        //printf("%s\n",pat_buf);
        //printf("##\n");

      }
      */
    }
    //printf("@%s\n%s\n+\t%s\n%s\n",&nls[i][1],reads[i],eds[i],eds[i]);
    return;
  }
}


void ola(char* txt, char* pat, char * txt_buf, char * pat_buf, char * ali_buf){
  int N;
  {
    int t_len = strlen(txt);
    int p_len = strlen(pat);
    if(t_len > p_len){
      N = 2*t_len+1;
    }
    else{
      N = 2*p_len+1;
      char * tmp;
      tmp = txt;
      txt = pat;
      pat = tmp;
    }
  }
  int ** m_distance = (int**)malloc(sizeof(int*)*N);
  int ** m_back_track = (int**)malloc(sizeof(int*)*N);
  {
    int i;
    for(i=0; i<N; ++i){
      m_distance[i] = (int*)malloc(sizeof(int)*N);
      if(m_distance[i] == NULL){
        fprintf(stderr, "malloc m_distance[%d]\n",i);
        exit(1);
      }
      m_back_track[i] = (int*)malloc(sizeof(int)*N);
      if(m_back_track[i] == NULL){
        fprintf(stderr, "malloc m_back_track[%d]\n",i);
        exit(1);
      }
    }
  }
  int * bt_str = (int*)malloc(sizeof(int)*N);
  if(bt_str == NULL){
    fprintf(stderr,"malloc bt_str\n");
    exit(1);
  }

  m_distance[0][0]=0;
  m_back_track[0][0]=-1;
  // TODO
  //  *bt=1;
  //  return diag;
  //  *bt=2;
  //  return vert;
  //  *bt=4;
  //  return hori;
  {
    int i;
    int loop=strlen(pat);
    for(i=1; i<=loop; ++i){
      m_distance[i][0]=i;
      m_back_track[i][0]=2;// vert
    }
    /*
    for(i=1; i<=AW; ++i){
      m_distance[i][0]=0;
      m_back_track[i][0]=-1;
    }
    for(i=AW+1; i<=loop; ++i){
      m_distance[i][0]=i-AW;
      m_back_track[i][0]=-1;
    }
    */
  }
  {
    int i;
    int loop=strlen(txt);
    for(i=1; i<=loop; ++i){
      m_distance[0][i]=i;
      m_back_track[0][i]=4;// hori
    }
    /*
    for(i=1; i<=AW; ++i){
      m_distance[0][i]=0;
      m_back_track[0][i]=-1;
    }
    for(i=AW+1; i<=loop; ++i){
      m_distance[0][i]=i-AW;
      m_back_track[0][i]=-1;
    }
    */
  }
  {
    int i,j;
    int l_pat=strlen(pat);
    int l_txt=strlen(txt);
    int MATCH=0;
    int MISMA=1;
    int GAPPE=1;
    for(i=1; i<=l_pat; ++i){
      for(j=1; j<=l_txt; ++j){
        int tmpScore = (pat[i-1]==txt[j-1]) ? MATCH : MISMA;
        m_distance[i][j] = min_of_the_3(tmpScore+m_distance[i-1][j-1], GAPPE+m_distance[i-1][j],GAPPE+m_distance[i][j-1], &m_back_track[i][j]);
        if(m_back_track[i][j] == 1 && m_distance[i][j] > m_distance[i-1][j-1]) m_back_track[i][j]=0; // mismatch
      }
    }
  }

  // back track
  int pat_is_suffix=0; // pat = i
  {
    int i=strlen(pat); int j=strlen(txt);
    // find i and j
    //printf("i,j=%d,%d\n",i,j);
    {
      int mi,mj,tmp; mi=i;mj=j;
      int min_ed=m_distance[i][j];
      int txt_is_prefix=0;
      for(tmp=1; tmp<i; ++tmp){
        if(m_distance[tmp][j] < min_ed){
          mi=tmp;
          min_ed=m_distance[tmp][j];
        }
      }
      for(tmp=1; tmp<j; ++tmp){
        if(m_distance[i][tmp] < min_ed){
          mj=tmp;
          min_ed=m_distance[i][tmp];
          txt_is_prefix=1;
          //pat_is_suffix=1;
        }
      }
      if(txt_is_prefix==1){
        j=mj;
        if(j-i+1<0) pat_is_suffix=1;
      }
      else{
        i=mi;
        if(!(i-j+1<0)) pat_is_suffix=1;
      }
    }
    //printf("i,j=%d,%d\n",i,j);
  
    int index=0;
    while(m_back_track[i][j] != -1){
      bt_str[index++]=m_back_track[i][j];
      if(m_back_track[i][j] == 1 || m_back_track[i][j] == 0){
        --i;
        --j;
      }
      else if(m_back_track[i][j] == 2){
        --i;
      }
      else{
        --j;
      }
    }
    {
      int loop=index/2;
      int i;
      for(i=0; i<loop; ++i){
        swap(&bt_str[i],&bt_str[index-1-i]);
      }
    }

    int i_pat=0; int i_txt=0;
    int pos=0;
    int idx;
    if(pat_is_suffix==0){
      for(idx=0; idx<j; ++idx,++pos){
        txt_buf[pos]=txt[i_txt++];
        pat_buf[pos]=' ';
        ali_buf[pos]=' ';
      }
    }
    else{
      for(idx=0; idx<i; ++idx,++pos){
        pat_buf[pos]=pat[i_pat++];
        txt_buf[pos]=' ';
        ali_buf[pos]=' ';
      }
    }
    for(idx=0; idx<index; ++idx,++pos){
      if(bt_str[idx] == 0){ // mismatch
        pat_buf[pos] = pat[i_pat++];
        txt_buf[pos] = txt[i_txt++];
        ali_buf[pos] = 'x';
      }
      else if(bt_str[idx] == 1){ // match
        pat_buf[pos] = pat[i_pat++];
        txt_buf[pos] = txt[i_txt++];
        ali_buf[pos] = '|';
      }
      else if(bt_str[idx] == 2){ // gap
        pat_buf[pos] = pat[i_pat++];
        txt_buf[pos] = '-';
        ali_buf[pos] = '@';
      }
      else if(bt_str[idx] == 4){ // gap
        pat_buf[pos] = '-';
        txt_buf[pos] = txt[i_txt++];
        ali_buf[pos] = '*';
      }
      else{
        fprintf(stderr,"sth strange when back tracing.\n");
        exit(1);
      }
    }
    pat_buf[pos] = '\0';
    txt_buf[pos] = '\0';
    ali_buf[pos] = '\0';
    {
      int tmp;
      tmp=strlen(pat);
      for(;i_pat<tmp; ++pos,++i_pat){
        pat_buf[pos]=pat[i_pat];
      }
      pat_buf[pos]='\0';
      tmp=strlen(txt);
      for(;i_txt<tmp; ++pos,++i_txt){
        txt_buf[pos]=txt[i_txt];
      }
      txt_buf[pos]='\0';
    }
  }

  {
    int i;
    for(i=0; i<N; ++i){
      free(m_distance[i]);
      free(m_back_track[i]);
    }
    free(m_distance);
    free(m_back_track);
  }
  free(bt_str);
}

void malloc_h64(unsigned long long ** D0s,unsigned long long ** HPs,unsigned long long ** VPs,unsigned long long ** HNs,unsigned long long ** VNs){
  *D0s = (unsigned long long*)malloc(sizeof(unsigned long long)*(max_read_length*2+1));
  if(D0s == NULL){
    fprintf(stderr,"cannot allocate memory: %s\n","D0s");
    exit(1);
  }
  *HPs = (unsigned long long*)malloc(sizeof(unsigned long long)*(max_read_length*2+1));
  if(HPs == NULL){
    fprintf(stderr,"cannot allocate memory: %s\n","HPs");
    exit(1);
  }
  *VPs = (unsigned long long*)malloc(sizeof(unsigned long long)*(max_read_length*2+1));
  if(VPs == NULL){
    fprintf(stderr,"cannot allocate memory: %s\n","VPs");
    exit(1);
  }
  *HNs = (unsigned long long*)malloc(sizeof(unsigned long long)*(max_read_length*2+1));
  if(HNs == NULL){
    fprintf(stderr,"cannot allocate memory: %s\n","HNs");
    exit(1);
  }
  *VNs = (unsigned long long*)malloc(sizeof(unsigned long long)*(max_read_length*2+1));
  if(VNs == NULL){
    fprintf(stderr,"cannot allocate memory: %s\n","VNs");
    exit(1);
  }
}

void free_h64(unsigned long long ** D0s,unsigned long long ** HPs,unsigned long long ** VPs,unsigned long long ** HNs,unsigned long long ** VNs){
  free(*D0s);
  free(*HPs);
  free(*VPs);
  free(*HNs);
  free(*VNs);
}

