#include "depthMerger.h"
#include <pthread.h> 

depthReader::depthReader(const char *depth_fname) {
  //  cerr<<"Opening "<<depth_fname<<endl;
  fp = gzopen(depth_fname, "r");
  if(fp==NULL) {
    cerr<<"problem opening "<< depth_fname<<endl;
    exit(1);
  }
  open=true;
}

depthReader::~depthReader() {
  gzclose(fp);
}

int depthReader::next() {

  if(gzread(fp,buf,20)) {
    chrom=buf[0];
    line.start=buf[1];
    line.stop=buf[2];
    line.depth=buf[3];
    line.gq=buf[4];
    return(1);
  }
  else {
    open=false;
    return(0);
  }
}

int depthMerger::getCurChr() { return(curr_chrom); };
int depthMerger::getCurPos() { return(cur_pos); };
bcf_hdr_t *depthMerger::getHeader() {return(_hdr);};

bcf_hdr_t *depthMerger::makeDepthHeader() {
  _hdr = bcf_hdr_init("w");
  bcf_hdr_append(_hdr, "##source=agg");
  bcf_hdr_append(_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">");
  bcf_hdr_append(_hdr,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");

  for(int i=0;i<_nfile;i++) {
    bcf_hdr_t *src_hdr = sr->readers[i].header;
    if(i==0) copyContigs(src_hdr,_hdr);
    if(bcf_hdr_nsamples(src_hdr)!=1) {
      cerr<<"ERROR: file "<<sr->readers[i].fname<<" had >1 sample.  This cannot be output from agg ingest1!"<<endl;
      exit(1);
    }
    char *sample_name=src_hdr->samples[0];
    if ( bcf_hdr_id2int(_hdr, BCF_DT_SAMPLE, sample_name)!=-1 )
      die("duplicate sample names. use --force-samples if you want to merge anyway");
    bcf_hdr_add_sample(_hdr,sample_name);      
  }
  bcf_hdr_sync(_hdr);
  cerr<<  bcf_hdr_nsamples(_hdr) << " samples in merged chunk"<<endl;
  return(_hdr);
}

depthMerger::~depthMerger() {
  bcf_sr_destroy(sr);   
}

depthMerger::depthMerger(vector<string> & files) {
  _files = files;
  _eof_warn=false;
  _nsample=_nfile=files.size();

  sr=bcf_sr_init() ;
  sr->require_index=1;
  for(int i=0;i<_nfile;i++) {
    if(!bcf_sr_add_reader (sr, files[i].c_str())) {
      fprintf(stderr,"problem opening bcf");
      exit(1);
    }
  }

  cerr << "nsample="<<_nsample<<endl;
  for(int i=0;i<_nfile;i++)  {
    string tmp = files[i];
    tmp = tmp.substr(0,tmp.size()-3) + "tmp";
    r.push_back(new depthReader(tmp.c_str()));
  }
  cerr << "opened "<<_nsample<<" .tmp files for merging..."<<endl;    
  dp_buf.resize(_nsample);
  for(int i=0;i<_nsample;i++)  {
    assert(r[i]->next());
    dp_buf[i].push_back( depthInterval(r[i]->line) );
    //    cerr <<  i<<": "<<dp_buf[i].front().start << endl;
  }
  curr_chrom=r[0]->chrom;
  cur_pos = -1;
  makeDepthHeader();
}

int depthMerger::findCurrPos() {
  cur_pos = -1;
  for(int i=0;i<_nsample;i++) {    
    if(!dp_buf[i].empty() && dp_buf[i].front().start < cur_pos)
      cur_pos = dp_buf[i].front().start;
  }
  cerr << "bcf_hdr_nsamples(_hdr) = "<<bcf_hdr_nsamples(_hdr)<<endl;
  cerr << curr_chrom<< ":"<<cur_pos << " line103"<<  endl;
  cerr << "chromosome start = " << bcf_hdr_id2name(_hdr,curr_chrom)<<":"<< cur_pos<<endl;
  return(cur_pos);
}

extern "C" {
  void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
  //  void kt_for(int n_threads, int n_items, void (*func)(void*,int,int), void *data);
}

//work function for updating all the dp_buffers
static void buffer_update(void *_g, long i, int tid) {
  depthMerger *dm = (depthMerger *) _g;
  //  cerr << "buffer_update "<<i<<endl;
  dm->fillBuffer(i);
}

void depthMerger::fillBuffer(int i) {
  dp[i]=bcf_int32_missing; //default.
  gq[i]=bcf_int32_missing; //default.
  if(r[i]->chrom==curr_chrom) {//add to buffer until we pass cur_pos
    while(r[i]->chrom<=curr_chrom && (dp_buf[i].size() < buf_size || dp_buf[i].back().stop < cur_pos)) {
      if(r[i]->next() && r[i]->chrom==curr_chrom)
	dp_buf[i].push_back( depthInterval(r[i]->line));
      else
	break;
    }
  }
  //discard intervals behind pos
  while(!dp_buf[i].empty() && dp_buf[i].front().stop < cur_pos) 
    dp_buf[i].pop_front();	      
  
  //update depth array.
  if(!dp_buf[i].empty() && dp_buf[i].front().start <= cur_pos && dp_buf[i].front().stop >= cur_pos){
    dp[i] = dp_buf[i].front().depth;
    gq[i] = dp_buf[i].front().gq;
  }
  //  cerr << "buf_size["<<i<<"] = "<<dp_buf[i].size()<<endl;
}

int depthMerger::next() {
  if(cur_pos==-1) findCurrPos();
  else cur_pos++;

  bool buf_refresh = false;
  for(int i=0;i<_nsample;i++)  {
    fillBuffer(i);
    if(dp_buf[i].back().stop < cur_pos) buf_refresh=true;
  }
  if(buf_refresh)
    kt_for(4, buffer_update, (void *)this,_nsample);

  //check if we have reached end of chromosome
  bool all_empty=true;
  for(int i=0;i<_nsample;i++)
    if(!dp_buf[i].empty()||(r[i]->chrom==curr_chrom && r[i]->open))
      all_empty=false;

  //if chromosome is finished, move to next one
  if(all_empty) {

    cerr << "ALL EMPTY pos="<<cur_pos+1<<endl;
    for(int i=0;i<_nsample;i++)
      cerr<< r[i]->chrom << "/" << r[i]->open << " ";
    cerr <<endl;
    
    int nopen=0;
    for(int i=0;i<_nsample;i++)
      if(r[i]->open)
	nopen++;

    if(nopen==0) 
      return(0);
    cerr << "chromosome end   = " << bcf_hdr_id2name(_hdr,curr_chrom) <<":"<< cur_pos<<endl;//

    //sweeps for the new chrom (one with the smallest index)
    int prev_chrom=curr_chrom;
    curr_chrom=r[0]->chrom;
    for(int i=1;i<_nsample;i++) 
      if(r[i]->chrom < curr_chrom && r[i]->open) 
	curr_chrom=r[i]->chrom;
    assert(curr_chrom!=prev_chrom);
    for(int i=0;i<_nsample;i++)  {
      if(r[i]->open)  {
	if(r[i]->chrom == curr_chrom)
	  dp_buf[i].push_back( depthInterval(r[i]->line));
      }
      else {
	if(!_eof_warn) {
	  cerr << "WARNING: " << _files[i] << " is finished but "<<nopen<<" .tmp files claim to have more data!"<<endl;
	  _eof_warn=true;
	}
      }
    }
    cur_pos=-1;
    return(next());
  }
  return(1);
}

pthread_mutex_t teh_mutex = PTHREAD_MUTEX_INITIALIZER;

typedef struct{
  depthMerger *dm;
  bcf1_t **buf;
  int buf_size,buf_off,buf_n;
  htsFile *out_fh;
  int fin;
  //  pthread_mutex_t mutex;
  pthread_cond_t less,more;//signals buffer was incremented/decrermented.
} thread_args;

void *producer_func(void *ptr) {
  thread_args *args = (thread_args *)ptr;
  assert(args->buf_off==0);
  assert(args->buf_size>0);
  assert(args->buf_n<args->buf_size);
  depthMerger *dm =   args->dm;
  bcf_hdr_t *hdr = dm->getHeader();
  cerr<<"\n producer has "<<dm->_nsample<<" samples"<<endl;
  cerr<<"dm->getCurPos() = "<<dm->getCurPos()<<endl;;
  while(dm->next()) {//gets a new line
    //    cerr<<"dm->getCurPos() = "<<dm->getCurPos()<<endl;;
    pthread_mutex_lock (&(teh_mutex));
    while(args->buf_n>=args->buf_size) 
      pthread_cond_wait(&args->less, &(teh_mutex));
    assert(args->buf_n<args->buf_size);
    bcf1_t *line= (args->buf)[(args->buf_off + args->buf_n)%args->buf_size];
    args->buf_n++;//increment number of items
    //    pthread_mutex_unlock (&(teh_mutex));
    bcf_clear1(line) ;
    line->rid = dm->getCurChr();
    line->pos = dm->getCurPos();
    bcf_update_alleles_str(hdr, line, "N,.");
    bcf_update_format_int32(hdr, line,"DP",dm->dp,dm->_nsample);
    bcf_update_format_int32(hdr, line,"GQ",dm->gq,dm->_nsample);      

    if(line->pos%10000000==0)        cerr<<bcf_hdr_id2name(hdr,line->rid)<<":"<<line->pos+1<<endl;    
    //    pthread_mutex_lock (&(teh_mutex));
    pthread_cond_signal(&(args->more));
    pthread_mutex_unlock (&(teh_mutex));
    //    cerr<<"args->buf_n="<<args->buf_n<<endl;
  }
  pthread_mutex_lock (&(teh_mutex));
  args->fin = 1;
  args->buf_n=2;
  pthread_cond_signal(&(args->more));
  pthread_mutex_unlock (&(teh_mutex));
  cerr << "exiting consumer thread"<<endl;
  pthread_exit((void*) 0);
}

void *consumer_func(void *ptr) {
  thread_args *args = (thread_args *)ptr;
  bcf_hdr_t *hdr = args->dm->getHeader();
  assert(args->buf_off==0);
  assert(args->buf_size>0);
  assert(args->buf_n<args->buf_size);
  while(1) {
    pthread_mutex_lock (&(teh_mutex));
    while(args->buf_n<2)//we do <2 rather <1 because the producer might be writing to buf[1]
      pthread_cond_wait(&args->more, &(teh_mutex));
    if(  args->fin )   
      pthread_exit((void*) 0);
    assert(args->buf_n>=2); 
    bcf1_t *line = (args->buf)[args->buf_off]; 
    bcf_write1(args->out_fh, hdr, line) ;
    args->buf_off++;//move offset forward
    args->buf_off %= args->buf_size; //reset offset when it surpasses buffer lenght
    args->buf_n--;//decrement number of items
    pthread_cond_signal(&(args->less));
    pthread_mutex_unlock (&(teh_mutex));  
  }
}


int depthMerger::writeDepthMatrix(const char *output_file,int nthreads) {
  pthread_t threads[2];//1 thread to produce. 1 thread to consume 
  cerr << "Writing out "<<output_file<<endl;
  dp = new int32_t[_nsample];
  gq = new int32_t[_nsample];

  thread_args ta;
  ta.fin=0;
  ta.buf_size=10;
  ta.buf_off=0;
  ta.buf_n=0;
  ta.buf = new bcf1_t*[ta.buf_size];
  for(int i=0;i<ta.buf_size;i++) ta.buf[i]=bcf_init();
  cerr << "first pos was "<<cur_pos<<endl;
  ta.dm=this;
  ta.out_fh =   hts_open(output_file,"wb9");
  htsFile *out_fh = ta.out_fh;
  if(nthreads>0)  hts_set_threads(ta.out_fh,max(nthreads/2,1) );
  bcf_hdr_write(ta.out_fh, _hdr);
  ta.less = PTHREAD_COND_INITIALIZER;
  ta.more = PTHREAD_COND_INITIALIZER;

  pthread_create(&threads[0], NULL, producer_func, (void *)&ta);
  pthread_create(&threads[1], NULL, consumer_func, (void *)&ta);
  pthread_join(threads[0], NULL);
  pthread_join(threads[1], NULL);

  // bcf1_t *line=bcf_init();
  // while( next() ) {
  //   //producer
  //   bcf_clear1(line) ;
  //   line->rid = curr_chrom;
  //   line->pos = cur_pos;
  //   bcf_update_alleles_str(_hdr, line, "N,.");
  //   bcf_update_format_int32(_hdr, line,"DP",dp,_nsample);
  //   bcf_update_format_int32(_hdr, line,"GQ",gq,_nsample);
  //   //consumer
  //   bcf_write1(out_fh, _hdr, line) ;
  // }

  cerr << "finished."<<endl;
  hts_close(ta.out_fh);
  bcf_hdr_destroy(_hdr);
  delete[] dp;
  delete[] gq;
  return(0);
}

