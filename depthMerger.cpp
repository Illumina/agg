#include "depthMerger.h"
#include <pthread.h> 
#define THREADED 0
void depthReader::open(const char *depth_fname) { 
  //  cerr<<"Opening "<<depth_fname<<endl;
  _fp = gzopen(depth_fname, "r");
  if(_fp==NULL) {
    cerr<<"problem opening "<< depth_fname<<endl;
    exit(1);
  }
  assert(  gzread(_fp,_buf,20) );
  //  for(int i=0;i<5;i++) cerr << _buf[i] << " "; cerr << endl;
  _open=true;
}

bool depthReader::next() {
  if(gzread(_fp,_buf,20)) {
    return(1);
  }
  else {
    _open=false;
    return(0);
  }
}

int depthMerger::getCurChr() { return(_cur_chr); };
int depthMerger::getCurPos() { return(_cur_pos); };
bcf_hdr_t *depthMerger::getHeader() {return(_hdr);};

bcf_hdr_t *depthMerger::makeDepthHeader() {
  _hdr = bcf_hdr_init("w");
  bcf_hdr_append(_hdr, "##source=agg");
  bcf_hdr_append(_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">");
  bcf_hdr_append(_hdr,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  char clash_prefix[10]; //prefix for duplicated samples (borrowed from vcfmerge.c)
  for(int i=0;i<_nsample;i++) {
    snprintf(clash_prefix,10,"%d",i+1);
    bcf_hdr_t *src_hdr = sr->readers[i].header;
    if(i==0) copyContigs(src_hdr,_hdr);
    if(bcf_hdr_nsamples(src_hdr)!=1) {
      cerr<<"ERROR: file "<<sr->readers[i].fname<<" had >1 sample.  This cannot be output from agg ingest1!"<<endl;
      exit(1);
    }
    char *name=src_hdr->samples[0];
    if ( bcf_hdr_id2int(_hdr, BCF_DT_SAMPLE, name)!=-1 ) {
      if(!_force_samples)
	die("duplicate sample names.\n Try --allow-duplicates if you want to merge anyway");
      
      int len = strlen(src_hdr->samples[0]) + strlen(clash_prefix) + 1;
      name = (char*) malloc(sizeof(char)*(len+1));
      sprintf(name,"%s:%s",clash_prefix,src_hdr->samples[0]);
      bcf_hdr_add_sample(_hdr,name);
      free(name);     
    }
    else
      bcf_hdr_add_sample(_hdr,name);      
  }
  bcf_hdr_sync(_hdr);
  cerr<<  bcf_hdr_nsamples(_hdr) << " samples in merged chunk"<<endl;
  return(_hdr);
}

depthMerger::~depthMerger() {
  bcf_sr_destroy(sr);   
  delete[]  _dp_buf_mutex;
  delete[] _dp_rdr;
}

int depthMerger::setThreads(int nthreads) {
  _nthreads = min(_nsample,nthreads);
  _threads = new pthread_t[_nthreads];
  return(_nthreads);
}

void depthMerger::setForceSamples(int f) {
  _force_samples=f;
}

depthMerger::depthMerger(vector<string> & files,bool force_samples) {
  setForceSamples(force_samples);
  _nthreads=1;
  _files = files;
  _eof_warn=false;
  _nsample=files.size();
  _dp_buf_mutex = new pthread_mutex_t[_nsample];
  sr=bcf_sr_init() ;
  sr->require_index=1;
  for(int i=0;i<_nsample;i++) {
    _dp_buf_mutex[i] = PTHREAD_MUTEX_INITIALIZER;
    if(!bcf_sr_add_reader (sr, files[i].c_str())) {
      fprintf(stderr,"problem opening bcf");
      exit(1);
    }
  }
  cerr << "nsample="<<_nsample<<endl;
  _dp_rdr = new depthReader[_nsample];
  for(int i=0;i<_nsample;i++)  {
    string tmp = files[i];
    tmp = tmp.substr(0,tmp.size()-3) + "tmp";
    _dp_rdr[i].open(tmp.c_str());
  }
  cerr << "opened "<<_nsample<<" .tmp files for merging..."<<endl;    

  _cur_chr= -1;
  _cur_pos = -1;
  makeDepthHeader();
}

int depthMerger::startNewChromosome() {
  cerr << "starting new chromosome..."<<endl;
  _cur_pos = -1;
  int prev_chrom=_cur_chr;
  _cur_chr=-1;
  for(int i=0;i<_nsample;i++) {    
#if THREADED==1
    lock(i);
    //do something to update the reader
#endif
    if(_dp_rdr[i].isOpen()) {
      if(_cur_chr==-1 || _dp_rdr[i].getChrom() < _cur_chr) 
	_cur_chr=_dp_rdr[i].getChrom();
    }
#if THREADED==1
    unlock(i);
#endif
  }

  for(int i=0;i<_nsample;i++) {    
#if THREADED==1
    lock(i);
    //do something to update the reader
#endif
    if(_dp_rdr[i].isOpen()) {
      if(_cur_pos==-1 || _dp_rdr[i].getChrom() < _cur_pos) 
	_cur_pos=_dp_rdr[i].getPos();
    }
#if THREADED==1
    unlock(i);
#endif
  }

  assert(_cur_chr!=prev_chrom);
  cerr << "bcf_hdr_nsamples(_hdr) = "<<bcf_hdr_nsamples(_hdr)<<endl;
  cerr << "chromosome start = " << bcf_hdr_id2name(_hdr,_cur_chr)<<":"<< _cur_pos<<endl;
  return(_cur_pos);
}

void depthMerger::lockDepthBuffer() {
  for(int i=0;i<_nsample;i++)  pthread_mutex_lock(&_dp_buf_mutex[i]);
}

void depthMerger::unlockDepthBuffer() {
  for(int i=0;i<_nsample;i++)  pthread_mutex_unlock(&_dp_buf_mutex[i]);
}

typedef struct next_args {
  depthMerger *dm;
  int offset;
} next_args;


void depthMerger::lock(int i) {
  pthread_mutex_lock(&_dp_buf_mutex[i]);  
}
void depthMerger::unlock(int i) {
  pthread_mutex_unlock(&_dp_buf_mutex[i]);  
}
void depthMerger::wait_less(int i) {
  //  pthread_cond_wait(&_less[0],&_dp_buf_mutex[i]);
    pthread_cond_wait(&_less[i],&_dp_buf_mutex[i]);
}
void depthMerger::wait_more(int i) {
  pthread_cond_wait(&_more[i],&_dp_buf_mutex[i]);
}

void depthMerger::more(int i){
  pthread_cond_signal(&_more[i]);
}

void depthMerger::less(int i){
  pthread_cond_signal(&_less[i]);
  //  pthread_cond_broadcast(&_less[0]);
}

void *buffer_updater(void *_g) {
  next_args *a = (next_args *)_g;
  depthMerger *d = a->dm;
  int nopen=1;
  cerr << "buffer_updater "<<a->offset<<endl;
  while(nopen>0) {
    nopen=0;
    for(int i=a->offset;i<d->_nsample;i+=d->_nthreads) {
      d->lock(i);
      if(d->_dp_rdr[i].isOpen()) nopen++;
      //do stuff to update dp_rdr
      // while(d->dp_buf[i].full()) {
      // 	//	cerr << "buffer_updater "<<a->offset<<" is full size="<<d->dp_buf[i].size()<<endl;
      // 	d->wait_less(i);
      // }

      // if(!d->dp_buf[i].full())  {
      // 	d->fillBuffer(i);
      // 	d->more(i);
      // }      
      d->unlock(i);
    }
  }
  cerr << "buffer_updater " << a->offset<<" finished"<<endl;
  pthread_exit((void*) 0);
}

void depthMerger::startReadBuffer() {
  _dp_buf_args = new next_args[_nthreads];
  assert(_nthreads>0);
  _less = new pthread_cond_t[_nsample];
  _more = new pthread_cond_t[_nsample];
  for(int i=0;i<_nsample;i++) {
    _less[i]=PTHREAD_COND_INITIALIZER;
    _more[i] = PTHREAD_COND_INITIALIZER;
  }
  for(int i=0;i<_nthreads;i++) {
    _dp_buf_args[i].dm=this;
    _dp_buf_args[i].offset=i;
    pthread_create(&_threads[i], NULL, buffer_updater, (void *)&_dp_buf_args[i]);
  }
};

//kt_for funciton - this is too slow
// static void update_buffer(void *_data, long i, int tid) {
//   depthMerger *d = (depthMerger *)_data;
//   d->dp[i]=bcf_int32_missing; //default.
//   d->gq[i]=bcf_int32_missing; //default.
//   int _cur_pos=_cur_pos;
//   //remove intervals behnd _cur_pos
//   while(!d->dp_buf[i].empty() && (d->dp_buf[i].front()->stop < _cur_pos || d->dp_buf[i].front()->chrom<d->_cur_chr) )      
//     d->dp_buf[i].pop_front();  

//   if(d->dp_buf[i].empty())
//     d->fillBuffer(i);

//   if(!d->dp_buf[i].empty() && d->dp_buf[i].front()->start <= _cur_pos && d->dp_buf[i].front()->stop >= _cur_pos && d->dp_buf[i].front()->chrom == d->_cur_chr){
//     d->dp[i] = d->dp_buf[i].front()->depth;
//     d->gq[i] = d->dp_buf[i].front()->gq;
//   }
// }

int depthMerger::next() {

  if(_cur_pos==-1) {
    startNewChromosome();
  }
  else _cur_pos++;
  assert(_cur_pos != -1);

  //  cerr << "next " << bcf_hdr_id2name(_hdr,_cur_chr) <<":"<< _cur_pos<<endl;//
  int n_on_chromosome=0;
  int nopen=0;

  //kt_for version - slow didnt bother finishing running it
  // for(int i=0;i<_nsample;i++)
  //   if( r[i]->open ) 
  //     nopen++;
  // kt_for(_nthreads,update_buffer,(void *)this,_nsample);

  // for(int i=0;i<_nsample;i++)
  //   if(!dp_buf[i].empty() && dp_buf[i].front()->chrom == _cur_chr)
  //     n_on_chromosome++;    

  //single threaded version - reasonably fast. 90 seconds with 8 threaded bgzip compression.
#if THREADED==0
  for(int i=0;i<_nsample;i++)  {              

    //remove intervals behnd _cur_pos
    while( _dp_rdr[i].isBehind(_cur_chr,_cur_pos) )
      _dp_rdr[i].next();

    if(_dp_rdr[i].contains(_cur_chr,_cur_pos) ) {
      dp[i] = _dp_rdr[i].getDepth();
      gq[i] = _dp_rdr[i].getGQ();
    } 
    else  {//set missing
      dp[i]=bcf_int32_missing; 
      gq[i]=bcf_int32_missing;
    }
    //cerr << _dp_rdr[i].getDepth() << " ";
    if(_dp_rdr[i].isOpen()) {
      if(_dp_rdr[i].getChrom()==_cur_chr)
	n_on_chromosome++;
      nopen++;
    }
  }
#else
  //producer-consumer implementation. unmitigated disaster. 320 seconds using 16 threads. nthreads must equal nsample or you get a race condition.
  for(int i=0;i<_nsample;i++)  {              
    lock(i);
    if( r[i]->open )       nopen++;

    dp[i]=bcf_int32_missing; //default.
    gq[i]=bcf_int32_missing; //default.
    //remove intervals behnd _cur_pos
    while(!dp_buf[i].empty() && (dp_buf[i].front()->stop < _cur_pos || dp_buf[i].front()->chrom<_cur_chr) )      
      dp_buf[i].pop_front();  
    less(i);
    
    // i is behind _cur_pos. needs data.
    while(dp_buf[i].empty())
      wait_more(i);

    if(!dp_buf[i].empty() && dp_buf[i].front()->start <= _cur_pos && dp_buf[i].front()->stop >= _cur_pos && dp_buf[i].front()->chrom == _cur_chr){
      dp[i] = dp_buf[i].front()->depth;
      gq[i] = dp_buf[i].front()->gq;
    }
    if(!dp_buf[i].empty() && dp_buf[i].front()->chrom == _cur_chr)
      n_on_chromosome++;    
    unlock(i);
  }
#endif
  //all done.
  if(nopen==0)    
    return(0);

  //if chromosome is finished, move to next one
  if(n_on_chromosome==0) {
    cerr << "chromosome end   = " << bcf_hdr_id2name(_hdr,_cur_chr) <<":"<< _cur_pos<<endl;//
    cerr << "nopen="<<nopen<<endl;
    _cur_pos = -1;
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
    bcf_clear1(line) ;
    line->rid = dm->getCurChr();
    line->pos = dm->getCurPos();
    bcf_update_alleles_str(hdr, line, "N,.");
    // cerr<<line->pos+1<<"\t";
    // for(int i=0;i<dm->_nsample;i++)      cerr<<dm->dp[i]<< " ";    cerr<<endl;
    bcf_update_format_int32(hdr, line,"DP",dm->dp,dm->_nsample);
    bcf_update_format_int32(hdr, line,"GQ",dm->gq,dm->_nsample);      
    if(line->pos%1000000==0)        cerr<<bcf_hdr_id2name(hdr,line->rid)<<":"<<line->pos+1<<endl;    
    pthread_cond_signal(&(args->more));
    pthread_mutex_unlock (&(teh_mutex));
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


int depthMerger::writeDepthMatrix(const char *output_file) {

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
  cerr << "first pos was "<<_cur_pos<<endl;
  ta.dm=this;
  ta.out_fh =   hts_open(output_file,"wb9");
  //  htsFile *out_fh = ta.out_fh;
  if(_nthreads>0)  
    hts_set_threads(ta.out_fh,max(_nthreads,1) );
  bcf_hdr_write(ta.out_fh, _hdr);
  ta.less = PTHREAD_COND_INITIALIZER;
  ta.more = PTHREAD_COND_INITIALIZER;

#if THREADED==1
  startReadBuffer(); 
#endif

  bcf1_t* line =  bcf_init();
  while(next()) {
    line->rid = getCurChr();
    line->pos = getCurPos();
    if(line->pos%1000000==0)        
      cerr<<bcf_hdr_id2name(_hdr,line->rid)<<":"<<line->pos+1<<endl;    

    bcf_update_alleles_str(_hdr, line, "N,.");
    bcf_update_format_int32(_hdr, line,"DP",dp,_nsample);
    bcf_update_format_int32(_hdr, line,"GQ",gq,_nsample);      
    bcf_write1(ta.out_fh, _hdr, line);
    bcf_clear1(line);
  }

  // cerr << "debug2"<<endl;
  // pthread_t consumer;
  // pthread_t producer;
  // pthread_create(&producer, NULL, producer_func, (void *)&ta);
  // pthread_create(&consumer, NULL, consumer_func, (void *)&ta);
  // cerr << "debug3"<<endl;
  // pthread_join(producer, NULL);
  // pthread_join(consumer, NULL);


  cerr << "finished creating dpt file"<<endl;
  hts_close(ta.out_fh);
  bcf_hdr_destroy(_hdr);
  delete[] dp;
  delete[] gq;

  return(0);
}


