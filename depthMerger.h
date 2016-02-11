#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <limits>
#include "agg.h"
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <map>
#include "string.h"

#include "utils.h"

extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
}

using namespace std;

struct depthInterval {
   uint32_t start,stop,depth,gq; 
   int chrom; 
};

class depthReader {
 public:

  depthReader() { _open=false; }
  depthReader(const char *depth_fname) {       open(depth_fname) ;  }
  ~depthReader() {  gzclose(_fp);}
  void open(const char *depth_fname);
  bool next();
  int getDepth() {return(_buf[3]);}
  int getGQ() {return(_buf[4]);}
  int getChrom() {return(_buf[0]);}
  int getPos() {return(_buf[1]);}
  bool isAhead(int rid,int pos)  { if(_open) return(_buf[0]>rid || (_buf[0]==rid && _buf[1]>pos)) ; else return(false);}
  bool isBehind(int rid,int pos) { if(_open) return(_buf[0]<rid || (_buf[0]==rid && _buf[2]<pos)) ; else return(false);}
  bool contains(int rid,int pos) { if(_open) return(_buf[0]==rid && _buf[1]<= pos && _buf[2]>=pos); else return(false);}
  bool isOpen() {return _open;}
  int nread,nsample;
 private:
  gzFile _fp;
  int _buf[5];
  depthInterval _line;
  bool _open;
};

class depthMerger {
 public:
  int _nsample;
  int32_t *dp,*gq;
  depthMerger(vector<string> & files,bool force_samples=0);
  depthMerger();
  ~depthMerger();
  int writeDepthMatrix(const char *output_file);
  int next();
  int startNewChromosome();
  bcf_hdr_t *makeDepthHeader();
  int getCurChr();
  int getCurPos();
  bcf_hdr_t *getHeader();
  void fillBuffer(int i);
  int setThreads(int nthreads);
  bool anyOpen();
  int _nthreads;
  void lock(int i);
  void unlock(int i);
  void less(int i);
  void wait_less(int i);
  void more(int i);
  void wait_more(int i);
  void setForceSamples(int f);

  pthread_t *_threads;
  const static  int buf_size=1000;

  int _cur_pos, _cur_chr;
  //  vector< deque< depthInterval > > dp_buf;//stores the depth intervals for nsamples.
  depthReader* _dp_rdr;
  bcf_srs_t *sr;
  vector<string> _files;
  bcf_hdr_t *_hdr;
  bool _eof_warn;
  pthread_mutex_t *_dp_buf_mutex;
  pthread_cond_t *_less,*_more;//signals buffer was incremented/decrermented.
  void  unlockDepthBuffer();
  void lockDepthBuffer();
  bool checkBufferIsOkayToRead();
  void startReadBuffer();
  bool _force_samples;
  struct  next_args *_dp_buf_args;
};

