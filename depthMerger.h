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
#include "string.h"

#include "utils.h"

extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"
}

using namespace std;

class depthInterval {
 public:
  depthInterval() { };
  depthInterval(uint32_t a,uint32_t b,int c,int d) {
    start=a;
    stop=b;
    depth=c;
    gq=d;
  };
  uint32_t start,stop,depth,gq;
};

class depthReader {
 public:
  depthReader(const char *depth_fname);
  ~depthReader();
  int next();
  int nread,nsample;
  depthInterval line;
  int chrom;
  bool open;
 private:
  gzFile fp;
  int buf[5];
};



class depthMerger {
 public:
  int _nsample;
  int32_t *dp,*gq;
  depthMerger(vector<string> & files);
  ~depthMerger();
  int writeDepthMatrix(const char *output_file);
  int next();
  int findCurrPos();
  bcf_hdr_t *makeDepthHeader();
  int getCurChr();
  int getCurPos();
  bcf_hdr_t *getHeader();
  void fillBuffer(int i);
  int setThreads(int nthreads);
  bool anyOpen();
  int _nthreads;

 private:
  pthread_t *_threads;
  const static  int buf_size=10000;

  int cur_pos;
  int curr_chrom;
  vector< deque< depthInterval > > dp_buf;//stores the depth intervals for nsamples.
  vector< depthReader* > r;
  bcf_srs_t *sr;
  vector<string> _files;
  bcf_hdr_t *_hdr;
  bool _eof_warn;
  pthread_mutex_t *_dp_buf_mutex;
  pthread_cond_t _less,_more;//signals buffer was incremented/decrermented.
  void  unlockDepthBuffer();
  void lockDepthBuffer();
  bool checkBufferIsOkayToRead();
  void startReadBuffer();
  
  struct  next_args *_dp_buf_args;
};
