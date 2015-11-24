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
void kt_for(int n_threads, void (*func)(void*,long,int), void *data, long n);
}

using namespace std;

struct depthInterval {
   uint32_t start,stop,depth,gq; 
   int chrom; 
};

class depthReader {
 public:
  depthReader(const char *depth_fname);
  ~depthReader();
  int next();
  int nread,nsample;
  depthInterval line;
  bool open;
 private:
  gzFile fp;
  int buf[5];
};

class circularBuffer {
 public:

  circularBuffer() {
    resize(0);
  }

  circularBuffer(int buffer_size) {
    resize(buffer_size);
  }

  void resize(int buffer_size) {
    _bufsize=buffer_size;
    _offset=_bufn=0;
    _buf = new depthInterval[_bufsize];
    //_buf.resize(_bufsize);
  }

  ~circularBuffer() {    delete[] _buf;  }

  int size() {    return(_bufn);  }
  bool  empty() {    return(_bufn==0);  }
  bool  full() {    return(_bufn==_bufsize);  }

  bool push_back(depthInterval o) {
    int idx=(_offset + _bufn)%_bufsize;
    //    cerr << "push_back " << _offset << " " <<_bufn<<" "<<idx<<" "<<endl;
    assert(idx<_bufsize&&idx>=0);
    if(_bufn<_bufsize) _buf[idx]=o;
    else die("buffer overflow");
    _bufn++;
    return(true);
  }

  depthInterval *front() { 
    assert(_bufn>0);
    assert(_offset>=0 && _offset<_bufsize);      
    return(&_buf[_offset]);
  }

  depthInterval *back() { 
    assert(_bufn>0);
    return(&_buf[(_offset+_bufn-1)%_bufsize]);
  }

  void pop_front() {
    if(_bufn>0) {
      _offset++;
      _offset %= _bufsize;
      _bufn--;
    }
    else 
      die("buffer empty (pop_front)");
  }
 private:
  //  vector<depthInterval> _buf;
  depthInterval *_buf;
  int _bufsize,_offset,_bufn;  
};


class depthMerger {
 public:
  int _nsample;
  int32_t *dp,*gq;
  depthMerger(vector<string> & files);
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


  pthread_t *_threads;
  const static  int buf_size=1000;

  int cur_pos;
  int curr_chrom;
  //  vector< deque< depthInterval > > dp_buf;//stores the depth intervals for nsamples.
  circularBuffer *dp_buf;//stores the depth intervals for nsamples.  
  vector< depthReader* > r;
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

  struct  next_args *_dp_buf_args;
};

