#pragma once
#include "utils.h"

#include <math.h>
#include <map>
#include <vector>


extern "C" {
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
}

struct tabLine {
  int a,b;
  string chrom,id;
  int dp;
};

class tabReader {
 public:
  tabReader(const string& bed_fname,const  string& region);
  ~tabReader();
  int next();
  int nread;
  tabLine line;
  bool open;
 private:
  tbx_t *tbx;
  BGZF *fp;
  hts_itr_t *itr;  
  kstring_t s;
  string sampleid;
};
