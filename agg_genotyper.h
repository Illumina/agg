#pragma once
#include "agg.h"
#include "htslib/hts.h"

#include "utils.h"
//#include "tabReader.h"
#include <deque>

class depthInterval {
 public:
  depthInterval(uint32_t a,uint32_t b,int c,int d) {
    start=a;
    stop=b;
    depth=c;
    gq=d;
  };

  uint32_t start,stop;
  int32_t depth,gq;
};

class variantRow {
 public:
  int  clear();
  variantRow(int nsample);
  ~variantRow();
  int n,n_allele;
  int32_t *dp,*dpf,*ad,*gt,*pf;
  int32_t  *gq;
};

class aggReader {
 public:
  aggReader(const vector<string>& input_file,const string &region);
  ~aggReader();
  int writeVcf(const char *output_file,char *output_type,int n_threads=0);//writes out a (homref called) vcf/bcf
  //  int writeDepthMatrix();//prints a depth matrix to stdout.
  int _nsample,_nsample2;

 private:
  int setSites(string region);
  int openFiles();
  int next();//moves forward one variant. dp_buf is updated accordingly.
  uint32_t dp_pos,var_start,var_stop;//start and end+1 positions of current variant. ie. a snp at pos 100 has var_start=100 and var_end=101; 
  int var_type;//0:snp 1:insertion 2:deletion
  //handles the depth intervals
  vector< deque< depthInterval > > dp_buf;//stores the depth intervals for nsamples.
  int syncBuffer();
  int setDepth();
  bcf_hdr_t *dp_hdr;
  bcf_srs_t *dp_rdr;    
  vector<bcf1_t *> dp_line;       //dp line
  void annotate_line();

  //handles the variants
  bcf1_t *out_line;//output line
  vector<bcf1_t *> line;       //current bcf record.
  int32_t *work;
  bcf_hdr_t *out_hdr;
  bcf_srs_t *var_rdr;    
  int nreader;
  variantRow *vr;
  int32_t*_dp,*_out_dp,ndp;
  int32_t*_gq,*_out_gq,ngq;
  int _current_rid;//current chromosome (bcf key)
  int dp_chr;
  int line_count;
  int moveDepthForward();
  int interval_start,interval_end;//the region specified.
};
