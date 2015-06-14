#include <math.h>
#include <map>
#include <getopt.h>
#include <vector>

#include "utils.h"
#include "tabReader.h"

extern "C" {
#include "htslib/synced_bcf_reader.h"
}

#define DEBUG 0

using namespace std;

tabReader::tabReader(const string& bed_fname,const string& region) {
  nread=0;
  if ((tbx = tbx_index_load(bed_fname.c_str())) == 0) 
    die("problem loading index");

  if(region!="") {
    if ((itr = tbx_itr_querys(tbx, region.c_str())) == 0) 
      die("problem with region "+region);
  }
  else {
    if((itr = tbx_itr_queryi(tbx, HTS_IDX_START, 0, 0))==0)
      die("problem with region "+region);
  }
  if ((fp = bgzf_open(bed_fname.c_str(), "r")) == 0) 
    die("problem opneing bed");
  s.s = 0; s.l = s.m = 0;
  sampleid=bed_fname;
  open=true;
  if(DEBUG>0)  cerr << "opened "<< bed_fname << endl;
}

int tabReader::next() {
  if (tbx_bgzf_itr_next(fp, tbx, itr, &s) >= 0) {
    stringstream ss(s.s);
    ss >> line.chrom;
    ss >> line.a;
    ss >> line.b;
    if(!(ss >> line.dp)) line.dp=30;
    if(!(ss >> line.id)) line.id=sampleid;
    //    cerr << line.chrom<<":"<<line.a<<"-"<<line.b<<" "<<line.dp<<" "<<line.id<<endl;
    return(1);
  }
  else {
    open=false;
    return(0);
  }
}
  
tabReader::~tabReader() {
  if(DEBUG>0)cerr<< "destructing "<<sampleid<<endl;
  tbx_itr_destroy(itr);
  bgzf_close(fp);
  tbx_destroy(tbx);
}
