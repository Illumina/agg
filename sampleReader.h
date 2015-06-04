#pragma once

//multi-file genome. this will solve all your problems.
#include "agg.h"
#include <omp.h>

extern "C" {
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/vcfutils.h>
}

//stores a FORMAT field for one sample/variant. builds this from a bcf1_t record
class Genotype {
 public:
  Genotype();
  ~Genotype();
  void init();
  void homref();//sets this as a passing homref genotype
  void missing();//set this as missing genotype

  int update(bcf1_t *line,bcf_srs_t *reader);
  int32_t *ad,*dp;
  int *gt;
  float *gq;
  bool pass;
  int ret,ngt_arr;//working variables
};

class sampleReader {
public:
    sampleReader(const string& vcf_filename,const  string& bed_filename);
    ~sampleReader();
    int setSites(vector<marker> * site_list, const string & chromosome);     //set list of sites we are genotyping.
    int next(marker & m,int *gt,float *gq,int32_t *dp,int32_t *ad);//fills the respective FORMAT fields for marker m
    bcf_srs_t *reader;
    int close();
    Genotype g;//stores the genotype information at the next ALT genotype.


private:
    bool closed;
    int var_pos;        //position of current variant. -1 if out of region;
    vector<marker> * sites;  
    vector< pair<uint32_t,uint32_t> > bed;
    int bed_idx;
    string vcf_fname, bed_fname;
    bcf1_t *line;       //current bcf record.
    int readBed(const string& region);     //dumps the bed region into memory
    int move_vcf_forward();         //moves the vcf forwards to the next variant.
    int ngt,mgt_arr, *gt_arr;       //htslib working varialbes.
};
