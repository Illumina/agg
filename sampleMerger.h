#pragma once

//multi-file genome. this will solve all your problems.
#include "agg.h"
#include "sampleReader.h"
#include <omp.h>

extern "C" {
#include <htslib/bgzf.h>
#include <htslib/tbx.h>
#include <htslib/vcfutils.h>
}


class sampleMerger {
public:
    sampleMerger(const string& database_location);
    ~sampleMerger();
    int setSites(vector<marker> * site_list,const string& chromosome);     //list of sites we are genotyping.
    int siteSummary(vector<marker> * site_list,const string& chromosome);     
    //  int next(marker & m, vector<byte> & out);
    int next();
    int nsample, nsample_in_db;      //samples being included, total samples in database
    int writeVcf(char *out_filename,int output_type);
    int setSamples(char *samples,bool isfile);      //subset samples.
    string dbdir;
    int countGenotypes(vector<int> & counts);

private:
    void    wipeRow();
    bool samples_initialised;
    int32_t *new_gts,*new_ad,*new_dp;       //stores a line of genotypes.
    float *new_gq;
    vector<string> sampleids;
    vector<bool> included;      //vector of nsample bools telling if we should include this sample or not.
    vector<sampleReader> genomes;
    vector<marker> *sites;
    string chrom;  
    int marker_idx;
};

