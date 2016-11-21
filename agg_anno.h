#define __STDC_LIMIT_MACROS
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>
#include <math.h>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <vector>

#include "utils.h"
extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "filter.h"
}


class BinResidualiser
{
public:
    BinResidualiser();
    int initialise(int nbin, int nsample);
    int fit(vector<float> &x, vector<float> &y);
    float residual(float x, float y);

private:
    int _nbin, _nsample;
    vector<float> _median, _mad;
    vector<int> _bins;
};

class Residualiser
{
public:
    Residualiser(){};
    int fit(vector<float> &x, vector<float> &y);

    float residual(float x, float y);

private:
    float _median, _mad;
    vector<float> _beta;
};

//class reads in a vcf and annotates features in INFO.
//features can be used in downstream filtering
//1. normalised depth.
class Standardiser
{
public:
    Standardiser(const string &vcf1, const char *include);
    int estimateParameters();
    int standardise(char *output_file,char *output_type);
    
private:
    bcf_srs_t *_sr;
    bcf_hdr_t *_hdr;
    string _vcf;
    filter_t *_filter;
    BinResidualiser _dp_residuals, _ab_residuals;
    double _alpha, _beta, _p_dpf;
    int _nsample;
};
