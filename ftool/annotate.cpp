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

using namespace std;

//calculates the median and the mad
//note this modifies x
int madian(vector<float> & x,float & median,float & mad) {
  assert(x.size()>1);
  sort(x.begin(),x.end());
  vector<float> r;
  r.resize(x.size());
  if(x.size()%2==1) 
    median = x[x.size()/2];
  else 
    median = ( (float)x[x.size()/2] + (float)x[x.size()/2+1] ) / 2;
  for(int i=0;i<x.size();i++)
    r[i] = fabs((float)(x[i]-median));
  sort(r.begin(),r.end()); 
  if(x.size()%2==1) 
    mad = r[r.size()/2];
  else 
    mad = ( (float)r[r.size()/2] + (float)r[r.size()/2+1] ) / 2;
  mad *= 1.4826;
  return(0);
}

class Residualiser  {
public:
  Residualiser() {  };
  int fit(vector<float> & x,vector<float> & y) {
    //model fitting.
    vector<float> r;
    madian(r,_median,_mad);
    cerr << "n="<<y.size()<<"\nmedian="<<_median<<"\nmad="<<_mad<<endl;
    return(0);
  }

  float residual(float x,float y) {
    float ret =  y - (_beta[0] + _beta[1]*x + _beta[2]*x*x);
    ret -=  _median;
    ret /= _mad;
    return(ret);
  }

private:
  float _median,_mad;
  vector<float> _beta;  
};

class BinResidualiser  {
public:
  BinResidualiser() {
    _nsample=0;
    _nbin=0;
  }

  int  initialise(int nbin,int nsample) {  
    _nsample=nsample;
    _bins.clear();
    _nbin=nbin;
    for(int i=0;i<_nbin;i++)      
      _bins.push_back( ceil(nsample * (float)i / (float)nbin ));
    _bins.push_back(nsample+1);
    //    for(int i=0;i<_nbin;i++)            cerr << i << " " << _bins[i]<<endl;
    return(0);
  };

  int fit(vector<float> & x,vector<float> & y) {
    assert(_nbin>0);
    vector<float> x1;
    x1.reserve(x.size()/_nbin);
    _median.clear();
    _mad.clear();
    _median.resize(_nbin);
    _mad.resize(_nbin);
    for(int b=0;b<_nbin;b++) {
      for(int i=0;i<x.size();i++) 
	if(x[i]>=_bins[b] && x[i]<_bins[b+1])
	  x1.push_back(y[i]);
      madian(x1,_median[b],_mad[b]);
      cerr <<"("<<_bins[b]<<","<<_bins[b+1]<<") "<< b<<" n="<<x1.size()<<"\nmedian="<<_median[b]<<"\nmad="<<_mad[b]<<endl;
      x1.clear();
    }
    return(0);
  }

  float residual(float x,float y) {
    int b=0;
    while(b<_nbin&&_bins[b]<=x) b++;
    //    cerr << x << " "<<b<<endl;
    b--;
    if(!(b>=0&&b<_nbin)) {
      cerr <<"x="<<x<<endl;
      cerr <<"b="<<b<<" _nbin="<<_nbin<<endl;
      exit(1);
    }
    return(    (y - _median[b])/_mad[b] );
  }

private:
  int _nbin,_nsample;
  vector<float> _median,_mad;
  vector<int> _bins;
};


//class reads in a vcf and annotates features in INFO.
//features can be used in downstream filtering
//1. normalised depth.
class Standardiser {
public:
  Standardiser(const string & vcf1,const string & include);
  int estimateParameters();
  int standardise();
private:
  bcf_srs_t *_sr;
  bcf_hdr_t *_hdr;
  string _vcf;
  filter_t *_filter;
  BinResidualiser _res;
  double _alpha,_beta,_p_dpf;
  int _nsample;
};
  
Standardiser::Standardiser(const string & vcf1,const string & include) {
  _sr =  bcf_sr_init() ;      
  _vcf = vcf1;

  if(!bcf_sr_add_reader (_sr, _vcf.c_str())) 
    die("Problem opening vcf1");

  _hdr=_sr->readers[0].header;
  _nsample =   bcf_hdr_nsamples(_hdr);
  vector<float> af_bins;
  for(int i=0;i<20;i++)
    af_bins.push_back( (float)i/20. );
  _filter=NULL;
  if(!include.empty())    _filter = filter_init(_hdr, include.c_str());
}

//first pass 
//1. read in (filtered) depth and AF
//2. regress DP ~ AF + ALT_COUNT^2
//3. find mean/sd of residuals.
int Standardiser::estimateParameters() {
  bcf1_t *line;
  int32_t info_dp,ac,an,dpf,dpa;
  int nval=1;
  vector<float> depth;
  vector<float> alt_count;

  // binomial_data dpa_dpf;
  // dpa_dpf.n=0;
  // dpa_dpf.k.reserve(100e6);
  // dpa_dpf.size.reserve(100e6);
  double n_dpf=0;
  double n_dpa=0;
  depth.reserve(100e6);
  alt_count.reserve(100e6);
  int maxn=0;
  while(bcf_sr_next_line (_sr)) { 
    assert(bcf_sr_has_line(_sr,0));  
    line = bcf_sr_get_line(_sr,0) ;
    bcf_unpack(line, BCF_UN_INFO);
    if(_filter==NULL||filter_test(_filter,line,NULL)) {
      int32_t *ptr=&info_dp;
      assert(      bcf_get_info_int32(_hdr,line,"DP",&ptr,&nval) == 1);
      ptr=&ac;
      assert(      bcf_get_info_int32(_hdr,line,"AC",&ptr,&nval) == 1);
      ptr=&an;
      assert(      bcf_get_info_int32(_hdr,line,"AN",&ptr,&nval) == 1);
      ptr=&dpf;
      assert(      bcf_get_info_int32(_hdr,line,"DPF",&ptr,&nval) == 1);
      ptr=&dpa;
      assert(      bcf_get_info_int32(_hdr,line,"DPA",&ptr,&nval) == 1);

      //if an==0 then it is not a valid site.
      if(an>0) {
	alt_count.push_back( (float)ac );
	depth.push_back((float)info_dp);
	// dpa_dpf.k.push_back(dpa);
	// dpa_dpf.size.push_back(dpf+dpa);
	// dpa_dpf.n++;
	n_dpf+=dpf;
	n_dpa+=dpa;
      }
      if(an>maxn) 
	maxn=an;
    }
  }
  _p_dpf = n_dpf/(n_dpf+n_dpa);
  cerr <<"_p_dpf="<<_p_dpf<<endl;
  //  fit_betab(&dpa_dpf,_alpha,_beta);
  bcf_sr_destroy(_sr);	
  _res.initialise(20,maxn);
  _res.fit(alt_count,depth);

  return(0);
}

//second pass.
//re-open file. append INFO/S_DP - our normalised depth measure.
int Standardiser::standardise() { 
  bcf1_t *line;
  int32_t info_dp,ac,an,dpf,dpa;
  int nval=1;
  _sr =  bcf_sr_init() ;  
  if(!bcf_sr_add_reader (_sr, _vcf.c_str()))    
    die("Problem opening vcf1");
  _hdr=_sr->readers[0].header;
  bcf_hdr_append(_hdr, "##INFO=<ID=S_DP,Number=1,Type=Float,Description=\"normalised depth (mean 0 and sd 1)\">");
  bcf_hdr_append(_hdr, "##INFO=<ID=S_DPF,Number=1,Type=Float,Description=\"-log10(p-value) that DPA/(DPA+DPF) is extremely low. Assumes beta-binomial distribution \">");

  htsFile *out_fh  = hts_open("-", "wv");
  bcf_hdr_write(out_fh, _hdr);
  while(bcf_sr_next_line (_sr)) { 
    assert(bcf_sr_has_line(_sr,0));  
    line = bcf_sr_get_line(_sr,0) ;
    bcf_unpack(line, BCF_UN_INFO);    
    int32_t *ptr=&info_dp;
    assert(      bcf_get_info_int32(_hdr,line,"DP",&ptr,&nval) == 1);
    ptr=&ac;
    assert(      bcf_get_info_int32(_hdr,line,"AC",&ptr,&nval) == 1);
    ptr=&an;    assert(      bcf_get_info_int32(_hdr,line,"AN",&ptr,&nval) == 1);
    ptr=&dpf;
    assert(      bcf_get_info_int32(_hdr,line,"DPF",&ptr,&nval) == 1);
    ptr=&dpa;
    assert(      bcf_get_info_int32(_hdr,line,"DPA",&ptr,&nval) == 1);

    float s_dpf = 0;
    if(dpa>0&&dpf>0)      {
      //s_dpf = -log10(pbetab(dpa,dpa+dpf,_alpha,_beta,false,1000));
      s_dpf = (float)dpf/(float)(dpa+dpf);      
      s_dpf = ( s_dpf - _p_dpf ) / sqrt ( _p_dpf*(1-_p_dpf)/(float)(dpa+dpf) );
    }

    //    cerr << line->pos+1<<" "<<dpa<<"/"<<dpf<<" = "<<p_dpf<<endl;
    float x = (float)ac;
    if(an==0) x=0.;
    float ndp = _res.residual(x,(float)info_dp);
    //    float call = (float)an/(2. * (float)_nsample);
    bcf_update_info_float(_hdr, line, "S_DP", &ndp, 1);
    //    bcf_update_info_float(_hdr, line, "CALL", &call, 1);
    bcf_update_info_float(_hdr, line, "S_DPF", &s_dpf, 1);
    //    bcf_update_info_float(_hdr, line, "CALLRATE", &callrate, 1);
    bcf_write1(out_fh, _hdr, line) ;
    bcf_clear1(line) ;
  }
  hts_close(out_fh);
  bcf_sr_destroy(_sr);	  
  return(0);
}


static void usage(){
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   ftool annotate input.bcf\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -i, --include      filters to apply eg. -i 'QUAL>=10 && DP<100000 && HWE<10' \n");  
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(1);
}

int annotate1(int argc,char **argv) {
  int c;
  static struct option loptions[] =    {
    {"include",required_argument,NULL,'i'},
    {0,0,0,0}
  };
  char *vcf1=NULL,*include=NULL;
  if(argc<2) usage();
  while ((c = getopt_long(argc, argv, "i:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'i': include = optarg; break;
      default: usage();
      }
  }
  optind++;
  if(optind==argc)
    usage();
  vcf1=argv[optind];  
  if(!vcf1)
    die("missing input");
  Standardiser s(vcf1,include);
  s.estimateParameters();
  s.standardise();
  return(0);
}
