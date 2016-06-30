#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <getopt.h>
#include <math.h>
#include <map>
#include <algorithm>
#include <vector>
#include "utils.h"
extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "filter.h"
}

using namespace std;


inline bool transition(char r,char a) {
  return ( (r=='A'&&a=='G') || (r=='G'&&a=='A') || (r=='C'&&a=='T') || (r=='T'&&a=='C') );
}


class FilterEvaluator {
public:
  FilterEvaluator(const string & vcf1,string annotation,const string & include);
  int setQuantile(int nquantile) {_nquant=nquantile;}
  int evaluate();
private:
  bcf_srs_t *_sr;
  bcf_hdr_t *_hdr;
  string _vcf;
  filter_t *_filter;
  int _nquant;
  const char *_annotation;
  vector< pair<float,vector<int> > > _values_snps; //vector of mendel counts.
  vector< pair<float,vector<int> > > _values_indels; //same again but for indels.
  int  printQuantiles(  vector< pair<float,vector<int> > >  & values,const string & variant_type);
};

  
FilterEvaluator::FilterEvaluator(const string & vcf1,string annotation,const string & include) {
  _sr =  bcf_sr_init() ;      
  _vcf = vcf1;
  _nquant=1000;
  _annotation = annotation.c_str();
  if(!bcf_sr_add_reader (_sr, _vcf.c_str())) 
    die("Problem opening vcf1");

  _hdr=_sr->readers[0].header;
  _filter=NULL;
  if(!include.empty())    _filter = filter_init(_hdr, include.c_str());
}

bool pairCompare(const std::pair<double, bool>& firstElem, const std::pair<double, bool>& secondElem) {
  return firstElem.first < secondElem.first;
}

int  FilterEvaluator::printQuantiles(  vector< pair<float,vector<int> > >  & values,const string & variant_type) {
  sort(values.begin(),values.end());
  int w = values.size()/_nquant;
  int current_quantile=w;
  int n = values.begin()->second.size();
  vector<int> counts(n,0);
  for(int i=0;i<values.size();i++) {
    if(i==current_quantile || i+1==values.size()) {
      cout << variant_type<<"\t"<< i << "\t"<< values[i].first << "\t";
      for(int j=0;j<n;j++)
	cout<<counts[j]<<"\t";
      cout<<endl;
      counts.assign(n,0);
      current_quantile+=w;
    }
    for(int j=0;j<n;j++)
      counts[j] += values[i].second[j];
  }
  return(0);
}

int FilterEvaluator::evaluate() {
  bcf1_t *line;
  int32_t info_dp,ac,an;
  int nval;
  int32_t *trio = new int32_t[27];
  int maxn=0;
  float val;
  //trio[dad*9 + mum*3 + kid] where dad/mum/kid are repsective genotypes from {0,1,2};
  while(bcf_sr_next_line (_sr)) { 
    assert(bcf_sr_has_line(_sr,0));
    line = bcf_sr_get_line(_sr,0);
    bcf_unpack(line, BCF_UN_INFO);
    if(_filter==NULL||filter_test(_filter,line,NULL)) {
      nval=1;
      float *ptr=&val;
      if(_annotation==NULL)
	val=line->qual;
      else
	assert(  bcf_get_info_float(_hdr,line,_annotation,&ptr,&nval) == 1);
      nval=27;
      assert(      bcf_get_info_int32(_hdr,line,"TRIO",&trio,&nval) == 27);
      if(strlen(line->d.allele[0])==1&&strlen(line->d.allele[1])==1) {//snp
	_values_snps.push_back( pair<float,vector<int> >(val, vector<int>(29,0) ));
	for(int i=0;i<27;i++) 	  (_values_snps.end()-1)->second[i]+=trio[i];
	if(transition(line->d.allele[0][0],line->d.allele[1][0])) (_values_snps.end()-1)->second[27]++;
	else (_values_snps.end()-1)->second[28]++;
      }
      else {//indel
	_values_indels.push_back( pair<float,vector<int> >(val, vector<int>(29,0) ));
	for(int i=0;i<27;i++) 	  (_values_indels.end()-1)->second[i]+=trio[i];
      }
    }
  }
  bcf_sr_destroy(_sr);	

  //header
  cout << "type\tN\tVALUE\t";
  string gt[3] = {"RR","RA","AA"};	
  for(int i=0;i<3;i++) 
    for(int j=0;j<3;j++) 
      for(int k=0;k<3;k++) 
  	cout<<	gt[i]<<"_"<<gt[j]<<"_"<<gt[k]<<"\t";
  cout <<"transition\ttransversion" << endl;
  printQuantiles(_values_snps,"snp");
  printQuantiles(_values_indels,"indel");
  return(0);
}

static void usage(){
  fprintf(stderr, "\n"); 
  fprintf(stderr, "Usage:   ftool eval input.bcf\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -a, --annotation  annotation to stratify on (defaults to QUAL)\n");
  fprintf(stderr, "    -i, --include     bcftools style filter to apply\n");
  fprintf(stderr, "    -q, --quantile    number of quantiles(=1000)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(1);
}

int evaluate1(int argc,char **argv) {
  int c;
  static struct option loptions[] =    {
    {"annotation",0,NULL,'a'},
    {"include",required_argument,NULL,'i'},
    {0,0,0,0}
  };
  string vcf1,ann;
  char *include=NULL;
  int qua=1000;
  if(argc<2) usage();
  while ((c = getopt_long(argc, argv, "a:i:q:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'a': ann = optarg; break;
      case 'i': include = optarg; break;
      case 'q': qua = atoi(optarg); break;
      default: die("Unrecognised argument");
      }
  }
  optind++;
  vcf1=argv[optind];
  if(vcf1.empty())
    die("missing input");
  FilterEvaluator f(vcf1,ann,include);
  f.setQuantile(qua);
  f.evaluate();
  return(0);
}
