#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include <getopt.h>
#include <math.h>
#include <map>
#include <vector>
#include "utils.h"

extern "C" {
#include "htslib/synced_bcf_reader.h"
}


class qbin {
public:
  qbin() {
    _n=0;
    _nvar=0;
  }

  qbin(int n) {
    _n=n;
    _counts.assign(_n,0);
    _nvar=0;
  }
  vector<int> _counts;
  int _n,_nvar;
  
};

int getQuantiles(char *vcf,char *tag,vector<float> & quantiles)
{
  bcf_srs_t *sr =  bcf_sr_init() ;  
  sr->require_index=1;
  if(!bcf_sr_add_reader (sr, vcf))    die("Problem opening vcf1");
  while(bcf_sr_next_line (sr))     {
    
  }
  return(0);
}

int power(char *vcf1,char *vcf2,char *tag,float min_val,float max_val,int nquantile) {
  bcf_srs_t *sr =  bcf_sr_init() ;  
  sr->require_index=1;
  if(!bcf_sr_add_reader (sr, vcf1))    die("Problem opening vcf1");
  if(!bcf_sr_add_reader (sr, vcf2))    die("Problem opening vcf2");
  vector<float> qual;
  qual.push_back(-10e6);
  float v=min_val;
  float w=(float)(max_val-min_val)/(float)nquantile;
  for(int i=0;i<=nquantile;i++) {
    v=min_val+i*w;
    qual.push_back(v);
    cerr << i<<":"<<v<<endl;
  }
  qual.push_back(10e6);
  map<float, qbin> bins;

  bcf1_t *line1,*line2;
  float *af=NULL;

  int naf=0;
  while(bcf_sr_next_line (sr)) { 
    if( bcf_sr_has_line(sr,1) ) {
      int ret,nval=1;
      float value;

      if( bcf_sr_has_line(sr,0) )  {
	line1 = bcf_sr_get_line(sr,0) ;
	assert(line1->n_allele==2);
	value=line1->qual;
	float *ptr=&value;
	if(tag!=NULL)
	  assert(bcf_get_info_float(sr->readers[0].header,line1,tag,&ptr,&nval)==1);
      }
	
      line2 = bcf_sr_get_line(sr,1) ;
      assert(line2->n_allele==2);
      ret=bcf_get_info_float(sr->readers[1].header,line2,"AF",&af,&naf);
      //      cerr<<line2->pos<<" "<<ret<<" "<<af[0]<<endl;
      if(ret!=1) die("no AF tag");
      if(!bins.count(af[0])) 
	bins[af[0]] =  qbin(qual.size());
      else {
	bins[af[0]]._nvar++;    
      }
      if( bcf_sr_has_line(sr,0) ) {
	for(size_t i=0;i<qual.size();i++)
	  if(value<=qual[i])
	    bins[af[0]]._counts[i]++;	    
      }
    }
  }
  cout <<"AF\tVALUE\tVCF1\tVCF2"<<endl;
  for( map<float,qbin >::iterator it1=bins.begin();it1!=bins.end();it1++ )  {
    for(size_t i=0;i<qual.size();i++) 
      cout << it1->first <<"\t"<<qual[i]<<"\t"<<it1->second._counts[i]<<"\t"<<it1->second._nvar<<endl;
  }

  return(0);
}

static void usage(){
  fprintf(stderr, "\n");
  fprintf(stderr, "Usage:   ftool power -1 test.vcf.gz -2 truth.vcf.gz\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Required options:\n");
  fprintf(stderr, "    -1, --vcf1      list of nominal variants\n");
  fprintf(stderr, "    -2, --vcf2      list of high quality variants\n");
  fprintf(stderr, "    -a, --annotation  annotation to stratify on (defaults to QUAL)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(1);
}

int power1(int argc,char **argv) {
  int c;
  static struct option loptions[] =    {
    {"vcf1",required_argument,NULL,'1'},
    {"vcf2",required_argument,NULL,'2'},
    {"annotation",0,NULL,'a'},
    {0,0,0,0}
  };
  char *vcf1=NULL,*vcf2=NULL,*ann=NULL;
  if(argc<2) usage();
  while ((c = getopt_long(argc, argv, "1:2:a:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'a': ann = optarg; break;
      case '1': vcf1 = optarg; break;
      case '2': vcf2 = optarg; break;
      default: usage();
      }
  }
  optind++;
  if(!vcf1&&!vcf2)
    die("missing arg");
  power(vcf1,vcf2,ann,-6,14,10);
  return(0);
}
