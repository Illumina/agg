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

int power(char *vcf1,char *vcf2) {
  bcf_srs_t *sr =  bcf_sr_init() ;  
  sr->require_index=1;
  if(!bcf_sr_add_reader (sr, vcf1))    die("Problem opening vcf1");
  if(!bcf_sr_add_reader (sr, vcf2))    die("Problem opening vcf2");
  vector<int> qual;
  for(int i=0;i<=7;i++) qual.push_back(i*5);
  map<float, qbin> bins;

  bcf1_t *line1,*line2;
  float *af=NULL;

  int naf=0;
  while(bcf_sr_next_line (sr)) { 
    if( bcf_sr_has_line(sr,1) ) {

      if( bcf_sr_has_line(sr,0) )  {
	line1 = bcf_sr_get_line(sr,0) ;
	assert(line1->n_allele==2);
      }

      //      int ret=bcf_get_info_float(sr->readers[0].header,line2,tag,&af,&naf)
      line2 = bcf_sr_get_line(sr,1) ;
      assert(line2->n_allele==2);
      int ret=bcf_get_info_float(sr->readers[1].header,line2,"EUR_AF",&af,&naf);
      //      cerr<<line2->pos<<" "<<ret<<" "<<af[0]<<endl;

      assert(ret==1);
      if(!bins.count(af[0])) 
	bins[af[0]] =  qbin(qual.size());
      else {
	bins[af[0]]._nvar++;    
      }
      if( bcf_sr_has_line(sr,0) ) {
	for(int i=0;i<qual.size();i++)
	  if(line1->qual>=qual[i])
	    bins[af[0]]._counts[i]++;	    
      }
    }
  }

  for( map<float,qbin >::iterator it1=bins.begin();it1!=bins.end();it1++ )  {
    for(int i=0;i<qual.size();i++) 
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
  fprintf(stderr, "    -2, --vcf2      list of high quality variants from population requencing study\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(1);
}

int power1(int argc,char **argv) {
  int c;
  static struct option loptions[] =    {
    {"vcf1",required_argument,NULL,'1'},
    {"vcf2",required_argument,NULL,'2'},
    {0,0,0,0}
  };
  char *vcf1=NULL,*vcf2=NULL;
  if(argc<2) usage();
  while ((c = getopt_long(argc, argv, "1:2:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case '1': vcf1 = optarg; break;
      case '2': vcf2 = optarg; break;
      default: usage();
      }
  }
  optind++;
  if(!vcf1&&!vcf2)
    die("missing arg");
  power(vcf1,vcf2);
  return(0);
}
