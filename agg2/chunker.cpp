#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <limits>
#include <vector>
#include <string>
#include <iostream>
#include <getopt.h>
#include <iomanip>
#include <math.h>

extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "htslib/hts.h"
}

using namespace std;

inline void die(const string& s) {
    cerr << "ERROR: " << s << "\nExiting..." << endl;
    exit(1);
}

static void usage(){
  fprintf(stderr, "\n");
  fprintf(stderr, "About:   prints a sensible set of regions with roughly n variants per region \n");
  fprintf(stderr, "Usage:   chunker input.bcf\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Required options:\n");
  fprintf(stderr, "    -n, --nchunk <int>              number of chunks\n");
  fprintf(stderr, "\n");
  exit(1);
}

void printChunks(const string&chrom,vector<int> & pos,int n) {
  float nwin = floor( (float)pos.size() / (float)n);
  n = ceil((float)pos.size()/nwin);
    
  
  for( int i=0;i<pos.size();i+=n) {
    int a = pos[i]+1;
    int b;
    if(i+n<pos.size())    
      b = pos[i+n];
    else
      b = pos[pos.size()-1];
    cout << chrom << ":" << setw(9) << setfill('0') << a <<"-"<< setw(9) << setfill('0') <<b<<endl;
  }
};


int main(int argc, char **argv) {
  int c,nchunk=0;
  //  if(argc<3) usage();
  static struct option loptions[] =    {
    {"nchunk",1,0,'n'},
    {0,0,0,0}
  };

  while ((c = getopt_long(argc, argv, "n:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'n': nchunk = atoi(optarg); break;
      case '?': usage();
      default: die("Unknown argument:"+(string)optarg+"\n");
      }
  }

  if( nchunk<=0 ) 
    die("the --nchunk argument (-n) is required");

  string site_list = argv[optind];

  int32_t curr_chrom,prev_chrom;


  bcf_srs_t *rdr=bcf_sr_init() ;     
  if(!(bcf_sr_add_reader (rdr, site_list.c_str()) ))
    die("problem opening "+site_list);
  vector<int> positions;
  bcf1_t *line; 
  int nline=0;
  while(bcf_sr_next_line(rdr)) {
    line = bcf_sr_get_line(rdr,0);
    if(line->rid!=prev_chrom && nline>0) {//print chunks
      printChunks(bcf_hdr_id2name(rdr->readers[0].header,prev_chrom),positions,nchunk);
      prev_chrom=line->rid;
      positions.clear();
    }
    positions.push_back(line->pos);
    nline++;
    prev_chrom=line->rid;
  } 
   
  return(0);
}

