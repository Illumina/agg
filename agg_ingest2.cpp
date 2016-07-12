#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <limits>
#include "agg.h"
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include "string.h"
#include "utils.h"
#include "depthMerger.h"

extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"

//this is a dummy definition for a function pilfered directly from bcftools.
int dummy_main_vcfmerge(int argc, char *argv[],char *file_list, char *output_fname,int nthreads);
}


using namespace std;


static void usage(){
  fprintf(stderr, "\n");
  fprintf(stderr, "About:   merges single sample agg files into an agg chunk\n");
  fprintf(stderr, "Usage:   agg ingest2 input1 [input2 [...]] -o output_prefix\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Required options:\n");
  fprintf(stderr, "    -o, --output <output_prefix>       agg will output output_prefix.bcf and output_prefix.dpt\n");
  fprintf(stderr, "Optional options:\n");
  fprintf(stderr, "    -@, --thread INT                   number of compression threads [0]\n");
  fprintf(stderr, "    -l, --list   files.txt             plain text file listing agg chunks to merge]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(1);
}


int merge_main(int argc,char **argv) {

  int c;
  char *output=NULL,*file_list_fname=NULL;
  int n_threads=0;
  if(argc<3) usage();
  static struct option loptions[] =    {
    {"output",1,0,'o'},
    {"thread",1,0,'@'},
    {"list",1,0,'l'},
    {0,0,0,0}
  };

  while ((c = getopt_long(argc, argv, "l:@:o:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'l': file_list_fname = optarg; break;
      case 'o': output = optarg; break;
      case '@': n_threads = atoi(optarg); break;
      default: die("Unknown argument:"+(string)optarg+"\n");
      }
  }
  if(!output)    
    die("the -o option is required");
  optind++;
  cerr << "output = " <<output<<endl;

  vector<string> file_list;
  if(file_list_fname) {
    if(optind!=argc)
      die("cannot use -l option AND list input files on the command line");
    else
      readTextFile(file_list_fname,file_list);
  }
  else {
    if(optind==argc)
      die("no input files provided!");
    else 
      for(int i=optind;i<argc;i++)
	file_list.push_back(argv[i]);
  }
  cerr << file_list.size() << " files to merge"<<endl;

  //merge variants.
  char *output_bcf=(char *)malloc(strlen(output)+5);  strcat(strcpy(output_bcf,output),".bcf");
  cerr << "Merging variants..." <<output_bcf<<endl;
  dummy_main_vcfmerge(argc,argv,file_list_fname,output_bcf,n_threads);
  output_bcf=(char *)malloc(strlen(output)+5);  strcat(strcpy(output_bcf,output),".bcf");
  cerr << "Indexing " <<output_bcf<<endl;
  bcf_index_build(output_bcf, BCF_LIDX_SHIFT);

  ///build the depth tract
  char *dp_out_fname=(char *)malloc(strlen(output)+5);
  strcat(strcpy(dp_out_fname,output),".dpt");
  depthMerger d(file_list);
  d.setThreads(n_threads);
  d.writeDepthMatrix(dp_out_fname);
  cerr << "Indexing " <<dp_out_fname<<endl;
  bcf_index_build(dp_out_fname, BCF_LIDX_SHIFT);

  return(0);
}
