#include "agg.h"
#include "version.h"
#include <omp.h>

int main(int argc, char **argv) {


  if(argc < 2) {
    cerr << "\nProgram:\tagg (aggregation tool for multiple gvcfs)" << endl;
    cerr << "Version:\t" << VERSION <<endl;
    cerr << "Contact:\tjoconnell@illumina.com\n" << endl;
    cerr << "Copyright (c) 2015, Illumina, Inc. All rights reserved. See LICENSE for further details.\n"<<endl;
    cerr << "Usage:\tagg <command> [options]\n" << endl;
    cerr << "Commands:" << endl;
    cerr << "\tingest1\t\tconverts gvcfs to input suitable for agg ingest2" << endl;
    cerr << "\tingest2\t\tuses output files from ingest1 to build an agg chunk" << endl;
    cerr << "\tgenotype\tgenotypes and merges agg chunks into a multi-sample bcf/vcf" << endl;
    return(1);
  }
  else if(((string)argv[1]) == "ingest2") {
    merge_main(argc, argv);
    //    die("merge not implemented");
  }
  else if(((string)argv[1]) == "ingest1") {
    //    ingest1(argc, argv);
    ingest_main(argc,argv);
  }
  //    count1(argc,argv);
  else if(((string)argv[1]) == "genotype") {
    view1(argc, argv);
  }
  else if(((string)argv[1]) == "merge") {
    die("agg merge is not yet implemented. you can (carefully) merge dpt/bcf files using bcftools merge");
  }
  else {
    cerr << "Invalid command: " << argv[1] << endl;
  }
}
