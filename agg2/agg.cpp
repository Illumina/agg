#include "agg.h"
#include "version.h"
#include <omp.h>

int main(int argc, char **argv) {


  if(argc < 2) {
    cerr << "\nProgram:\tagg " << VERSION << " (aggregation tool for multiple samples)" << endl;
    cerr << "Contact:\tjoconnell@illumina.com\n" << endl;
    cerr << "Copyright (c) 2015, Illumina, Inc. All rights reserved. See LICENSE.pdf for further details.\n"<<endl;
    cerr << "Usage:\tagg <command> [options]\n" << endl;
    cerr << "Commands:" << endl;
    cerr << "\tcollate\t\tinitialises the database's variant list" << endl;
    cerr << "\tupdate\t\tadds new samples to the database" << endl;
    cerr << "\tcount\t\tsummary statistics (genotype counts, passrate, etc)" << endl;
    cerr << "\tgenotype\tproduce multisample vcf from the database" << endl;
    //  ": aggregates variants and summary information from a list of gvcfs\nContact: joconnell@illumina.com\nUsage: agg file_list.txt output.vcf.gz [chrom]\n"<<endl;
    return(1);
  }
  else if(((string)argv[1]) == "collate") {
    collate1(argc, argv);
  }
  else if(((string)argv[1]) == "count") {
    count1(argc, argv);
  }
  //    count1(argc,argv);
  else if(((string)argv[1]) == "genotype") {
    genotype1(argc, argv);
  }
  else if(((string)argv[1]) == "update") {
    update1(argc, argv);
  }
  else {
    cerr << "Invalid command: " << argv[1] << endl;
  }
}
