#include "ftool.h"

using namespace std;

int main(int argc, char **argv) {


  if(argc < 2) {
    cerr << "\nProgram:\tftool (some post-hoc filtering for agg output)" << endl;
    //    cerr << "Version:\t" << VERSION <<endl;
    cerr << "Contact:\tjoconnell@illumina.com\n" << endl;
    cerr << "Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.\n"<<endl;
    cerr << "Usage:\tftool <command> [options]\n" << endl;
    cerr << "Commands:" << endl;
    cerr << "\nannotate        adds some INFO fields that are useful for filtering" << endl;
    cerr << "eval        tabulates mendel rates/tstv across quantiles of an INFO field" << endl;
    cerr << "power           looks at proportion of \"true\" variants maintained under filtering routines" <<endl << endl;
    return(1);
  }
  else if(((string)argv[1]) == "annotate") {
    annotate1(argc, argv);
  }
  else if(((string)argv[1]) == "eval") {
    evaluate1(argc,argv);
  }
  else if(((string)argv[1]) == "power") {
    power1(argc, argv);
  }
  else {
    cerr << "Invalid command: " << argv[1] << endl;
  }
}
