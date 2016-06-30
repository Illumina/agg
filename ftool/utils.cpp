#include "utils.h"

void die(string s) {
  cerr << "ERROR: "<<s<<"\nExiting..."<<endl;
  exit(1);
}
