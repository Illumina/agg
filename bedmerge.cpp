#include <math.h>
#include <map>
#include <getopt.h>
#include <vector>

#include "utils.h"
#include "tabReader.h"

extern "C" {
#include "htslib/synced_bcf_reader.h"
}

using namespace std;


int merge(vector<string> & files,string region) {
  vector< tabReader* > r;
  int n=files.size();
  for(int i=0;i<n;i++) 
    r.push_back(new tabReader(files[i],region) );
  
  cerr << "opened "<<n<<" bed files for merging..."<<endl;


  for(int i=0;i<n;i++)
    r[i]->next();

  int nopen=n;
  string curr_chrom=r[0]->line.chrom;
  while(nopen>0) {
    int minidx=-1;
    nopen=0;
    for(int i=0;i<n;i++) {
      if(r[i]->open) {
	if(minidx==-1||r[i]->line.a<r[minidx]->line.a) {
	  if(r[i]->line.chrom==curr_chrom) 
	    minidx=i;
	}
	nopen++;
      }
    }
    if(nopen>0) {
      cout << r[minidx]->line.chrom<<"\t"<<r[minidx]->line.a<<"\t"<<r[minidx]->line.b<<"\t"<<r[minidx]->line.dp<<"\t"<<r[minidx]->line.id<<endl;
      r[minidx]->next();
    }
    bool change_chrom=true;
    for(int i=0;i<n;i++)  {
      if(curr_chrom==r[i]->line.chrom) {
	change_chrom=false;
	break;
      }
    }
    if(change_chrom)
      curr_chrom=r[0]->line.chrom;
  }
 
  for(int i=0;i<n;i++) delete(r[i]);
  cerr << "finished."<<endl;
  return(0);
}


int main(int argc,char **argv) {
  string reg="";
  int c;
  static struct option loptions[] =    {
    {"regions",1,0,'r'},
    {0,0,0,0}
  };
  
  while ((c = getopt_long(argc, argv, "r:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'r': reg=optarg; break;    
      case '?': die("lol");
      default: die("Unknown argument:"+(string)optarg+"\n");
      }
  }
  
  vector<string> files;
  while(optind<argc) files.push_back((string)argv[optind++]);   

  merge(files,reg);
  return(0);
}
