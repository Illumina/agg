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

extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "htslib/hts.h"

//this is a dummy definition for a function pilfered directly from bcftools.
int main_vcfmerge(int argc, char *argv[],char *file_list, char *output_fname,int nthreads);
}


using namespace std;

class depthInterval {
public:
  depthInterval() {
  };
  depthInterval(uint32_t a,uint32_t b,int c,int d) {
    start=a;
    stop=b;
    depth=c;
    gq=d;
  };

  //considering using 8bit depths, but htslib does not clearly support it yet.
  // void  setDepth(int x) {
  //   if(x<0)
  //     depth=0;
  //   else if(x>UINT8_MAX)
  //     depth=UINT8_MAX;
  //   else
  //     depth=x;
  // }

  uint32_t start,stop,depth,gq;
};

class depthReader {
public:
  depthReader(const char *depth_fname);
  ~depthReader();
  int next();
  int nread,nsample;
  depthInterval line;
  int chrom;
  bool open;
private:
  gzFile fp;
  int buf[5];
};
 

depthReader::depthReader(const char *depth_fname) {
  //  cerr<<"Opening "<<depth_fname<<endl;
  fp = gzopen(depth_fname, "r");
  open=true;
}

depthReader::~depthReader() {
  gzclose(fp);
}

int depthReader::next() {

  if(gzread(fp,buf,20)) {
    chrom=buf[0];
    line.start=buf[1];
    line.stop=buf[2];
    line.depth=buf[3];
    line.gq=buf[4];
    return(1);
  }
  else {
    open=false;
    return(0);
  }
}


class depthMerger {
  
public:
  int n,nsample;
  int32_t *dp,*gq;
  depthMerger(vector<string> & files);
  ~depthMerger();
  int writeDepthMatrix(const char *output_file,int nthreads=0);
  int syncBuffer();
  int findCurrPos();
  bcf_hdr_t *makeDepthHeader();
private:
  int cur_pos;
  int curr_chrom;
  vector< deque< depthInterval > > dp_buf;//stores the depth intervals for nsamples.
  vector< depthReader* > r;
  bcf_srs_t *sr;
  int _nfile;//how many input fiels are there?
  vector<string> _files;
  bcf_hdr_t *_hdr;
  bool _eof_warn;
};

bcf_hdr_t *depthMerger::makeDepthHeader() {
  _hdr = bcf_hdr_init("w");
  bcf_hdr_append(_hdr, "##source=agg");
  bcf_hdr_append(_hdr,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">");
  bcf_hdr_append(_hdr,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");

  bcf_hdr_t *src_hdr;
  htsFile *hfp;
  for(int i=0;i<_nfile;i++) {
    bcf_hdr_t *src_hdr = sr->readers[i].header;
    if(i==0) copyContigs(src_hdr,_hdr);
    if(bcf_hdr_nsamples(src_hdr)!=1) {
      cerr<<"ERROR: file "<<sr->readers[i].fname<<" had >1 sample.  This cannot be output from agg ingest1!"<<endl;
      exit(1);
    }
    char *sample_name=src_hdr->samples[0];
    if ( bcf_hdr_id2int(_hdr, BCF_DT_SAMPLE, sample_name)!=-1 )
      die("duplicate sample names. use --force-samples if you want to merge anyway");
    bcf_hdr_add_sample(_hdr,sample_name);      
  }
  bcf_hdr_sync(_hdr);
  cerr<<  bcf_hdr_nsamples(_hdr) << " samples in merged chunk"<<endl;
  return(_hdr);
}

depthMerger::~depthMerger() {
  bcf_sr_destroy(sr);
   
}

int depthMerger::syncBuffer() {
  int buf_size=1000;//maximum number of intervals to buffer.

  for(int i=0;i<nsample;i++) {
    dp[i]=bcf_int32_missing; //default.
    gq[i]=bcf_int32_missing; //default.
    if(r[i]->chrom==curr_chrom) {      
      //add to buffer until we pass cur_pos
      while(r[i]->chrom<=curr_chrom && (dp_buf[i].size() < buf_size || dp_buf[i].back().stop < cur_pos)) {
	if(r[i]->next() && r[i]->chrom==curr_chrom)
	  dp_buf[i].push_back( depthInterval(r[i]->line));//pretty sure we dont need to create a new deptinterval here?
	else
	  break;
      }
    }
    //discard intervals behind pos
    while(!dp_buf[i].empty() && dp_buf[i].front().stop < cur_pos) 
      dp_buf[i].pop_front();	      
    
    //update depth array.
    if(!dp_buf[i].empty() && dp_buf[i].front().start <= cur_pos && dp_buf[i].front().stop >= cur_pos){
      dp[i] = dp_buf[i].front().depth;
      gq[i] = dp_buf[i].front().gq;
    }
    //n   cerr << "buf_size["<<i<<"] = "<<dp_buf[i].size()<<endl;
  }

  //check if we have reached end of chromosome
  bool all_empty=true;
  for(int i=0;i<nsample;i++)
    if(!dp_buf[i].empty()||(r[i]->chrom==curr_chrom && r[i]->open))
      all_empty=false;


  //if chromosome is finished, move to next one
  if(all_empty) {

    // cerr << "ALL EMPTY pos="<<cur_pos+1<<endl;
    // for(int i=0;i<nsample;i++)
    //   cerr<< r[i]->chrom << "/" << r[i]->open << " ";
    // cerr <<endl;
    
    int nopen=0;
    for(int i=0;i<nsample;i++)
      if(r[i]->open)
	nopen++;

    if(nopen==0) 
      return(0);
    cerr << "chromosome end   = " << bcf_hdr_id2name(_hdr,curr_chrom) <<":"<< cur_pos<<endl;//

    //sweeps for the new chrom (one with the smallest index)
    int prev_chrom=curr_chrom;
    curr_chrom=r[0]->chrom;
    for(int i=1;i<nsample;i++) 
      if(r[i]->chrom < curr_chrom && r[i]->open) 
	curr_chrom=r[i]->chrom;
    assert(curr_chrom!=prev_chrom);
    for(int i=0;i<nsample;i++)  {
      if(r[i]->open)  {
	if(r[i]->chrom == curr_chrom)
	  dp_buf[i].push_back( depthInterval(r[i]->line));
      }
      else {
	if(!_eof_warn) {
	  cerr << "WARNING: " << _files[i] << " is finished but "<<nopen<<" .dpt files claim to have more data!"<<endl;
	  _eof_warn=true;
	}
      }
    }

    findCurrPos();
    return(syncBuffer());
  }

  return(1);
}

depthMerger::depthMerger(vector<string> & files) {
  _files = files;
  _eof_warn=false;
  nsample=n=_nfile=files.size();


  sr=bcf_sr_init() ;
  sr->require_index=1;
  for(int i=0;i<_nfile;i++) {
    if(!bcf_sr_add_reader (sr, files[i].c_str())) {
      fprintf(stderr,"problem opening bcf");
      exit(1);
    }
  }

  cerr << "nsample="<<nsample<<endl;
  for(int i=0;i<_nfile;i++)  {
    string tmp = files[i];
    tmp = tmp.substr(0,tmp.size()-3) + "dpt";
    r.push_back(new depthReader(tmp.c_str()));
  }
  cerr << "opened "<<n<<" .dpt files for merging..."<<endl;    
  dp_buf.resize(n);
  for(int i=0;i<n;i++)  {
    assert(r[i]->next());
    dp_buf[i].push_back( depthInterval(r[i]->line) );
    //    dp_buf[i].push_back( depthInterval(r[i]->line.start,r[i]->line.stop,r[i]->line.depth) );
  }
  curr_chrom=r[0]->chrom;
  makeDepthHeader();
}

int depthMerger::findCurrPos() {
  cur_pos = -1;
  for(int i=0;i<nsample;i++)
    if(dp_buf[i].front().start < cur_pos)
      cur_pos = dp_buf[i].front().start;
  cerr << "chromosome start = " << bcf_hdr_id2name(_hdr,curr_chrom)<<":"<< cur_pos<<endl;
  return(cur_pos);
}

int depthMerger::writeDepthMatrix(const char *output_file,int nthreads) {
  cerr << "Writing out "<<output_file<<endl;
  dp = new int32_t[nsample];
  gq = new int32_t[nsample];
  findCurrPos();

  htsFile *out_fh =   hts_open(output_file,"wb1");
  if(nthreads>0)  hts_set_threads(out_fh,nthreads);

  bcf_hdr_write(out_fh, _hdr);
  bcf1_t *line = bcf_init();

  while(    syncBuffer() ) {
    bcf_clear1(line) ;
    line->rid = curr_chrom;
    line->pos = cur_pos;
    bcf_update_alleles_str(_hdr, line, "N,.");
    bcf_update_format_int32(_hdr, line,"DP",dp,nsample);
    bcf_update_format_int32(_hdr, line,"GQ",gq,nsample);
    bcf_write1(out_fh, _hdr, line) ;
    cur_pos++;
    if(cur_pos%10000000==0)    cerr<<bcf_hdr_id2name(_hdr,curr_chrom)<<":"<<cur_pos<<endl;
  }
  cerr << "finished."<<endl;
  hts_close(out_fh);
  bcf_hdr_destroy(_hdr);
  delete[] dp;
  delete[] gq;
  return(0);
}

static void usage(){
  fprintf(stderr, "\n");
  fprintf(stderr, "About:   merges single sample agg files into an agg chunk\n");
  fprintf(stderr, "Usage:   agg ingest2 <input1> <input2> ... <inputN> -o output_prefix\n");
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
  for(vector<string>::iterator it1=file_list.begin();it1!=file_list.end();it1++)
    cerr << *it1<<endl;

  char *output_bcf=(char *)malloc(strlen(output)+5);  strcat(strcpy(output_bcf,output),".bcf");
  cerr << "Merging variants..." <<output_bcf<<endl;
  main_vcfmerge(argc,argv,file_list_fname,output_bcf,n_threads);
  output_bcf=(char *)malloc(strlen(output)+5);  strcat(strcpy(output_bcf,output),".bcf");
  cerr << "Indexing " <<output_bcf<<endl;
  bcf_index_build(output_bcf, BCF_LIDX_SHIFT);


  depthMerger d(file_list);
  char *dp_out_fname=(char *)malloc(strlen(output)+5);
  strcat(strcpy(dp_out_fname,output),".dpt");
  d.writeDepthMatrix(dp_out_fname , n_threads);
  cerr << "Indexing " <<dp_out_fname<<endl;
  bcf_index_build(dp_out_fname, BCF_LIDX_SHIFT);


  return(0);
}
