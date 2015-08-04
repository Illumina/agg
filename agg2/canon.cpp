#include "agg.h"
#include <vector>
#include <algorithm>

typedef std::pair<string,int> mypair;

bool comparator ( const mypair& l, const mypair& r) { 
  return l.first < r.first; 
}

int main(int argc,char **argv) {
  int bufsize=1000;
  bcf_srs_t *rdr = bcf_sr_init() ; 

  if(!(bcf_sr_add_reader (rdr, argv[1])))
    die("problem opening file");
  char *outfile = "-";
  bcf_hdr_t *out_hdr = rdr->readers[0].header; 
  htsFile *out = hts_open(outfile, "wb");
  if ( !out ) die("problem opening output");

  bcf_hdr_write(out, out_hdr);
  vector<mypair> tosort;
  bcf1_t **buf = new bcf1_t*[bufsize];

  if (!bcf_sr_next_line (rdr)) {
    die("Empty bcf file");
  }
  bcf1_t *line = rdr->readers[0].buffer[0];    
  bool okay=true;

  while (okay) {
    int buf_idx=0;
    buf[buf_idx] = bcf_dup(line);
    bcf_unpack(buf[buf_idx],BCF_UN_STR);
    buf_idx++;
    if(line->n_allele!=2) {
      cerr << endl << line->pos+1 <<endl;
      die("sites had a line with n_allele!=2");
    }
    bcf_sr_next_line (rdr);
    line = rdr->readers[0].buffer[0]; 
    while(buf[buf_idx-1]->pos == line->pos) {//unpacking.
      buf[buf_idx] = bcf_dup(line);
      bcf_unpack(buf[buf_idx], BCF_UN_STR);
      buf_idx++;
      if(!bcf_sr_next_line (rdr)) {
	okay = false ;
	break;
      }
      line = rdr->readers[0].buffer[0]; 
      assert(buf_idx < bufsize);
    }
    tosort.clear();
    for(int i = 0; i < buf_idx; i++) {
      string ref_alt = ((string)buf[i]->d.allele[0]+","+(string)buf[i]->d.allele[1]);
      mypair m(ref_alt, i);
      tosort.push_back(m);
      //      cerr << m.first << " "<<m.second <<endl;
    }
    sort(tosort.begin(), tosort.end());
    for(int i = 0; i < buf_idx; i++) {///flussshhh
      //      cerr << tosort[i].first << " "<<tosort[i].second <<endl;
      int idx = tosort[i].second;
      bcf_write1(out, out_hdr, buf[idx]);
      bcf_destroy(buf[idx]);
    }
    buf_idx=0;
  }
  hts_close(out);
  delete[] buf;
  return(0);
}
