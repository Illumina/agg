#include "agg.h"

int bcf_hdr_sync(bcf_hdr_t *h);

void bcf_hdr_merge(bcf_hdr_t *hw, const bcf_hdr_t *hr, const char *clash_prefix, int force_samples)
{
    // header lines
    int ret = bcf_hdr_combine(hw, hr);
    if ( ret!=0 ) die("Error occurred while merging the headers.\n");

    // samples
    int i;
    for (i=0; i<bcf_hdr_nsamples(hr); i++)
    {
        char *name = hr->samples[i];
        if ( bcf_hdr_id2int(hw, BCF_DT_SAMPLE, name)!=-1 )
        {
            // there is a sample with the same name
	  if ( !force_samples ) die("Error: Duplicate sample names  use --force-samples to proceed anyway.\n");

            int len = strlen(hr->samples[i]) + strlen(clash_prefix) + 1;
            name = (char*) malloc(sizeof(char)*(len+1));
            sprintf(name,"%s:%s",clash_prefix,hr->samples[i]);
            bcf_hdr_add_sample(hw,name);
            free(name);
        }
        else
            bcf_hdr_add_sample(hw,name);
    }
}

int main(int argc,char **argv) {
  assert(argc>2);
  bcf_srs_t *sr =  bcf_sr_init() ; //htslib synced reader.
  int nfile = argc-1;
  htsFile *out_file  = hts_open("-", "wu");
  bcf_hdr_t *  out_hdr = bcf_hdr_init("w");
  sr->collapse = COLLAPSE_NONE;
  sr->require_index=0;
  sr->streaming=1;
  for(int i=1;i<argc;i++) 
    if(!(bcf_sr_add_reader (sr, argv[i]))) 
      cerr<<"Problem opening "<<i<< " "<<argv[i]<<endl;
  
  for (int i=0; i<sr->nreaders; i++)  {
    char buf[10]; snprintf(buf,10,"%d",i+1);
    bcf_hdr_merge(out_hdr, sr->readers[i].header,buf,0);
  }
  bcf_hdr_write(out_file, out_hdr);

  int nsample=0;
  int *idx = new int[nfile+1];
  for(int i=0;i<nfile;i++){
    idx[i]=nsample;
    nsample+=bcf_hdr_nsamples(sr->readers[i].header);
  }
  idx[nfile]=nsample;
  cerr << nsample << " samples in output" <<endl;
  int *gt = new int[nsample * 2];
  int ngt;
  bcf1_t **lines = new bcf1_t*[nfile];
  for(int i=0;i<nfile;i++) lines[i] = bcf_init1();
  bcf1_t *out_line = bcf_init1();
  while(bcf_sr_next_line (sr)) { 
    for(int i=0;i<nfile;i++) {
      if(bcf_sr_has_line(sr,i)) {//fill out_line
	lines[i] = bcf_sr_get_line(sr,i);
	bcf_get_genotypes(sr->readers[i].header,lines[i],gt+(2*idx[i]),&ngt);
      }
      else //make empty
	for(int j=2*idx[i];j<(idx[i+1]*2);j++) gt[j]=bcf_gt_missing;
    }
    bcf_write1(out_file, out_hdr, out_line);
  }
  return(0);
}
