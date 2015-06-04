#include "sampleReader.h"
#define DEBUG 0


sampleReader::sampleReader(const string& vcf_filename,const string& bed_filename) {
    vcf_fname = vcf_filename.c_str();
    bed_fname = bed_filename.c_str();
    gt_arr = NULL;
    reader=NULL;
    mgt_arr=0;
    closed=true;
}

int sampleReader::close() {
  if(!closed) {
    if(reader!=NULL) bcf_sr_destroy(reader);
    bed.clear();
    closed=true;
  }
  return(0);
}

sampleReader::~sampleReader() {//close files etc.
  close();
}

int sampleReader::setSites(vector<marker> * site_list,const string & chromosome) {
  closed=false;
  g.init();
  sites = site_list;
  int nvar=(*site_list).size();
  reader=bcf_sr_init() ;
  reader->require_index = 1;
  string reg;
  if(nvar>0) {
    stringstream sstm;
    sstm<<chromosome<<":"<<(*site_list)[0].pos<<"-"<<((*site_list)[(*site_list).size()-1].pos+1);
    reg=sstm.str();
    //  if(DEBUG>0)  cerr << reg << endl;
    //check that we are okay to process
    if(bcf_sr_set_regions(reader,reg.c_str(),0)!=0)
      die("problem setting regions");
  }

  if(!bcf_sr_add_reader (reader, vcf_fname.c_str()))
    die("problem with "+(string)vcf_fname + " (is the file indexed?)");
  if(bcf_hdr_nsamples(reader->readers[0].header)!=1)
    die("nsample!=1 in "+string(vcf_fname)+" ie. this is not a single sample vcf");

  if(nvar==0) 
    return(-1);

  move_vcf_forward();
  readBed(reg);
  bed_idx=0;
  return(0);
}

int sampleReader::move_vcf_forward() {

  if (bcf_sr_next_line (reader)) {
    line = bcf_sr_get_line(reader,0);
    var_pos=line->pos;
    return(1);
  }
  else {
    var_pos=-1;
    close();
    return(0);
  }
}

int sampleReader::next(marker & m,int *gt,float *gq,int32_t *dp,int32_t *ad) {  //out1/2 are the alleles. m is the marker we are genotyping.

  while(bed_idx<bed.size() && bed[bed_idx].second < m.pos) bed_idx++;
  if(DEBUG>2)  cerr <<" "<< var_pos << " " << line->d.allele[0] <<","<<line->d.allele[1]<<" ";
  if(!(var_pos==-1 || var_pos>=m.pos)) {
    cerr << "\nsample "<<  reader->readers[0].header->samples[0] << " out of sync " << var_pos+1 << " < " << m.pos+1 <<endl;
    cerr << "This can occur when a sample has a variant not present in the sites list\n"<<endl;
    die("sync error");
  }
  if( var_pos==m.pos && !m.ref.compare(line->d.allele[0]) && !m.alt.compare(line->d.allele[1])  ) {//has varaint?
    g.update(line,reader);
    //    cerr<<  bcf_gt_allele(g.gt[0])<<"/" << bcf_gt_allele(g.gt[1])<<"\t";
    move_vcf_forward();
  }
  else if( m.pos >= bed[bed_idx].first && m.pos < bed[bed_idx].second ) {//homref block.
    //    cerr<<"0/0"<<"\t";
    g.homref();
  }
  else {
    //    cerr<<"./."<<"\t";
    g.missing();
  }

  gt[0]=g.gt[0];
  gt[1]=g.gt[1];
  *gq=g.gq[0];
  *dp=g.dp[0];
  ad[0] = g.ad[0];
  ad[1] = g.ad[1];
  return(1);
}


// ichr,ifrom,ito are 0-based;
// returns -1 on error, 0 if the line is a comment line, 1 on success
int _regions_parse_line(char *line, int ichr,int ifrom,int ito, char **chr,char **chr_end,int *from,int *to)
{
  *chr_end = NULL;

  if ( line[0]=='#' ) return 0;

  int k,l;    // index of the start and end column of the tab-delimited file
  if ( ifrom <= ito )
    k = ifrom, l = ito;
  else
    l = ifrom, k = ito;

  int i;
  char *se = line, *ss = NULL; // start and end
  char *tmp;
  for (i=0; i<=k && *se; i++)
    {
      ss = i==0 ? se++ : ++se;
      while (*se && *se!='\t') se++;
    }
  if ( i<=k ) return -1;
  if ( k==l )
    {
      *from = *to = strtol(ss, &tmp, 10);
      if ( tmp==ss ) return -1;
    }
  else
    {
      if ( k==ifrom )
	*from = strtol(ss, &tmp, 10);
      else
	*to = strtol(ss, &tmp, 10);
      if ( ss==tmp ) return -1;

      for (i=k; i<l && *se; i++)
        {
	  ss = ++se;
	  while (*se && *se!='\t') se++;
        }
      if ( i<l ) return -1;
      if ( k==ifrom )
	*to = strtol(ss, &tmp, 10);
      else
	*from = strtol(ss, &tmp, 10);
      if ( ss==tmp ) return -1;
    }

  ss = se = line;
  for (i=0; i<=ichr && *se; i++)
    {
      if ( i>0 ) ss = ++se;
      while (*se && *se!='\t') se++;
    }
  if ( i<=ichr ) return -1;
  *chr_end = se;
  *chr = ss;
  return 1;
}

int  sampleReader::readBed(const string & region) {
  tbx_t *tbx;
  BGZF *fp;
  kstring_t s;
  int i;
  int nread=0;
  if ((tbx = tbx_index_load(bed_fname.c_str())) == 0) die("problem loading index");
  if ((fp = bgzf_open(bed_fname.c_str(), "r")) == 0) die("problem opneing bed");
  int ichr=0, ifrom=1, ito=2,from,to;
  char *chr, *chr_end;
  s.s = 0; s.l = s.m = 0;
  for (i =2; i < 3; ++i) {
    hts_itr_t *itr;
    if ((itr = tbx_itr_querys(tbx, region.c_str())) == 0) die("problem with region "+region);
    while (tbx_bgzf_itr_next(fp, tbx, itr, &s) >= 0) {
      //      puts(s.s);
      _regions_parse_line(s.s, ichr,ifrom,abs(ito), &chr,&chr_end,&from,&to);
      bed.push_back(pair<uint32_t,uint32_t> (from,to));
    }
    tbx_itr_destroy(itr);
  }
  if(DEBUG>0) cerr<<"stored "<<bed.size()<<" intervals ("<<bed[0].first <<","<< bed[bed.size()-1].second<<")"<<endl;
  free(s.s);
  bgzf_close(fp);
  tbx_destroy(tbx);
  
  return(nread);
}

Genotype::Genotype() {
  gt=NULL;
  gq=NULL;
  dp=NULL;
  ad=NULL;
  //  init()
}

void Genotype::init() {
  gt = (int*)malloc(2*sizeof(int));
  ad = (int*)malloc(2*sizeof(int32_t));
  dp = (int*)malloc(sizeof(int32_t));
  gq = (float*)malloc(sizeof(float));
  ngt_arr=2;
}

Genotype::~Genotype() {
  free(gt);
  free(ad);
  free(dp);
  free(gq);
}

void Genotype::homref() {
  gt[0]=bcf_gt_unphased(0);
  gt[1]=bcf_gt_unphased(0);
  pass=1;
  dp[0]=bcf_int32_missing;
  ad[0]=bcf_int32_missing;
  ad[1]=bcf_int32_missing;
  bcf_float_set_missing(   gq[0] );
}

void Genotype::missing() {
  gt[0]=bcf_gt_missing;
  gt[1]=bcf_gt_missing;
  pass=0;
  dp[0]=bcf_int32_missing;
  ad[0]=bcf_int32_missing;ad[1]=bcf_int32_missing;
  bcf_float_set_missing(   gq[0] );
}

int Genotype::update(bcf1_t *line,bcf_srs_t *reader) {
  if(line->n_allele!=2)  {
    cerr << line->pos + 1 << " " <<  line->n_allele << endl;
    die("variant file incorrectly normalised (should be bi-allelic)");
  }

  bcf_hdr_t *hdr = reader->readers[0].header;
  int mgt=2;
  ret = bcf_get_genotypes(hdr, line, &gt, &mgt);
  assert(ret==2);
  pass = bcf_has_filter(hdr, line, ".");  

  //  if(bcf_gt_allele(gt[0])>0 || bcf_gt_allele(gt[1])>0) {
    bcf_unpack(line, BCF_UN_FMT);
    int nval=2;
    ret = bcf_get_format_int32(hdr, line, "DP", &dp, &nval);
    if(ret<=0)//"handle dpi"
      ret = bcf_get_format_int32(hdr, line, "DPI", &dp, &nval);  
    ret = bcf_get_format_float(hdr, line, "GQ", &gq, &nval); 
    if(ret!=1)   bcf_float_set_missing(   gq[0] );
    ret =  bcf_get_format_int32(hdr, line, "AD", &ad, &nval);    
    if(ret!=2) {
      ad[0]=0;
      ad[1]=0;
    }
    //  }

  return(0);
}
