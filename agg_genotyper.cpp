#include "agg_genotyper.h"

#define DEBUG 0

int fillHeader(bcf_hdr_t *hdr) {//fills in the standard stuff for an agg header.

  bcf_hdr_append(hdr, ((string)"##source=agg "+VERSION).c_str());
  bcf_hdr_append(hdr, "##INFO=<ID=PF,Number=A,Type=Float,Description=\"proportion of genotypes containing an ALT that passed the original single sample gvcf filter\">");
  bcf_hdr_append(hdr, "##INFO=<ID=GN,Number=G,Type=Integer,Description=\"count of each genotype.\">"); 
  bcf_hdr_append(hdr, "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"sum of allele depths for ALL individuals\">"); 
  bcf_hdr_append(hdr, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"sum of depth  across all samples\">");
  bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
  bcf_hdr_append(hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Filtered basecall depth used for site genotyping\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=DPF,Number=1,Type=Integer,Description=\"Basecalls filtered from input prior to site genotyping\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed.\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
  bcf_hdr_append(hdr, "##FORMAT=<ID=PF,Number=A,Type=Integer,Description=\"variant was PASS filter in original sample gvcf\">");
  return(0);
}


int variantRow::clear() {
  for(int i=0;i<n;i++) {
    pf[i]= bcf_int32_missing;
    dp[i]= bcf_int32_missing;
    dpf[i]= bcf_int32_missing;
    gt[2*i]= bcf_int32_missing;
    gt[2*i+1]= bcf_int32_missing;
    ad[2*i]= bcf_int32_missing;
    ad[2*i+1]= bcf_int32_missing;
    gq[i]= bcf_int32_missing;
  }
  return(0);
}

variantRow::variantRow(int nsample) {
  n_allele=2;
  n=nsample;
  pf = (int32_t *)malloc(n*sizeof(int32_t));
  dp = (int32_t *)malloc(n*sizeof(int32_t));
  dpf = (int32_t *)malloc(n*sizeof(int32_t));
  gt = (int32_t *)malloc(2*n*sizeof(int32_t));
  gq = (int32_t *)malloc(n*sizeof(int32_t));
  ad = (int32_t *)malloc(2*n*sizeof(int32_t));
}

variantRow::~variantRow() {
  free(dp);
  free(ad);
  free(gq);
  free(gt);
}

bcf_srs_t *vcf_ropen(const vector<string>& input_files,const string &region) {
  bcf_srs_t *sr=bcf_sr_init() ;
  sr->require_index = 1;
  sr->collapse = COLLAPSE_NONE;

  if(region!="")
    if(bcf_sr_set_regions(sr,region.c_str(),0)!=0)
      die("problem setting regions");

  for(int i=0;i<input_files.size();i++)
    if(!bcf_sr_add_reader (sr, (input_files[i]).c_str()))
      die("problem with "+input_files[i] + " (is the file indexed?)");
  
  return(sr);
}

aggReader::aggReader(const vector<string>& input_files,const string &region) {
  nreader = input_files.size();
  _nsample = 0;  
  vector<string> tmp(nreader);
  for(int i=0;i<nreader;i++)    tmp[i]=input_files[i];
  var_rdr = vcf_ropen(tmp,region);


  if(region.find(":")<region.size()) {//interval specified, we need to buffer behind the start due to indel issues.
    string tmp1 = region.substr(region.find(":")+1);
    if(tmp1.find("-")<tmp1.size()) {
      cerr << tmp1.substr(0,tmp1.find("-"))<<endl;
      interval_start = stoi( tmp1.substr(0,tmp1.find("-")))-1;
      interval_end = stoi(tmp1.substr(tmp1.find("-")+1))-1;
    }
    else {
      interval_start=stoi(tmp1)-1;
      interval_end=interval_start;
    }
    interval_start=max(0,interval_start);
    cerr << "interval_start = "<<interval_start<<" interval_end="<<interval_end<<endl;
  } else {
    interval_start=0;
    interval_end=INT32_MAX;
  }

  for(int i=0;i<nreader;i++)    tmp[i]=input_files[i].substr(0,input_files[i].size()-4)+".dpt"  ;
  dp_rdr = vcf_ropen(tmp,region);

  for(int i=0;i<nreader;i++)
    assert(bcf_hdr_nsamples(var_rdr->readers[i].header) == bcf_hdr_nsamples(dp_rdr->readers[i].header));

  for(int i=0;i<nreader;i++)    _nsample += bcf_hdr_nsamples(var_rdr->readers[i].header);
  dp_buf.resize(_nsample);
  ndp=0;

  _nsample2 = 2*_nsample;
  _dp = (int32_t *)malloc(_nsample*sizeof(int32_t));
  _gq = (int32_t *)malloc(_nsample*sizeof(int32_t));
  _out_dp = (int32_t *)malloc(_nsample*sizeof(int32_t));
  _out_gq = (int32_t *)malloc(_nsample*sizeof(int32_t));
  for(int i=0;i<_nsample;i++) {
    _dp[i] = -1;
    _gq[i] = -1;
    _out_dp[i] = -1;
    _out_gq[i] = -1;
  }
  vr = new variantRow(_nsample);
  line.resize(nreader);
  dp_line.resize(nreader);
  cerr << _nsample<< " samples"<<endl;
  dp_pos=-1;
  moveDepthForward();
  line_count=0;
}

aggReader::~aggReader() {
  free(_out_dp);
  free(_out_gq);
  free(_dp);
  free(_gq);
  bcf_sr_destroy(dp_rdr);
  bcf_sr_destroy(var_rdr);
  delete vr;
}

int aggReader::moveDepthForward() {
  int32_t*dp_ptr=_dp;
  int32_t*gq_ptr=_gq;
  if(  bcf_sr_next_line (dp_rdr) ) {
    for(int j=0;j<nreader;j++) {
      int ntmp=bcf_hdr_nsamples(dp_rdr->readers[j].header);
      if( bcf_sr_has_line(dp_rdr,j) ) {
	dp_line[j] = bcf_sr_get_line(dp_rdr,j);	  
	ngq=ndp=bcf_hdr_nsamples(dp_rdr->readers[j].header);
	if(	bcf_get_format_int32(dp_rdr->readers[j].header, dp_line[j], "DP",&dp_ptr , &ndp) <= 0)
	  die("Problem with a .dpt file at position "+dp_line[j]->pos+1);
	if(	bcf_get_format_int32(dp_rdr->readers[j].header, dp_line[j], "GQ",&gq_ptr , &ngq) <= 0)	  
	  die("Problem with a .dpt file at position "+dp_line[j]->pos+1);//change this to gq

	dp_pos = dp_line[j]->pos;
	dp_chr = dp_line[j]->rid;
      } 
      else {//no line? all DP=missing
	for(int k=0;k<ntmp;k++)  {
	  dp_ptr[k]=bcf_int32_missing;
	  gq_ptr[k]=bcf_int32_missing;
	}
      }
      dp_ptr+=ntmp;
      gq_ptr+=ntmp;
    }
    //    cerr << "moveDepthForward: "<<dp_pos+1 <<"\t"<<dp[781]<<":"<<gq[781]<<endl;;
    return(1);
  }
  else
    return(0);
}


//this ensures dp_buf contains all the intervals covering the current position
int aggReader::syncBuffer() {
  int nsync=0;//number of samples in sync with var_pos
  if(DEBUG>1)  cerr << "var_pos = (" <<var_start<<","<<var_stop<<")  var_type="<<var_type<<endl;
  
  bool dp_open = true;
  assert(dp_buf.size()==_nsample);
  while(nsync!=_nsample) {
    nsync=0;
    for(int i=0;i<_nsample;i++) {

      if(!dp_buf[i].empty()) {
	//our dp interval contains the variant begin genotyped.
	if(dp_buf[i].front().stop >= var_start && dp_buf[i].back().stop >= var_stop)
	  nsync++;
	else
	  while(!dp_buf[i].empty() && dp_buf[i].front().stop < var_start)//discard intervals behind var_start
	    dp_buf[i].pop_front();
      }
    }

    if(dp_pos <= var_stop && dp_chr==cur_chr && dp_open) {//read in intervals up to end of variant position or chromosome ends.
      if(dp_pos>=var_start) { //are we interested in depth here?
	for(int i=0;i<_nsample;i++) { 
	  if(!dp_buf[i].empty() && _dp[i]==dp_buf[i].back().depth && _gq[i]==dp_buf[i].back().gq)
	    dp_buf[i].back().stop++;
	  else 
	    dp_buf[i].push_back(depthInterval(dp_pos,dp_pos,_dp[i],_gq[i]));
	}
      }   
      dp_open=moveDepthForward();
      if(DEBUG>1) cerr << "dp_pos="<<dp_pos<<endl;
    }
    if(dp_chr!=cur_chr || !dp_open) break;
  }

  if(DEBUG>3) 
    for(int i=0;i<_nsample;i++)
      if(!dp_buf[i].empty())   
	cerr << "sample"<<i<<" dp=("<<dp_buf[i].front().start<<","<<dp_buf[i].back().stop<<")"<<endl;  

  if(DEBUG>1) {
    float sum=0.;
    unsigned    int maxsize=0;
    for(int i=0;i<_nsample;i++)  {
      sum+=(float) dp_buf[i].size();
      if(dp_buf[i].size()>maxsize)
	maxsize=dp_buf[i].size();
    }
    cerr << "Mean buffer size = " << sum/_nsample <<"  Max size = "<< maxsize<< endl;
  }

  return(1);
}

int aggReader::setDepth() {
  if(DEBUG>1) cerr << "(var_start,var_stop) = ("<<var_start<<","<<var_stop<<")"<<endl;
  syncBuffer();
  for(int i=0;i<_nsample;i++) {
    deque<depthInterval>::iterator it1 = dp_buf[i].begin();
    if(var_type==0) {//snp. gets DP at inteval that contains var_start/var_stop. else DP=0
      while(it1!=dp_buf[i].end() && it1->stop < var_start)  {
	it1++;
      }
      if(it1->start <= var_start && it1->stop >= var_start) {
	_out_dp[i] = it1->depth;	
	_out_gq[i] = it1->gq;
      }
      else {
	_out_dp[i] = 0;
	_out_gq[i] = bcf_int32_missing;
      }
    }
    else {//indel. this calculates the avg DP across (var_start,var_stop)
      int a,b;//range we are getting depth for.
      if(var_type==1||var_type==3) a=var_start;//insertion
      else a=var_start+1;//deletion.
      b=var_stop;
      int   var_len=b-a+1;
      float num=0;
      int min_gq = 10000;
      while(it1!=dp_buf[i].end() && it1->stop < a) 
	it1++;

      while(it1!=dp_buf[i].end() && it1->start <= b) {
	if(it1->depth!=bcf_int32_missing) {	//summing DP
	  if(it1->start <= a && it1->stop >= b) 
	    num += it1->depth * var_len;	
	  else if(it1->start<=a)
	    num += (it1->stop - a + 1)*it1->depth;
	  else if(it1->stop>=b)
	    num += (b - it1->start + 1)*it1->depth;
	}
	//GQ. simply finding the mingq across the region.
	if( (a<=it1->start && b>=it1->start) || (a<=it1->stop && b>=it1->stop) )
	  if(it1->gq==bcf_int32_missing||it1->gq < min_gq)
	    min_gq = it1->gq;
	it1++;
      }
      _out_dp[i] = (int)(round(num/(float)var_len));
      _out_gq[i] = min_gq;
      assert(_out_dp[i]==bcf_int32_missing||_out_gq[i]==bcf_int32_missing||_out_dp[i]>=0 && _out_gq[i] >=00);
    }
  }

  if(DEBUG>3) {
    for(int i=0;i<_nsample;i++) cerr << _dp[i]<<":"<<_gq[i]<<"\t";
    cerr <<endl;
  }

  return(0);
}

int aggReader::next() {
  int n_allele;
  if (bcf_sr_next_line (var_rdr)) {
    vr->clear();
    //get variants for each reader.
    for(int i=0;i<nreader;i++) {

      if( bcf_sr_has_line(var_rdr,i) ) {
	line[i]=bcf_sr_get_line(var_rdr,i);

	if(line[i]->rid != cur_chr || line_count==0) {
	  cur_chr = line[i]->rid;
	  while(dp_chr!=cur_chr)   	    moveDepthForward();
	  if(line_count>0)
	    cerr <<line_count<<" variants genotyped." <<endl;

	  if(DEBUG>0) {
	    cerr << "rid = "<<cur_chr<<" ("<<	  bcf_hdr_int2id(var_rdr->readers[i].header,BCF_DT_CTG,cur_chr)<<")"<<endl;
	    cerr << "dp_chr = "<<dp_chr<<"    dp_pos = "<< dp_pos<<endl;
	  }

	  for(int i=0;i<_nsample;i++)
	    dp_buf[i].clear();
	}

	n_allele = line[i]->n_allele;
	if(n_allele>2)      die("multiealleic site found");

	var_start=line[i]->pos;
	//get var_stop
	if(line[i]->n_allele!=2)
	  die("multiallelic site found");
	else {
	  int l1 = strlen(line[i]->d.allele[0]);
	  int l2 = strlen(line[i]->d.allele[1]);
	  if(DEBUG>1) cerr << "(l1,l2) = "<< l1 << " " << l2 << endl;
	  if(l1==1 && l2==1) {
	    var_stop=var_start;
	    var_type=0;
	  }
	  else if(l1>l2) {//deletion. calculate avg DP across it
	    var_stop = var_start + l1 - 1;
	    var_type=2;
	  }
	  else if(l1<l2) {//insertion. check for coverage at +1 position. is this sensible??
	    var_stop = var_start + l1;	    
	    var_type=1;
	  }
	  else {
	    var_stop = var_start + l2 - 1;//complex variant.
	    var_type=3;
	  }
	}
      }
    }


    setDepth();
    if(DEBUG>1) cout << var_start+1 <<"-"<<var_stop+1 << endl;
    int offset=0;

    //fills in FORMAT for output line.
    out_line->qual=0;
    for(int i=0;i<nreader;i++) {
      bcf_hdr_t *hdr = var_rdr->readers[i].header;
      int32_t *gt=&(vr->gt[2*offset]);
      int32_t *dp1=&(vr->dp[offset]);
      int32_t *dpf=&(vr->dpf[offset]);
      int32_t *pf=&(vr->pf[offset]);
      int32_t *gq1= &(vr->gq[offset]);
      int32_t *ad=&(vr->ad[n_allele*offset]);
      bool has_line=bcf_sr_has_line(var_rdr,i);
      int ntmp = bcf_hdr_nsamples(hdr);
      if(has_line) {
	int nval = 2*ntmp;
	int ngt_read = bcf_get_genotypes(hdr, line[i], &gt, &nval);
	if(ngt_read!=(2*ntmp) && ngt_read!=ntmp)
	  die("incorrect ploidy at "+to_string(static_cast<long long>(line[i]->pos+1))+" ngt_read="+to_string(static_cast<long long>(ngt_read)));

	if(ngt_read==ntmp) {//this is a hack to deal with sites where everyone has been made haploid.
	  int *work =  new int[ntmp];
	  //	  cerr << "WARNING: position "<< line[i]->pos+1 <<" had an all-haploid site, this can be caused by gvcf hemizygous rules"<<endl;
	  for(int j=0;j<ntmp;j++)  
	    work[j]=gt[j];
	  for(int j=0;j<ntmp;j++){
	    if(work[j]!=bcf_gt_missing) {
	      gt[2*j]=bcf_gt_unphased(0);
	      gt[2*j+1]=work[j];
	    }
	    else {
	      gt[2*j]=bcf_gt_missing;
	      gt[2*j+1]=bcf_gt_missing;
	    }
	  }
	  delete[] work;
	}		       

	bcf_get_format_int32(hdr, line[i], "AD", &ad, &nval);

	nval=ntmp;
	bcf_get_format_int32(hdr, line[i], "DPF", &dpf, &nval);
	bcf_get_format_int32(hdr, line[i], "DP", &dp1, &nval);
	bcf_get_format_int32(hdr, line[i], "GQ", &gq1, &nval);
	bcf_get_format_int32(hdr, line[i], "PF", &pf, &nval);

	//updates pos chrom etc
	bcf_update_id(out_hdr, out_line, ".");      
	bcf_update_alleles(out_hdr, out_line, (const char**)line[i]->d.allele,line[i]->n_allele);
	// int32_t tmpi = bcf_hdr_id2int(out_hdr, BCF_DT_ID, "PASS");
	// bcf_update_filter(out_hdr, out_line, &tmpi, 1);
	out_line->rid = line[i]->rid;
	out_line->pos = line[i]->pos;
	if(line[i]->qual > out_line->qual)
	  out_line->qual =  line[i]->qual ;
      }

      for(int j=0;j<ntmp;j++) {
	if(!has_line || (gt[j*2]==bcf_gt_missing && gt[j*2+1]==bcf_gt_missing) )  {//missing vcf entry. fill from dp_buf
	  dp1[j]=_out_dp[j+offset];
	  gq1[j]=_out_gq[j+offset];
	  if(gq1[j]>0||dp1[j]>0) {
	    gt[j*2]=bcf_gt_unphased(0);
	    gt[j*2+1]=bcf_gt_unphased(0);
	  }
	  else {
	    gt[j*2]=bcf_gt_missing;
	    gt[j*2+1]=bcf_gt_missing;	      
	  }
	}
	else if(var_type!=0) {//not a SNP? fill DP from AD
	  dp1[j] = ad[j*2]+ad[j*2+1];
	}

	//fixes the sporadic hemizgyotes
	if((bcf_gt_is_missing(gt[j*2])||gt[j*2]==bcf_int32_vector_end) && !bcf_gt_is_missing(gt[j*2+1]))
	  gt[j*2]=bcf_gt_unphased(0);
	if((bcf_gt_is_missing(gt[j*2+1])||gt[j*2+1]==bcf_int32_vector_end) && !bcf_gt_is_missing(gt[j*2]))
	  gt[j*2+1]=bcf_gt_unphased(0);

      }
      offset+=ntmp;
    }


    //updates output-line
    assert(    bcf_hdr_nsamples(out_hdr)==_nsample );
    bcf_update_genotypes(out_hdr,out_line,vr->gt,_nsample*2); 
    bcf_update_format_int32(out_hdr,out_line,"GQ",vr->gq,_nsample);
    bcf_update_format_int32(out_hdr,out_line,"DP",vr->dp,_nsample );
    bcf_update_format_int32(out_hdr,out_line,"DPF",vr->dpf,_nsample );
    bcf_update_format_int32(out_hdr,out_line,"AD",vr->ad,_nsample*n_allele );
    bcf_update_format_int32(out_hdr,out_line,"PF",vr->pf,_nsample );
    line_count++;
    return(1);
  }

  return(0);//finished  
}

//adds sum(DP),sum(AD),mean(PF),AC,AN to INFO.
void aggReader::annotate_line() {
  int32_t ac=0,an=0,sum_dp=0;
  int32_t sum_ad[2] = {0,0};
  float pf = 0.;
  float nalt=0;// number of genotypes containing an ALT allele.
  for(int i=0;i<_nsample;i++) {
    //allele counts
    if(vr->gt[i*2]!=bcf_gt_missing) {
      ac+=bcf_gt_allele(vr->gt[i*2]);
      an++;
    }
    if(vr->gt[i*2+1]!=bcf_gt_missing) {
      ac+=bcf_gt_allele(vr->gt[i*2+1]);
      an++;
    }

    //depth
    if(vr->dp[i]!=bcf_int32_missing) sum_dp+=vr->dp[i];

    //ad
    if(vr->gt[i*2+1]!=bcf_gt_missing && vr->gt[i*2]!=bcf_gt_missing) {
      //hom
      if(vr->ad[i*2]!=bcf_int32_missing)
	sum_ad[0]+=vr->ad[i*2];
      else
	sum_ad[0]+=vr->dp[i];
      //alt
      if(bcf_gt_allele(vr->gt[i*2])>0 || bcf_gt_allele(vr->gt[i*2+1])>0)
	sum_ad[1] += vr->ad[i*2+1];
    }
  
    //pf
    if(vr->gt[i*2+1]!=bcf_gt_missing && vr->gt[i*2]!=bcf_gt_missing)  {
      if(bcf_gt_allele(vr->gt[i*2])>0 || bcf_gt_allele(vr->gt[i*2+1])>0) {
	pf+=vr->pf[i];
	nalt++;
      }
    }
  }
  if(ac>an)  {
    for(int i=0;i<_nsample;i++)
      cerr << bcf_gt_allele(vr->gt[i*2])<<"/"<<bcf_gt_allele(vr->gt[i*2+1])<<"\t";
    cerr<<endl;
    cerr<<bcf_gt_missing<<endl;
    cerr<<bcf_int32_vector_end<<endl;
    cerr<< bcf_int16_vector_end<<endl;
    die("problem with genotyping at "+to_string(static_cast<long long>(out_line->pos+1)));
  }

  if(pf>0.0 && nalt>0)  pf/=nalt;
  bcf_update_info_int32(out_hdr, out_line, "AN", &an, 1);
  bcf_update_info_int32(out_hdr, out_line, "AC", &ac, 1);
  bcf_update_info_float(out_hdr, out_line, "PF", &pf, 1);
  bcf_update_info_int32(out_hdr, out_line, "AD", sum_ad, 2);
  bcf_update_info_int32(out_hdr, out_line, "DP", &sum_dp, 1);
  bcf_update_filter(out_hdr,out_line,NULL,0);//just wipe the filters.
}

int aggReader::writeVcf(const char *output_file,char *output_type,int n_threads ) {
  int nwritten=0;
  if(!output_type) output_type="v";
  string mode = "w" + (string) output_type;
  htsFile *out_fh  = hts_open(output_file ? output_file : "-", mode.c_str());

  if(!out_fh)
    die("problem opening output file");
  else
    cerr << "Writing output to "<<output_file<<endl;
    
  if(n_threads>0)  
    hts_set_threads(out_fh,n_threads);

  out_line = bcf_init1();
  out_hdr = bcf_hdr_init("w");
  int repeat_count=0;
  bool force_samples=true;
  for(int i=0; i<nreader; i++) {
    bcf_hdr_t *hr = var_rdr->readers[i].header;
    for(int j=0;j<bcf_hdr_nsamples(hr);j++) {
      string sample_name=hr->samples[j];
      if ( bcf_hdr_id2int(out_hdr, BCF_DT_SAMPLE, sample_name.c_str())!=-1 ) {
	if(force_samples) {
	  cerr << "Warning duplicate sample found.\t" << sample_name;
	  sample_name += ":R"+to_string(static_cast<long long>(repeat_count++));
	  cerr << " -> "<< sample_name<<endl;
	}
	else
	  die("duplicate sample names. use --force-samples if you want to merge anyway");
      }
      bcf_hdr_add_sample(out_hdr,sample_name.c_str());
    }
  }
  fillHeader(out_hdr);
  //  bcf_hdr_combine(out_hdr, var_rdr->readers[0].header);
  copyContigs(var_rdr->readers[0].header,out_hdr);
  bcf_hdr_write(out_fh, out_hdr);

  while(next()) {
    if(var_start>=interval_start) {
      annotate_line();
      bcf_write1(out_fh, out_hdr, out_line) ;
      bcf_clear1(out_line) ;
      nwritten++;
    }
  }
  hts_close(out_fh);
  //  bcf_hdr_destroy(hdr);
  return(nwritten);
}


static void usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   genotypes samples from agg chunks\n");
    fprintf(stderr, "Usage:   agg genotype chunk1 [chunk2 [...]]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "    -r, --regions <region>              region to genotype eg. chr1 or chr20:5000000-6000000\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o,   --output-file <file>          output file name [stdout]\n");
    fprintf(stderr, "    -O,   --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -@, --thread INT                    number of threads [0]\n");
    fprintf(stderr, "\n");
    // fprintf(stderr, "Subset options:\n");
    // fprintf(stderr, "    -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n");
    // fprintf(stderr, "    -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n");
    // fprintf(stderr, "\n");
    // fprintf(stderr, "Filter options:\n");
    // fprintf(stderr, "    TBA");
    // fprintf(stderr, "\n");
    exit(1);
}


int view1(int argc,char **argv) {
    int c;
    string region="";
    int n_threads=0;
    if(argc<3) usage();
    static struct option loptions[] =    {
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"regions",1,0,'r'},
        {"output-type",1,0,'O'},
        {"output-file",1,0,'o'},
	{"thread",1,0,'@'},
        {0,0,0,0}
    };
    char *output=NULL,*output_type="b";
    string sample_names;
    bool sample_is_file=0;
    while ((c = getopt_long(argc, argv, "r:@:o:O:s:S:",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 'O':
            switch (optarg[0]) {
            case 'b': output_type = optarg; break;
            case 'u': output_type = optarg; break;
            case 'z': output_type = optarg; break;
            case 'v': output_type = optarg; break;
            default: die("The output type \""+(string)optarg+"\" not recognised\n");
            };
            break;
	case '@': n_threads = atoi(optarg); break;
        case 'o': output = optarg; break;
        case 's': sample_names = optarg; break;
        case 'S': sample_names = optarg; sample_is_file = 1; break;
        case 'r': region = optarg; break;    
        case '?': usage();
        default: die("Unknown argument:"+(string)optarg+"\n");
        }
    }

    argc    = argc; argv = argv;  
    optind++;
    if(optind==argc)
      die("no input provided");

    vector<string> input;
    while(optind<argc) {
      cerr << argv[optind]<<" ";
      input.push_back(argv[optind++]);
    }
    cerr<<endl;
    if(output==NULL)
      die("no output specified");
    if(fileexists(output))
      die(((string)output)+" exists! Will not overwrite");

    if(region.find(",")<region.size())
        die("only handles a single region. "+region);

    cerr << "n_threads = "<<n_threads<<endl;
    aggReader ag(input,region);
    ag.writeVcf(output,output_type,n_threads);
    
    return(0);
}

