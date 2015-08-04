#include "sampleMerger.h"
#define DEBUG 0

sampleMerger::sampleMerger(const string& database_location) {
  string filelist=database_location+"/samples.txt";
  dbdir=database_location;
  samples_initialised = false;
  ifstream is(filelist.c_str());
  if(!is) die("problem opening "+((string)filelist));
  string f1,f2,f3;
  nsample_in_db=0;
  while(is)  {
    is >> f1;
    is >> f2;
    is >> f3;
    if(!is) continue;
    sampleids.push_back(f1);
    genomes.push_back(sampleReader(database_location+"/variants/"+f1+".bcf",database_location+"/blocks/"+f1+".bed.gz"));
    nsample_in_db++;
  }
  cerr << "nsample_in_db = " << nsample_in_db<<endl;
}

sampleMerger::~sampleMerger() {
  delete[] new_gts;
  delete[] new_dp;
  delete[] new_ad;
  delete[] new_gq;
}

void sampleMerger::wipeRow() {
  for(int i=0;i<nsample;i++) {
    new_gts[i*2]=bcf_gt_missing;    new_gts[i*2+1]=bcf_gt_missing;
    new_ad[i*2] = bcf_int32_missing; new_ad[i*2+1] = bcf_int32_missing;
    new_dp[i] = bcf_int32_missing;
    bcf_float_set_missing( new_gq[i] ); 
  }
}

int sampleMerger::setSamples(char *samples,bool isfile) {
  samples_initialised=true;
  nsample=0;
  if(!samples) {
    nsample =    nsample_in_db;
    cerr << "Including all "<<nsample<<" samples in database."<<endl;
    included.assign(nsample_in_db,true);
  }
  else {
    included.assign(nsample_in_db,false);
    set<string> keep;
    if(isfile) {
      ifstream is(samples);
      string id;
      if(!is) die("problem opening "+((string)samples));
      while(is>>id)
	keep.insert(id);
    }
    else {
      vector<string> tmp;
      string splitme = samples;
      strsplit(splitme,',',tmp);
      for(int i=0;i<tmp.size();i++)
	keep.insert(tmp[i]);
    }
    for(int i=0;i<nsample_in_db;i++)
      if(keep.count(sampleids[i]))  {
	included[i]=true;
	nsample++;
      }
    cerr<<"Keeping "<<nsample<<" samples parsed from "<<samples<<endl;;
  }
  new_gts = new int32_t[nsample*2];
  new_dp = new int32_t[nsample];
  new_ad = new int32_t[nsample*2];
  new_gq = new float[nsample];
  if(nsample==0) die("no samples left to analyse!");
  return(nsample);
}

int sampleMerger::setSites(vector<marker> * site_list,const string & chromosome) {
  assert(samples_initialised);
  chrom = chromosome;
  sites = site_list;
  for(int i=0;i<nsample_in_db;i++)
    if(included[i])
      genomes[i].setSites(sites,chrom);
  marker_idx = -1;
  return(0);
}


int sampleMerger::countGenotypes(vector<int> & counts) {
  counts.assign(3,0);
  int32_t g0,g1;
  for(int i=0;i<(2*nsample);i+=2){
    g0 = new_gts[i];
    g1 = new_gts[i+1];

    if(g0!=bcf_gt_missing || g1!=bcf_gt_missing) {//at least one call.      
      if(g0==bcf_gt_missing) g0 = 0;
      else g0=bcf_gt_allele(g0);
      if(g1==bcf_gt_missing) g1=0;
      else g1=bcf_gt_allele(g1);
      counts[g0+g1]++;
    }
  }
  return(0);
}

int sampleMerger::next() {
  marker_idx++;
  if(marker_idx >= (*sites).size()) return(0);
  (*sites)[marker_idx].AC=0;
  (*sites)[marker_idx].AN=0;
  (*sites)[marker_idx].npass=0;
  if(DEBUG>0)  cerr <<  (*sites)[marker_idx].pos <<"\t"<<(*sites)[marker_idx].ref <<"\t"<<(*sites)[marker_idx].alt<<"\t";
  //  cerr <<  (*sites)[marker_idx].pos <<"\t"<<(*sites)[marker_idx].ref <<"\t"<<(*sites)[marker_idx].alt<<"\n";
  int count=0;//indexes new_gts
  wipeRow();
  //#pragma omp parallel for  
  for(int i=0;i<nsample_in_db;i++) {
    if(included[i]){
      int i1 = 2*count;int i2=2*count+1;
      //mfg::next(marker & m,int *gt,float  *gq,int32_t *dp,int32_t *ad)
      genomes[i].next((*sites)[marker_idx],new_gts+i1,new_gq+count,new_dp+count,new_ad+i1);
      new_gts[i1]=genomes[i].g.gt[0]; 
      new_gts[i2]=genomes[i].g.gt[1];
      new_gq[count]=genomes[i].g.gq[0];
      new_dp[count]=genomes[i].g.dp[0];
      new_ad[i1] = genomes[i].g.ad[0];
      new_ad[i2] = genomes[i].g.ad[1];
      if(bcf_gt_allele(new_gts[i1])>0 || bcf_gt_allele(new_gts[i2])>0)//only count ALT passes
	(*sites)[marker_idx].npass+=genomes[i].g.pass;
      count++;
      i1+=2;i2+=2;
    }
  }
    
  if(DEBUG>0)  cerr << endl;
  return(1);
}


int sampleMerger::writeVcf(char *out_filename,int output_type) {

  htsFile *out_fh  = vcf_wopen(out_filename,output_type);
  int nsample=0;
  bcf_hdr_t *out_hdr = bcf_hdr_init("w");
  for (int i=0; i < nsample_in_db; i++)    {
    if(included[i]) {
      bcf_hdr_add_sample(out_hdr,sampleids[i].c_str());
      nsample++;
    }
  }
    
  fillHeader(out_hdr,nsample,dbdir+"/sites.bcf");

  // bcf_hdr_sync(out_hdr);
  bcf_hdr_write(out_fh, out_hdr);
  bcf1_t *rec = bcf_init1() ;
  int nwrote=0;
  int *arr = NULL, marr = 0;//AN/AC working array

  vector<int> counts;
  while(next() && (*sites).size()>0) {
    if(marker_idx >= (*sites).size()) break;
    marker m = (*sites)[marker_idx];
    //    cerr << marker_idx<<"\t"<<m.pos <<"\t"<<m.ref<<"\t"<<m.alt <<endl;
    rec->rid = bcf_hdr_name2id(out_hdr, chrom.c_str());
    rec->pos = m.pos;
    bcf_update_id(out_hdr, rec, ".");
    string alleles = m.ref + "," + m.alt;
    bcf_update_alleles_str(out_hdr, rec, alleles.c_str());
    rec->qual = (*sites)[marker_idx].qual;
    bcf_update_genotypes(out_hdr,rec,(void *)new_gts,nsample*2);
    bcf_update_format_float(out_hdr,rec,"GQ",new_gq,nsample);
    bcf_update_format_int32(out_hdr,rec,"DP",new_dp,nsample);
    bcf_update_format_int32(out_hdr,rec,"AD",new_ad,nsample*2);

    hts_expand(int,rec->n_allele,marr,arr);
    bcf_calc_ac(out_hdr,rec,arr,BCF_UN_FMT);
    m.AN = arr[0] + arr[1];
    m.AC = arr[1];
    bcf_update_info_int32(out_hdr,rec,"AN",&m.AN,1);
    bcf_update_info_int32(out_hdr,rec,"AC",&m.AC,1);

    float af = (float)m.AC/(float)m.AN;
    bcf_update_info_float(out_hdr, rec, "AF", &af, 1);	
    countGenotypes(counts);
    float passrate=0.;
    if((counts[1]+counts[2])>0)
      passrate = (float)(m.npass) / (float)(counts[1]+counts[2]);

    assert(passrate<=1 && passrate>=0);

    bcf_update_info_float(out_hdr, rec, "PASSRATE", &passrate, 1);	

    float callrate=0.;
    for(int j=0;j<3;j++) callrate += (float)(counts[j]);
    callrate /= (float)nsample;	assert(callrate<=1 && callrate>=0);
    bcf_update_info_float(out_hdr, rec, "CALLRATE", &callrate, 1);

    bcf_write1(out_fh, out_hdr, rec) ;
    bcf_clear1(rec) ;
    nwrote++;
  }
  hts_close(out_fh);
  bcf_hdr_destroy(out_hdr);
  free(arr);
  bcf_destroy(rec);
  return(nwrote);
}

int sampleMerger::siteSummary(vector<marker> * site_list,const string& chromosome) {
  assert(samples_initialised);
  chrom = chromosome;
  sites = site_list;
  int gt[2],g0,g1;
  int32_t ad[2],dp;
  float gq;  
  bool pass;
  for(int i=0;i<nsample_in_db;i++) {
    //      cerr << "Sample " << i << endl;
    if(included[i]) {
      genomes[i].setSites(sites,chrom);//set site list
      for(marker_idx=0;marker_idx<(*sites).size();marker_idx++) {
	if(i==0) {//initialise marker statistics. (this should actually be a separate function)
	  (*sites)[marker_idx].AC=0;
	  (*sites)[marker_idx].AN=0;
	  (*sites)[marker_idx].npass=0;
	  (*sites)[marker_idx].counts.assign(3,0);
	}
	genomes[i].next((*sites)[marker_idx],gt,&gq,&dp,ad);
	g0 = gt[0];	  g1=gt[1];
	pass = genomes[i].g.pass;
	if(g0!=bcf_gt_missing || g1!=bcf_gt_missing) {//at least one call.

	  if(g0==bcf_gt_missing) g0 = 0;
	  else g0=bcf_gt_allele(g0);
	  if(g1==bcf_gt_missing) g1=0;
	  else g1=bcf_gt_allele(g1);

	  (*sites)[marker_idx].AN+=2;
	  if(g0>0||g1>0) (*sites)[marker_idx].npass+=pass;
	  assert( (g0+g1)<3 && (g0+g1)>=0);
	  (*sites)[marker_idx].counts[g0+g1]++;
	  (*sites)[marker_idx].AC+=(g0+g1);
	}
      }
      genomes[i].close();//close files etc
    }
  }

  if(DEBUG>0)    
    for(marker_idx=0;marker_idx<(*sites).size();marker_idx++)  {
      cerr << (*sites)[marker_idx].pos <<"\t"<< (*sites)[marker_idx].ref <<"\t"<< (*sites)[marker_idx].alt <<"\t"<< (*sites)[marker_idx].AN <<"\t"<< (*sites)[marker_idx].AC <<"\t"<< (*sites)[marker_idx].npass <<"\t"<< (*sites)[marker_idx].qual<<endl;
      cerr << (*sites)[marker_idx].counts[0]<<","<<(*sites)[marker_idx].counts[1]<<","<< (*sites)[marker_idx].counts[2]<<endl;
    }

  return(0);
}

