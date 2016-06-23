#include "agg_ingest1.h"
//#define DEBUG

static void remove_hdr_lines(bcf_hdr_t *hdr, int type)
{
    int i = 0, nrm = 0;
    while ( i<hdr->nhrec )
    {
        if ( hdr->hrec[i]->type!=type ) { i++; continue; }
        bcf_hrec_t *hrec = hdr->hrec[i];
        if ( type==BCF_HL_FMT )
        {
            // everything except FORMAT/GT
            int id = bcf_hrec_find_key(hrec, "ID");
            if ( id>=0 && !strcmp(hrec->vals[id],"GT") ) { i++; continue; }
        }
        nrm++;
        hdr->nhrec--;
        if ( i < hdr->nhrec )
            memmove(&hdr->hrec[i],&hdr->hrec[i+1],(hdr->nhrec-i)*sizeof(bcf_hrec_t*));
        bcf_hrec_destroy(hrec);
    }
    if ( nrm ) bcf_hdr_sync(hdr);
}

void remove_info(bcf1_t *line)
{
    // remove all INFO fields
    if ( !(line->unpacked & BCF_UN_INFO) ) bcf_unpack(line, BCF_UN_INFO);

    int i;
    for (i=0; i<line->n_info; i++)
    {
        bcf_info_t *inf = &line->d.info[i];
        if ( inf->vptr_free )
        {
            free(inf->vptr - inf->vptr_off);
            inf->vptr_free = 0;
        }
        line->d.shared_dirty |= BCF1_DIRTY_INF;
        inf->vptr = NULL;
    }
    line->n_info=0;
}

char *find_format(char *ptr,char *FORMAT) {
  char *fmt_ptr = strstr(ptr,FORMAT);
  if(fmt_ptr!=NULL) {
    int idx=0;
    for(int i=0;i<(fmt_ptr-ptr);i++)
      if(ptr[i]==':')
	idx++;	
    ks_tokaux_t aux;
    ptr = kstrtok(ptr,"\t",&aux);
    ptr = kstrtok(NULL,NULL,&aux);//sample column
    ptr = kstrtok(ptr,":",&aux);
    for(int i=0;i<idx;i++) 
      ptr = kstrtok(NULL,NULL,&aux);

    return(ptr);
  }
  else
    return(NULL);
}

//small class to buffer bcf1_t records and sort them as they are inserted.
class VarBuffer {
public:
  VarBuffer(int w) {
    _w = w;//window size (physical bp)
    _last_pos=-1;
    _ndup=0;
  }
  ~VarBuffer() {
    cerr << "Dropped "<<_ndup<<" duplicated variants after normalization."<<endl;
  }
  //add a new variant (and sort if necessary)
  int push_back(bcf1_t *v) {
    bcf_unpack(v, BCF_UN_ALL);
    bcf1_t *tmp=bcf_dup(v);
    bcf_unpack(tmp, BCF_UN_ALL);
    _buf.push_back(tmp);
    int i = _buf.size()-1;
    while(i>0 && _buf[i]->pos < _buf[i-1]->pos) {
      tmp=_buf[i-1];
      _buf[i-1]=_buf[i];
      _buf[i]=tmp;
      i--;
    }
    return(1);
  }

  //write out variants to out file
  int flush(int pos,htsFile *outf,bcf_hdr_t *hdr_out) {
    int n = 0;
    while(_buf.size()>0 && (pos - _buf.front()->pos) > _w ) {
      //      cerr << _last_pos<<"<="<<_buf.front()->pos<<endl;
      assert(_last_pos<=_buf.front()->pos);
      if(   _last_pos!=_buf.front()->pos )  
	_seen.clear();

      string variant=(string)_buf.front()->d.allele[0] +"."+ (string)_buf.front()->d.allele[1];

      if(_seen.count(variant)) {
	_ndup++;
      }
      else {
	_seen.insert(variant);
	bcf_write1(outf, hdr_out, _buf.front());
      }
      _last_pos=_buf.front()->pos;
      bcf_destroy1(      _buf.front() );
      _buf.pop_front();
      n++;
    }    
    return(n);
  }  

  //write out variants to out file
  int flush(htsFile *outf,bcf_hdr_t *hdr_out) {
    if(_buf.size()>0) {
      int pos=_buf.back()->pos+_w+1;
      int ret = flush(pos,outf,hdr_out) ;
      _last_pos= -1;
      return(ret);
    }
    else{
      return(0);
    }
  }
  
private:
  int _w,_last_pos,_ndup;
  deque<bcf1_t *> _buf;  
  set <string> _seen; //list of seen variants at this position.
};



//decomposes MNPs into multiple records and pushes them into the buffer.
int decompose(bcf1_t *rec,bcf_hdr_t *hdr,VarBuffer & buf) {
  assert(rec->n_allele  == 2);
  char *ref=rec->d.allele[0];
  char *alt=rec->d.allele[1];
  int refl = strlen(ref);
  int altl = strlen(alt);
  int n=0;
  if(refl>1 && refl==altl) {//is MNP
    char alleles[4] = "X,X";
    for(int i=0;i<refl;i++) {
      if(ref[i]!=alt[i]) {//new SNP
	bcf1_t *new_var = bcf_dup(rec);
	bcf_unpack(new_var, BCF_UN_ALL);
	alleles[0]=ref[i];
	alleles[2]=alt[i];
	new_var->pos+=i;
	bcf_update_alleles_str(hdr, new_var, alleles);	
	buf.push_back(new_var);
	n++;
      }
    }
  }
  else {
    buf.push_back(rec);    
  }  
  return(n);
}

int ingest1(const char *input,const char *output,char *ref,bool exit_on_mismatch=true) {
  cerr << "Input: " << input << "\tOutput: "<<output<<endl;

  kstream_t *ks;
  kstring_t str = {0,0,0};    
  gzFile fp = gzopen(input, "r");
  VarBuffer vbuf(1000);
  int prev_rid = -1;
  if(fp==NULL) {
    fprintf(stderr,"problem opening %s\n",input);
    exit(1);
  }

  char *out_fname = (char *)malloc(strlen(output)+5);
  strcpy(out_fname,output);
  strcat(out_fname,".tmp");
  if(fileexists(out_fname)) {
    fprintf(stderr,"%s file already exists. will not overwrite\n",out_fname);
    exit(1);
  }
  printf("depth: %s\n",out_fname);
  gzFile depth_fp = gzopen(out_fname, "wb1");
  strcpy(out_fname,output);
  strcat(out_fname,".bcf");
  if(fileexists(out_fname)) {
    fprintf(stderr,"%s file already exists. will not overwrite\n",out_fname);
    exit(1);
  }
  printf("variants: %s\n",out_fname);
  htsFile *variant_fp=hts_open(out_fname,"wb1");
  if(variant_fp==NULL) {
    fprintf(stderr,"problem opening %s\n",input);
    exit(1);    
  }

  ks = ks_init(fp);
  htsFile *hfp=hts_open(input, "r");
  bcf_hdr_t *hdr_in =  bcf_hdr_read(hfp);
  if(bcf_hdr_id2int(hdr_in, BCF_DT_ID, "GQX")==-1)
      die("FORMAT/GQX is not present");
  if(bcf_hdr_id2int(hdr_in, BCF_DT_ID, "DP")==-1)
      die("FORMAT/DP is not present");
  if(bcf_hdr_id2int(hdr_in, BCF_DT_ID,"BLOCKAVG_min30p3a")==-1)
     die("INFO/BLOCKAVG_min30p3a is not present");

  hts_close(hfp);
  //this is a hack to fix gvcfs where AD is incorrectly defined in the header. (vcf4.2 does not technically allow Number=R)
  bcf_hdr_remove(hdr_in,BCF_HL_FMT,"AD");
  assert(  bcf_hdr_append(hdr_in,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">") == 0);

  //this is a hack to fix broken gvcfs where GQ is incorrectly labelled as float (v4.3 spec says it should be integer)
  bcf_hdr_remove(hdr_in,BCF_HL_FMT,"GQ");
  assert(  bcf_hdr_append(hdr_in,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">") == 0);


  //  bcf_hdr_t  *hdr_out=hdr_in;
  bcf_hdr_t *hdr_out =  bcf_hdr_dup(hdr_in);
  remove_hdr_lines(hdr_out,BCF_HL_INFO);
  remove_hdr_lines(hdr_out,BCF_HL_FLT);
  bcf_hdr_sync(hdr_out);

  //here we add FORMAT/PF. which is the pass filter flag for alts.
  assert(  bcf_hdr_append(hdr_out,"##FORMAT=<ID=PF,Number=A,Type=Integer,Description=\"variant was PASS filter in original sample gvcf\">") == 0);

  args_t *norm_args = init_vcfnorm(hdr_out,ref);
  norm_args->check_ref |= CHECK_REF_WARN;
  bcf1_t *bcf_rec = bcf_init();
  bcf_hdr_write(variant_fp, hdr_out);
  kstring_t work1 = {0,0,0};            
  int buf[5];
  ks_tokaux_t aux;
  int ndec=0;
  int ref_len,alt_len;
  while(    ks_getuntil(ks, '\n', &str, 0) >=0) {
    //    fprintf(stderr,"%s\n",str.s);
    if(str.s[0]!='#')  {
      char *ptr = kstrtok(str.s,"\t",&aux);//chrom
      ptr = kstrtok(NULL,NULL,&aux);//pos
      work1.l=0;
      kputsn(str.s,ptr-str.s-1, &work1);   
      buf[0] =  bcf_hdr_name2id(hdr_in, work1.s);
      assert(      buf[0]>=0);
      buf[1]=atoi(ptr)-1;
      ptr = kstrtok(NULL,NULL,&aux);//ID
      ptr = kstrtok(NULL,NULL,&aux);//REF

      ref_len=0;
      while(ptr[ref_len]!='\t') ref_len++;

      ptr = kstrtok(NULL,NULL,&aux);//ALT

      bool is_variant=false;
      alt_len=0;
      while(ptr[alt_len]!='\t') alt_len++;
      if(ptr[0]!='.') 
	is_variant=true;
      

      for(int i=0;i<3;i++)  ptr = kstrtok(NULL,NULL,&aux);// gets us to INFO

      //find END if it is there
      char *end_ptr=strstr(ptr,"END=") ;
      if(end_ptr!=NULL) 
	buf[2]=atoi(end_ptr+4)-1;
      else
	buf[2]=buf[1]+alt_len-1;

      ptr  = kstrtok(NULL,NULL,&aux);//FORMAT
      //find index of DP (if present)
      //if not present, dont output anything (indels ignored)
      
      char *DP_ptr = find_format(ptr,"DP:");
      if(DP_ptr!=NULL) {
	buf[3]=atoi(DP_ptr);
	char *GQX_ptr = find_format(ptr,"GQX:");
	assert(GQX_ptr!=NULL);
	
	//trying to reduce entropy on GQ to get better compression performance.
	//1. rounds down to nearest 10. 
	//2. sets gq to min(gq,100). 
	buf[4]=atoi(GQX_ptr)/10;
	buf[4]*=10;
	if(buf[4]>100) buf[4]=100;
#ifdef DEBUG
	fprintf(stderr,"%d\t%d\t%d\t%d\t%d\n",buf[0],buf[1],buf[2],buf[3],buf[4]);
#endif 
	if(gzwrite(depth_fp,buf,5*sizeof(int))!=(5*sizeof(int)))
	  die("ERROR: problem writing "+(string)out_fname+".tmp");
      }
      if(is_variant) {//wass this a variant? if so write it out to the bcf
	norm_args->ntotal++;
	vcf_parse(&str,hdr_in,bcf_rec);
	//	cerr<<bcf_rec->rid<<":"<<bcf_rec->pos<<endl;
	if(prev_rid!=bcf_rec->rid) 
	  vbuf.flush(variant_fp,hdr_out);
	else
	  vbuf.flush(bcf_rec->pos,variant_fp,hdr_out);
	prev_rid=bcf_rec->rid;
	int32_t pass = bcf_has_filter(hdr_in, bcf_rec, ".");
	bcf_update_format_int32(hdr_out,bcf_rec,"PF",&pass,1);
	bcf_update_filter(hdr_out,bcf_rec,NULL,0);
	if(bcf_rec->n_allele>2) {//split multi-allelics (using vcfnorm.c from bcftools1.3
	  norm_args->nsplit++;
	  split_multiallelic_to_biallelics(norm_args,bcf_rec );
	  for(int i=0;i<norm_args->ntmp_lines;i++){
	    remove_info(norm_args->tmp_lines[i]);
	    if(realign(norm_args,norm_args->tmp_lines[i]) != ERR_REF_MISMATCH)
	      ndec+=decompose(norm_args->tmp_lines[i],hdr_out,vbuf);
	    else
	      if(exit_on_mismatch)
		die("vcf did not match the reference");
	      else
		norm_args->nskipped++;
	  }
	}
	else {
	  remove_info(bcf_rec);
	  if( realign(norm_args,bcf_rec) !=  ERR_REF_MISMATCH)
	    ndec+=decompose(bcf_rec,hdr_out,vbuf);
	  else
	    if(exit_on_mismatch)
	      die("vcf did not match the reference");
	    else
	      norm_args->nskipped++;
	}
	vbuf.flush(bcf_rec->pos,variant_fp,hdr_out);
      }
    }
  }
  vbuf.flush(variant_fp,hdr_out);
  bcf_hdr_destroy(hdr_in);
  bcf_hdr_destroy(hdr_out);
  bcf_destroy1(bcf_rec);
  ks_destroy(ks);
  gzclose(fp);
  gzclose(depth_fp);  
  free(str.s);
  free(work1.s);
  hts_close(variant_fp);
  destroy_data(norm_args);
  fprintf(stderr,"Variant lines   total/split/realigned/skipped:\t%d/%d/%d/%d\n", norm_args->ntotal,norm_args->nsplit,norm_args->nchanged,norm_args->nskipped);
  fprintf(stderr,"Decomposed %d MNPs\n", ndec);


  fprintf(stderr,"Indexing %s\n",out_fname);
  bcf_index_build(out_fname, BCF_LIDX_SHIFT);
  free(out_fname);
  return 0;
}

static void usage(){
  fprintf(stderr, "\n");
  fprintf(stderr, "About:   ingests a single sample gvcf into a variant-only .bcf and tempory depth interval (.tmp)\n");
  fprintf(stderr, "Usage:   agg ingest1 <input_gvcf>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Required options:\n");
  fprintf(stderr, "    -o, --output <output_prefix>      agg will output output_prefix.bcf and output_prefix.tmp\n");
  fprintf(stderr, "    -f, --fasta-ref <file>            reference sequence\n");
  fprintf(stderr, "        --ignore-non-matching-ref     skip non-matching ref alleles (will warn)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(1);
}

int ingest_main(int argc,char **argv) {
  int c;
  char *output=NULL;
  char *ref=NULL;
  if(argc<3) usage();

  static struct option loptions[] =    {
    {"output-file",1,0,'o'},
    {"fasta-ref",1,0,'f'},
    {"ignore-non-matching-ref",0,0,1},
    {0,0,0,0}
  };

  bool exit_on_mismatch=true;
  while ((c = getopt_long(argc, argv, "o:f:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'o': output = optarg; break;
      case 'f': ref = optarg; break;
      case 1: exit_on_mismatch=false;break;
      default: 
	if(optarg!=NULL) die("Unknown argument:"+(string)optarg+"\n");
	else die("unrecognised argument");
      }
  }
  if(!output)    die("the -o option is required");
  if(!ref)    die("the -f option is required");
  optind++;
  if(optind==argc)
    die("no input provided");
  ingest1(argv[optind],output,ref,exit_on_mismatch);
  cerr << "Done."<<endl;
  return(0);
}


