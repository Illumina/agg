#include "agg_ingest1.h"

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

int ingest1(const char *input,const char*output) {
  cerr << "Input: " << input << "\tOutput: "<<output<<endl;

  kstream_t *ks;
  kstring_t str = {0,0,0};    
  gzFile fp = gzopen(input, "r");
    
  if(fp==NULL) {
    fprintf(stderr,"problem opening %s\n",input);
    exit(1);
  }

  char *out_fname = (char *)malloc(strlen(output)+5);
  strcpy(out_fname,output);
  strcat(out_fname,".dpt");
  if(fileexists(out_fname)) {
    fprintf(stderr,"%s file already exists. will not overwrite\n",out_fname);
    exit(1);
  }
  printf("depth: %s\n",out_fname);
  gzFile depth_fp = gzopen(out_fname, "wb");
  strcpy(out_fname,output);
  strcat(out_fname,".bcf");
  if(fileexists(out_fname)) {
    fprintf(stderr,"%s file already exists. will not overwrite\n",out_fname);
    exit(1);
  }
  printf("variants: %s\n",out_fname);
  htsFile *variant_fp=hts_open(out_fname,"wb");



  bcf1_t *bcf_rec = bcf_init();
  ks = ks_init(fp);
  htsFile *hfp=hts_open(input, "r");
  bcf_hdr_t *hdr =  bcf_hdr_read(hfp);
  args_t *norm_args = init_vcfnorm(hdr);


  //this is a hack to fix broken gvcfs where AD is incorrectly defined in the header.
  bcf_hdr_remove(hdr,BCF_HL_FMT,"AD");
  assert(  bcf_hdr_append(hdr,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">") == 0);


  //this is a hack to fix broken gvcfs where GQ is incorrectly labelled as float (v4.3 spec says it should be integer)
  bcf_hdr_remove(hdr,BCF_HL_FMT,"GQ");
  assert(  bcf_hdr_append(hdr,"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">") == 0);

  //here we add FORMAT/PF. which is the pass filter flag for alts.
  assert(  bcf_hdr_append(hdr,"##FORMAT=<ID=PF,Number=A,Type=Integer,Description=\"variant was PASS filter in original sample gvcf\">") == 0);



  bcf_hdr_write(variant_fp, hdr);
  kstring_t work1 = {0,0,0};            
  int buf[5];
  ks_tokaux_t aux;
  while(    ks_getuntil(ks, '\n', &str, 0) >=0) {
    //    fprintf(stderr,"%s\n",str.s);
    if(str.s[0]!='#')  {
      char *ptr = kstrtok(str.s,"\t",&aux);//chrom
      ptr = kstrtok(NULL,NULL,&aux);//pos
      work1.l=0;
      kputsn(str.s,ptr-str.s-1, &work1);   
      buf[0] =  bcf_hdr_name2id(hdr, work1.s);
      assert(      buf[0]>=0);
      buf[1]=atoi(ptr)-1;
      for(int i=0;i<3;i++)  ptr = kstrtok(NULL,NULL,&aux);// gets us to ALT
      bool is_variant=ptr[0]!='.';

      for(int i=0;i<3;i++)  ptr = kstrtok(NULL,NULL,&aux);// gets us to INFO

      //find END if it is there
      char *end_ptr=strstr(ptr,"END=") ;
      if(end_ptr!=NULL) 
	buf[2]=atoi(end_ptr+4)-1;
      else
	buf[2]=buf[1];

      ptr  = kstrtok(NULL,NULL,&aux);//FORMAT
      //find index of DP (if present)
      //if not present, dont output anything (indels ignored)

      char *DP_ptr = find_format(ptr,"DP");
      if(DP_ptr!=NULL) {
	buf[3]=atoi(DP_ptr);
	char *GQX_ptr = find_format(ptr,"GQX");
	assert(GQX_ptr!=NULL);
	buf[4]=atoi(GQX_ptr);
	//	printf("%d\t%d\t%d\t%d\t%d\n",buf[0],buf[1],buf[2],buf[3],buf[4]);
	if(gzwrite(depth_fp,buf,5*sizeof(int))!=(5*sizeof(int)))
	  die("ERROR: problem writing "+(string)out_fname+".dpt");
      }
      if(is_variant) {//wass this a variant? if so write it out to the bcf
	vcf_parse(&str,hdr,bcf_rec);

	int32_t pass = bcf_has_filter(hdr, bcf_rec, ".");
	bcf_update_format_int32(hdr,bcf_rec,"PF",&pass,1);
	bcf_update_filter(hdr,bcf_rec,NULL,0);
	if(bcf_rec->n_allele>2) {//split multi-allelics (using vcfnorm.c from bcftools1.2)
	  split_multiallelic_to_biallelics(norm_args,bcf_rec );
	  for(int i=0;i<norm_args->ntmp_lines;i++){
	    remove_info(norm_args->tmp_lines[i]);
	    bcf_write1(variant_fp, hdr, norm_args->tmp_lines[i]);
	  }
	}
	else {
	  remove_info(bcf_rec);
	  bcf_write1(variant_fp, hdr, bcf_rec) ;
	}
	//	fprintf(stderr,"%s\n",str.s);
      }

    }
  }
  bcf_hdr_destroy(hdr);

  ks_destroy(ks);
  gzclose(fp);
  gzclose(depth_fp);  
  free(str.s);
  free(work1.s);
  hts_close(hfp);
  hts_close(variant_fp);

  fprintf(stderr,"Indexing %s\n",out_fname);
  bcf_index_build(out_fname, BCF_LIDX_SHIFT);
  return 0;
}

static void usage(){
  fprintf(stderr, "\n");
  fprintf(stderr, "About:   ingest a single sample gvcf into a variant-only bcf and depth track file (.dpt)\n");
  fprintf(stderr, "Usage:   agg ingest <input_gvcf>\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Required options:\n");
  fprintf(stderr, "    -o, --output <output_prefix>              agg will output output_prefix.bcf and output_prefix.dpt\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  exit(1);
}

int ingest_main(int argc,char **argv) {
  int c;
  char *output=NULL;
  if(argc<3) usage();

  static struct option loptions[] =    {
    {"output-file",1,0,'o'},
    {0,0,0,0}
  };


  while ((c = getopt_long(argc, argv, "o:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'o': output = optarg; break;
      default: die("Unknown argument:"+(string)optarg+"\n");
      }
  }
  if(!output)    die("the -o option is required");
  optind++;

  ingest1(argv[optind],output);
  cerr << "Done."<<endl;
  return(0);
}

