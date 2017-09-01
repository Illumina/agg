#include "agg_ingest1.h"
//#define DEBUG

//just a struct to count some things
struct Counts
{
    int mnp,complex;
};

static void remove_hdr_lines(bcf_hdr_t *hdr, int type)
{
    int i = 0, nrm = 0;
    while ( i<hdr->nhrec )
    {
        if ( hdr->hrec[i]->type!=type ) 
	{ 
	    i++; 
	    continue; 
	}
        bcf_hrec_t *hrec = hdr->hrec[i];
        if ( type==BCF_HL_FMT )
        {
            // everything except FORMAT/GT
            int id = bcf_hrec_find_key(hrec, "ID");
            if ( id>=0 && !strcmp(hrec->vals[id],"GT") ) 
	    {
		i++; continue;
	    }
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

char *find_format(char *ptr,char *FORMAT) 
{
    char *fmt_ptr = strstr(ptr,FORMAT);
    if(fmt_ptr!=NULL) 
    {
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
class VarBuffer 
{
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
	_buf.push_back(v);
	int i = _buf.size()-1;
	while(i>0 && _buf[i]->pos < _buf[i-1]->pos) {
	    v=_buf[i-1];
	    _buf[i-1]=_buf[i];
	    _buf[i]=v;
	    i--;
	}
	return(1);
    }

    //write out variants to out file
    int flush(int pos,htsFile *outf,bcf_hdr_t *hdr_out) {
	int n = 0;
	while(_buf.size()>0 && (pos - _buf.front()->pos) > _w ) {
	    bcf1_t *rec = _buf.front();	    
	    //      cerr << _last_pos<<"<="<<rec->pos<<endl;
	    assert(_last_pos<=rec->pos);
	    if(   _last_pos!=rec->pos )
	    {
		_seen.clear();		
	    }	

	    string variant=(string)rec->d.allele[0] +"."+ (string)rec->d.allele[1];
	    if(_seen.count(variant)) {
		_ndup++;
	    }
	    else {
		_seen.insert(variant);
//		cerr << rec->rid<<":"<<rec->pos+1<<":"<<rec->d.allele[0]<<":"<<rec->d.allele[1]<<endl;
		bcf_write1(outf, hdr_out, rec);
	    }
	    _last_pos=rec->pos;
	    bcf_destroy1(rec);
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

//vt triple structure
struct Triple
{
    int pos_ref;
    int pos_alt;
    int len_ref;
    int len_alt;

    Triple() : pos_ref(0), pos_alt(0), len_ref(0), len_alt(0)
    {}

    Triple(int pos_ref, int pos_alt, int len_ref, int len_alt) : pos_ref(pos_ref), pos_alt(pos_alt), len_ref(len_ref), len_alt(len_alt)
    {}
};

//this is a (slightly) modfified version of vt's aggressive decompose 
//see https://github.com/atks/vt
int vt_aggressive_decompose(bcf1_t *v,bcf_hdr_t *hdr,vector<bcf1_t *> & buf)
{


    char** allele = bcf_get_allele(v);


    int new_no_variants=0;
    kstring_t new_alleles= {0,0,0};
    // Use alignment for decomposition of substitutions where REF
    // and ALT have different lengths and the variant is not an
    // insertion or deletion.
    
    // Perform alignment of REF[1:] and ALT[1:]
    NeedlemanWunsch nw(true);
    nw.align(allele[0] + 1, allele[1] + 1);
    nw.trace_path();
    // Force-align first characters
    if (allele[0][0] == allele[1][0])
	nw.trace.insert(nw.trace.begin(), NeedlemanWunsch::CIGAR_M);
    else
	nw.trace.insert(nw.trace.begin(), NeedlemanWunsch::CIGAR_X);
    nw.read--;
    nw.ref--;

    // Break apart alignment
    std::vector<Triple> chunks;
    bool hasError = false;
    int pos_ref = 0, pos_alt = 0, k = 0;
    Triple nextChunk(pos_ref, pos_alt, 0, 0);
    while (pos_ref <= nw.len_ref || pos_alt <= nw.len_read)
    {
	switch ((int32_t)nw.trace.at(k++))
	{
	case NeedlemanWunsch::CIGAR_M:
	    if (hasError)
		chunks.push_back(nextChunk);
	    nextChunk = Triple(pos_ref++, pos_alt++, 1, 1);
	    hasError = false;
	    break;
	case NeedlemanWunsch::CIGAR_X:
	    if (hasError)
		chunks.push_back(nextChunk);
	    nextChunk = Triple(pos_ref++, pos_alt++, 1, 1);
	    hasError = true;
	    break;
	case NeedlemanWunsch::CIGAR_D:
	    nextChunk.len_ref++;
	    pos_ref++;
	    hasError = true;
	    break;
	case NeedlemanWunsch::CIGAR_I:
	    nextChunk.len_alt++;
	    pos_alt++;
	    hasError = true;
	    break;
	}
    }
    if (hasError)
	chunks.push_back(nextChunk);


    int32_t pos1 = bcf_get_pos1(v);
    char* ref = strdup(v->d.allele[0]);
    char* alt = strdup(v->d.allele[1]);

    // old_alleles.l = 0;
    // bcf_variant2string(hdr, v, &old_alleles);

    for (size_t i=0; i<chunks.size(); ++i)
    {
	bcf1_t *nv = bcf_dup(v);
	bcf_unpack(nv, BCF_UN_ALL);
	bcf_set_pos1(nv, pos1+chunks[i].pos_ref);
	std::vector<int32_t> start_pos_of_phased_block;

                        
	new_alleles.l=0;
	for (int j=chunks[i].pos_ref; j<chunks[i].pos_ref+chunks[i].len_ref; ++j)
	    kputc(ref[j], &new_alleles);
	kputc(',', &new_alleles);
	for (int j=chunks[i].pos_alt; j<chunks[i].pos_alt+chunks[i].len_alt; ++j)
	    kputc(alt[j], &new_alleles);

	bcf_update_alleles_str(hdr, nv, new_alleles.s);
//	bcf_update_info_string(hdr, nv, "OLD_CLUMPED", old_alleles.s);
                    
	buf.push_back(nv);
	kputc('\0', &new_alleles);

	++new_no_variants;
    }
    if(new_alleles.l)
    {
	free(new_alleles.s);
    }

    free(ref);
    free(alt);
    
    return(new_no_variants);
}


//1. decompose MNPs/complex substitutions 
//2. normalises the output using bcftools norm algorithm
vector<bcf1_t *> atomise(bcf1_t *rec,bcf_hdr_t *hdr,Counts & counts)
{
    assert(rec->n_allele  == 2);
    char *ref=rec->d.allele[0];
    char *alt=rec->d.allele[1];
    int ref_len = strlen(ref);
    int alt_len = strlen(alt);
    vector<bcf1_t *> ret;
    if(ref_len>1 && ref_len==alt_len) //is MNP
    {
	char alleles[4] = "X,X";
	for(int i=0;i<ref_len;i++) 
	{
	    if(ref[i]!=alt[i]) 
	    {//new SNP
		bcf1_t *new_var = bcf_dup(rec);
		bcf_unpack(new_var, BCF_UN_ALL);
		alleles[0]=ref[i];
		alleles[2]=alt[i];
		new_var->pos+=i;
		bcf_update_alleles_str(hdr, new_var, alleles);	
		ret.push_back(new_var);
	    }
	    counts.mnp++;	    
	}
    }
    else if((ref_len!=alt_len) && (ref_len!=1) && (alt_len>1)) //complex substitution
    {
	vt_aggressive_decompose(rec,hdr,ret);
	counts.complex++;
    }
    else //variant already is atomic
    {
	bcf1_t *new_var = bcf_dup(rec);
	bcf_unpack(new_var, BCF_UN_ALL);
	ret.push_back(new_var);
    }  
    return(ret);
}

int ingest1(const char *input,const char *output,char *ref,bool exit_on_mismatch=true) 
{
    cerr << "Input: " << input << "\tOutput: "<<output<<endl;
    Counts counts;
    counts.mnp=0;
    counts.complex=0;
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
    //this is a hack to fix gvcfs where AD is incorrectly defined in the header. (vcf4.1 does not technically allow Number=R)
    bcf_hdr_remove(hdr_in,BCF_HL_FMT,"AD");
    assert(  bcf_hdr_append(hdr_in,"##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed. For indels this value only includes reads which confidently support each allele (posterior prob 0.999 or higher that read contains indicated allele vs all other intersecting indel alleles)\">") == 0);

    //enforce GQ as integer (required)
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
	    if(DP_ptr!=NULL) 
	    {
		buf[3]=atoi(DP_ptr);
		char *GQX_ptr = find_format(ptr,"GQ:");
		if(GQX_ptr==NULL)
		{
		    GQX_ptr = find_format(ptr,"GQX:");
		}
		assert(GQX_ptr!=NULL);
	
		//trying to reduce entropy on GQ to get better compression performance.
		buf[4]=atoi(GQX_ptr)/10;
		buf[4]*=10;

#ifdef DEBUG
		fprintf(stderr,"%d\t%d\t%d\t%d\t%d\n",buf[0],buf[1],buf[2],buf[3],buf[4]);
#endif 
		if(gzwrite(depth_fp,buf,5*sizeof(int))!=(5*sizeof(int)))
		{
		    die("ERROR: problem writing "+(string)out_fname+".tmp");		    
		}
	    }
	    if(is_variant)
	    {//wass this a variant? if so write it out to the bcf
		norm_args->ntotal++;
		vcf_parse(&str,hdr_in,bcf_rec);
		//	cerr<<bcf_rec->rid<<":"<<bcf_rec->pos<<endl;
		if(prev_rid!=bcf_rec->rid)
		{
		    vbuf.flush(variant_fp,hdr_out);
		}		 
		else
		{
		    vbuf.flush(bcf_rec->pos,variant_fp,hdr_out);
		}		    
		prev_rid=bcf_rec->rid;
		int32_t pass = bcf_has_filter(hdr_in, bcf_rec, ".");
		bcf_update_format_int32(hdr_out,bcf_rec,"PF",&pass,1);
		bcf_update_filter(hdr_out,bcf_rec,NULL,0);
		bcf_update_id(hdr_out,bcf_rec,NULL);
		bcf1_t **split_records=&bcf_rec;
		int num_split_records=1;
		if(bcf_rec->n_allele>2)
		{//split multi-allelics (using vcfnorm.c from bcftools1.3
		    norm_args->nsplit++;
		    split_multiallelic_to_biallelics(norm_args,bcf_rec);
		    split_records=norm_args->tmp_lines;
		    num_split_records=norm_args->ntmp_lines;
		}
		
		for(int i=0;i<num_split_records;i++)
		{
		    remove_info(split_records[i]);
		    vector<bcf1_t *> atomised_variants = atomise(split_records[i],hdr_out,counts);
		    for(size_t j=0;j<atomised_variants.size();j++)
		    {
			if(realign(norm_args,atomised_variants[j]) == ERR_REF_MISMATCH)
			{
			    if(exit_on_mismatch)
			    {
				die("vcf did not match the reference");
			    }
			    else
			    {
				norm_args->nskipped++;
			    }			    
			}
			vbuf.push_back(atomised_variants[j]);
//			bcf_destroy1(atomised_variants[j]);
		    }		    
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

    fprintf(stderr,"Variant lines   total/split/realigned/skipped:\t%d/%d/%d/%d\n", norm_args->ntotal,norm_args->nsplit,norm_args->nchanged,norm_args->nskipped);
    fprintf(stderr,"Decomposed %d MNPs\n", counts.mnp);
    fprintf(stderr,"Decomposed %d complex substitutions\n", counts.complex);
    
    destroy_data(norm_args);
    free(norm_args);


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


