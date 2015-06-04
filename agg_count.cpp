#include "agg.h"
#include "sampleMerger.h"
#include "htslib/hts.h"

static void usage(){
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   summary statistics for samples from <database> at a region in the genome.\n");
    fprintf(stderr, "Usage:   agg count <database>\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Required options:\n");
    fprintf(stderr, "    -r, --regions <region>              region to genotype eg. chr1 or chr20:5000000-6000000\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Output options:\n");
    fprintf(stderr, "    -o,   --output-file <file>          output file name [stdout]\n");
    fprintf(stderr, "    -O,   --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Subset options:\n");
    fprintf(stderr, "    -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "    -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Filter options:\n");
    fprintf(stderr, "    TBA");
    fprintf(stderr, "\n");
    exit(1);
}

//args structure/parsing pilfered from bcftools view (vcfview.c)
typedef struct _args_t
{
    //    filter_t *filter;
    char *filter_str;
    int filter_logic;   // one of FLT_INCLUDE/FLT_EXCLUDE (-i or -e)

    bcf_srs_t *files;
    bcf_hdr_t *hdr, *hnull, *hsub; // original header, sites-only header, subset header
    char **argv, *format, *sample_names, *subset_fname, *targets_list, *regions_list;
    int argc, clevel, output_type, print_header, update_info, header_only, n_samples, *imap, calc_ac;
    int trim_alts, sites_only, known, novel, min_alleles, max_alleles, private_vars, uncalled, phased;
    int min_ac, min_ac_type, max_ac, max_ac_type, min_af_type, max_af_type, gt_type;
    int *ac, mac;
    float min_af, max_af;
    char *fn_ref, *fn_out, **samples;
    int sample_is_file, force_samples;
    char *include_types, *exclude_types;
    int include, exclude;
    htsFile *out;
}
args_t;

int count1(int argc,char **argv) {
    int c;
    args_t *args  = (args_t*) calloc(1,sizeof(args_t));
    if(argc<3) usage();
    static struct option loptions[] =    {
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {"regions",1,0,'r'},
        {"output-type",1,0,'O'},
        {"output-file",1,0,'o'},
        {0,0,0,0}
    };

    while ((c = getopt_long(argc, argv, "r:o:O:s:S:",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 'O':
            switch (optarg[0]) {
            case 'b': args->output_type = FT_BCF_GZ; break;
            case 'u': args->output_type = FT_BCF; break;
            case 'z': args->output_type = FT_VCF_GZ; break;
            case 'v': args->output_type = FT_VCF; break;
            default: die("The output type \""+(string)optarg+"\" not recognised\n");
            };
            break;
        case 'o': args->fn_out = optarg; break;
        case 's': args->sample_names = optarg; break;
        case 'S': args->sample_names = optarg; args->sample_is_file = 1; break;
        case 'r': args->regions_list = optarg; break;    
        case '?': usage();
        default: die("Unknown argument:"+(string)optarg+"\n");
        }
    }

    args->argc    = argc; args->argv = argv;  
    if( !args->regions_list ) 
        die("the --regions argument (-r) is required");
    optind++;
    string db = argv[optind];

    string region = args->regions_list;
    cerr << "Using database "<<db<<endl;
    if(args->fn_out)
      if(fileexists(args->fn_out))
	die(((string)args->fn_out)+" exists! Will not overwrite");

    if(args->fn_out)  cerr << "Writing output to "<<args->fn_out<<endl;
    else {
      cerr << "Writing to stdout" << endl;
      args->fn_out="-";
    }

    cerr << "Calculating summary information in region "<<region<<endl;
    if(region.find(",")<region.size())
        die("only handles a single region. "+region);
    string chromosome;
    if(region.find(":")<region.size())
        chromosome  = region.substr(0,region.find(":"));
    else
        chromosome = region;

    vector<marker> sites;
    sampleMerger  genomes(db);
    int nvar =  buildSiteList(db + "/sites.bcf",region,sites);
    if(nvar>0) {
      genomes.setSamples(args->sample_names,args->sample_is_file);
      cerr<<"Calculating summary statistics..."<<endl;
      genomes.siteSummary(&sites,chromosome);
    }
    else {
      cerr << "WARNING: no variants found in specified region. Writing header only vcf"<<endl;
    }
    cerr << "Writing output..."<<endl;
    writeSiteList(chromosome,sites,args->fn_out,args->output_type,genomes.nsample,db + "/sites.bcf");      
    return(0);
}

