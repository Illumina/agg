#include "agg_anno.h"

using namespace std;

inline float binomial_sd(int32_t a, int32_t b)
{
    if ((a + b) > 0 && a != bcf_int32_missing && b != bcf_int32_missing)
    {
	float n = a + b;
	return (((float)b / n - 0.5) / sqrt(.25 / n));
    }
    else
    {
	return (0.);
    }
}

//calculates the median and the mad
//note this modifies x
int madian(vector<float> &x, float &median, float &mad)
{
    sort(x.begin(), x.end());
    //  vector<float> r;
    //  r.resize(x.size());
    if (x.size() % 2 == 1)
	median = x[x.size() / 2];
    else
	median = ((float)x[x.size() / 2] + (float)x[x.size() / 2 + 1]) / 2;
    for (size_t i = 0; i < x.size(); i++)
	x[i] = fabs((float)(x[i] - median));
    sort(x.begin(), x.end());
    if (x.size() % 2 == 1)
	mad = x[x.size() / 2];
    else
	mad = ((float)x[x.size() / 2] + (float)x[x.size() / 2 + 1]) / 2;
    mad *= 1.4826;
    return (0);
}


int Residualiser::fit(vector<float> &x, vector<float> &y)
{
    if(x.size() <= 1)
    {
	die("no data to estimate parameters!");
    }

    //model fitting.
    vector<float> r;
    madian(r, _median, _mad);
//	    cerr << "n=" << y.size() << "\nmedian=" << _median << "\nmad=" << _mad << endl;
    return (0);
}

float Residualiser::residual(float x, float y)
{
    float ret = y - (_beta[0] + _beta[1] * x + _beta[2] * x * x);
    ret -= _median;
    ret /= _mad;
    return (ret);
}

BinResidualiser::BinResidualiser()
{
    _nsample = 0;
    _nbin = 0;
}

int BinResidualiser::initialise(int nbin, int nsample)
{
    _nsample = nsample;
    _bins.clear();
    _nbin = min(nsample, nbin);
    for (int i = 0; i < _nbin; i++)
	_bins.push_back(ceil(nsample * (float)i / (float)_nbin));
    _bins.push_back(nsample + 1);
	    //    for(int i=0;i<_nbin;i++)            cerr << i << " " << _bins[i]<<endl;
    return (0);
};

int BinResidualiser::fit(vector<float> &x, vector<float> &y)
{

    if(x.size() <= 1)
    {
	die("no data to estimate parameters!");
    }
    
    assert(_nbin > 0);
    vector<float> x1;
    x1.reserve(x.size() / _nbin);
    _median.clear();
    _mad.clear();
    _median.resize(_nbin);
    _mad.resize(_nbin);
    for (int b = 0; b < _nbin; b++)
    {
	for (size_t i = 0; i < x.size(); i++)
	    if (x[i] >= _bins[b] && x[i] < _bins[b + 1])
		x1.push_back(y[i]);
	madian(x1, _median[b], _mad[b]);
//		cerr << "(" << _bins[b] << "," << _bins[b + 1] << ") " << b << " n=" << x1.size() << "\nmedian=" << _median[b] << "\nmad=" << _mad[b] << endl;
	x1.clear();
    }
    return (0);
}

float BinResidualiser::residual(float x, float y)
{
    int b = 0;
    while (b < _nbin && _bins[b] <= x)
	b++;
    //    cerr << x << " "<<b<<endl;
    b--;
    if (!(b >= 0 && b < _nbin))
    {
	cerr << "x=" << x << endl;
	cerr << "b=" << b << " _nbin=" << _nbin << endl;
	exit(1);
    }
    return ((y - _median[b]) / _mad[b]);
}

Standardiser::Standardiser(const string &vcf1, const char *include)
{
    _sr = bcf_sr_init();
    _vcf = vcf1;

    if (!bcf_sr_add_reader(_sr, _vcf.c_str()))
    {
	die("Problem opening vcf1");
    }
    _hdr = _sr->readers[0].header;
    _nsample = bcf_hdr_nsamples(_hdr);
    vector<float> af_bins;
    for (int i = 0; i < 20; i++)
    {
	af_bins.push_back((float)i / 20.);
    }
    _filter = NULL;
    if (include != NULL)
    {
	cerr << "Using provided filter: "<< include<<endl;
        _filter = filter_init(_hdr, include);
    }
}

//first pass
//1. read in (filtered) depth and AF
//2. regress DP ~ AF (simple median per bin scheme)
//3. annodatate mean/sd of residuals.
int Standardiser::estimateParameters()
{
    cerr << "Estimating normalisation parameters..."<<endl;
    bcf1_t *line;
    int32_t info_dp, ac, an, dpf, dpa;
    int32_t *info_ab = (int32_t *)malloc(2 * sizeof(int32_t));
    int nval = 1;
    vector<float> depth;
    vector<float> alt_count;
    vector<float> ab;

    // binomial_data dpa_dpf;
    // dpa_dpf.n=0;
    // dpa_dpf.k.reserve(100e6);
    // dpa_dpf.size.reserve(100e6);
    double n_dpf = 0;
    double n_dpa = 0;
    depth.reserve(50e6);
    alt_count.reserve(50e6);
    int maxn = 0;
    int npass=0;
    while (bcf_sr_next_line(_sr))
    {
	assert(bcf_sr_has_line(_sr, 0));
	line = bcf_sr_get_line(_sr, 0);
	bcf_unpack(line, BCF_UN_INFO);
	if (_filter == NULL || filter_test(_filter, line, NULL))
	{
	    npass++;
	    int32_t *ptr = &info_dp;
	    if(bcf_get_info_int32(_hdr, line, "DP", &ptr, &nval) != 1)
	    {
		die("problem reading INFO/DP tag");
	    }
	    ptr = &ac;
	    if(bcf_get_info_int32(_hdr, line, "AC", &ptr, &nval) != 1)
	    {
		die("probleam reading INFO/AC");
	    }
	    ptr = &an;
	    if(bcf_get_info_int32(_hdr, line, "AN", &ptr, &nval) !=1)
	    {
		die("problem reading reading INFO/AN");
	    }
	    ptr = &dpf;
	    if(bcf_get_info_int32(_hdr, line, "DPF", &ptr, &nval) != 1)
	    {
		die("problem reading INFO/DPF");
	    }
	    ptr = &dpa;
	    if(bcf_get_info_int32(_hdr, line, "DPA", &ptr, &nval) != 1)
	    {
		die("problem reading INFO/DPA");
	    }
	    int n_ab = 2;
	    if(bcf_get_info_int32(_hdr, line, "AB", &info_ab, &n_ab) != 2)
	    {
		die("problem reading INFO/AB");
	    }
	    //if an==0 then it is not a valid site.
	    if (an > 0)
	    {
		alt_count.push_back((float)ac);
		depth.push_back((float)info_dp);
		ab.push_back(binomial_sd(info_ab[0], info_ab[1]));
		// dpa_dpf.k.push_back(dpa);
		// dpa_dpf.size.push_back(dpf+dpa);
		// dpa_dpf.n++;
		n_dpf += dpf;
		n_dpa += dpa;
	    }
	    if (an > maxn)
	    {
		maxn = an;
	    }
	}
    }
    cerr<<npass<<" sites passed provided filters"<<endl;
    _p_dpf = n_dpf / (n_dpf + n_dpa);
    cerr << "_p_dpf=" << _p_dpf << endl;
    //  fit_betab(&dpa_dpf,_alpha,_beta);
    bcf_sr_destroy(_sr);
    _dp_residuals.initialise(20, maxn);
    _dp_residuals.fit(alt_count, depth);
    _ab_residuals.initialise(20, maxn);
    _ab_residuals.fit(alt_count, ab);
    free(info_ab);
    return (0);
}

//second pass.
//re-open file. append INFO/S_DP - our normalised depth measure.
int Standardiser::standardise(char *output_file,char *output_type)
{
    cerr << "Annotating INFO fields..."<<endl;
    bcf1_t *line;
    int32_t info_dp, ac, an, dpf, dpa;
    int32_t *info_ab = (int32_t *)malloc(2 * sizeof(int32_t));
    int nval = 1;
    _sr = bcf_sr_init();
    if (!bcf_sr_add_reader(_sr, _vcf.c_str()))
	die("Problem opening vcf1");
    _hdr = _sr->readers[0].header;
    bcf_hdr_append(_hdr, "##INFO=<ID=S_AB,Number=1,Type=Float,Description=\"normalised alleleic balance measure (mean 0 sd 1)\">");
    bcf_hdr_append(_hdr, "##INFO=<ID=S_DP,Number=1,Type=Float,Description=\"normalised depth (mean 0 and sd 1)\">");
    //  bcf_hdr_append(_hdr, "##INFO=<ID=S_DPF,Number=1,Type=Float,Description=\"-log10(p-value) that DPA/(DPA+DPF) is extremely low. Assumes beta-binomial distribution \">");
    bcf_hdr_append(_hdr, "##INFO=<ID=S_DPF,Number=1,Type=Float,Description=\"proportion of bases filtered in samples with ALT genotypes\">");

    htsFile *out_fh = hts_open(output_file, output_type);
    bcf_hdr_write(out_fh, _hdr);
    while (bcf_sr_next_line(_sr))
    {
	assert(bcf_sr_has_line(_sr, 0));
	line = bcf_sr_get_line(_sr, 0);
	bcf_unpack(line, BCF_UN_INFO);
	int32_t *ptr = &info_dp;
	assert(bcf_get_info_int32(_hdr, line, "DP", &ptr, &nval) == 1);
	ptr = &ac;
	assert(bcf_get_info_int32(_hdr, line, "AC", &ptr, &nval) == 1);
	ptr = &an;
	assert(bcf_get_info_int32(_hdr, line, "AN", &ptr, &nval) == 1);
	ptr = &dpf;
	assert(bcf_get_info_int32(_hdr, line, "DPF", &ptr, &nval) == 1);
	ptr = &dpa;
	assert(bcf_get_info_int32(_hdr, line, "DPA", &ptr, &nval) == 1);
	int nab = 2;
	assert(bcf_get_info_int32(_hdr, line, "AB", &info_ab, &nab) == 2);
	//    cerr << line->pos+1<<" "<<info_ab[0]<<","<<info_ab[1]<<"="<<binomial_sd(info_ab[0],info_ab[1])<<endl;

	float s_dpf = 0;
	if (dpa > 0 && dpf > 0)
	{
	    //s_dpf = -log10(pbetab(dpa,dpa+dpf,_alpha,_beta,false,1000));
	    s_dpf = (float)dpf / (float)(dpa + dpf);
	    //      s_dpf = ( s_dpf - _p_dpf ) / sqrt ( _p_dpf*(1-_p_dpf)/(float)(dpa+dpf) );
	}

	//cerr << line->pos+1<<" "<<dpa<<"/"<<dpf<<" = "<<p_dpf<<endl;

	float x = (float)ac;
	if (an == 0)
	    x = 0.;
	float ndp = _dp_residuals.residual(x, (float)info_dp);
	float s_ab = _ab_residuals.residual(x, (float)binomial_sd(info_ab[0], info_ab[1]));
	//    cerr << line->pos+1<<" "<<info_ab[0]<<","<<info_ab[1]<<" = "<<s_ab<<endl;
	//    float call = (float)an/(2. * (float)_nsample);
	bcf_update_info_float(_hdr, line, "S_DP", &ndp, 1);
	bcf_update_info_float(_hdr, line, "S_AB", &s_ab, 1);
	//    bcf_update_info_float(_hdr, line, "CALL", &call, 1);
	bcf_update_info_float(_hdr, line, "S_DPF", &s_dpf, 1);
	//    bcf_update_info_float(_hdr, line, "CALLRATE", &callrate, 1);
	bcf_write1(out_fh, _hdr, line);
	bcf_clear1(line);
    }
    hts_close(out_fh);
    bcf_sr_destroy(_sr);
    free(info_ab);
    return (0);
}

static void usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   agg anno input.bcf\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -i, --include      filters to apply eg. -i 'QUAL>=10 && DP<100000 && HWE<10' \n");
    fprintf(stderr, "    -R, --regions      a set of variants that are trusted (eg. 1000G) \n");
    fprintf(stderr, "    -o, --output-file <file>          output file name [stdout]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed");
    fprintf(stderr, "\n");
    fprintf(stderr, "\n");
    exit(1);
}

int annotate1(int argc, char **argv)
{
    int c;
    char *regions = NULL;
    char *output_file="-";
    char output_type[3]="wv";  
    static struct option loptions[] = {
	{"include", required_argument, NULL, 'i'},
	{"regions", required_argument, NULL, 'R'},
	{"output-file", required_argument, NULL, 'o'},
	{"output-type", required_argument, NULL, 'O'},	
	{0, 0, 0, 0}};
    char *vcf1 = NULL, *include = NULL;
    if (argc < 2)
    {
	usage();
    }
    while ((c = getopt_long(argc, argv, "i:R:O:o:", loptions, NULL)) >= 0)
    {
	switch (c)
	{
	case 'i':
	    include = optarg;
	    break;
	case 'O':
	    output_type[1] =  optarg[0];
	    break;
	case 'o':
	    output_file = optarg;
	    break;	    
	case 'R':
	    regions = optarg;
	    break;
	default:
	    if (optarg != NULL)
	    {   
		die("Unknown argument:" + (string)optarg + "\n");
	    }
	    else
	    {
		die("unrecognised argument");
	    }
	}
    }
    optind++;
    if (optind == argc)
    {
	usage();
    }
    vcf1 = argv[optind];
    if (!vcf1)
    {
	die("missing input");
    }
    Standardiser s(vcf1, include);
    s.estimateParameters();
    s.standardise(output_file,output_type);
    return (0);
}
