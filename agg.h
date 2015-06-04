#pragma once

#include <math.h>
#include <map>
#include <getopt.h>
#include <vector>

#include "utils.h"


extern "C" {
#include "htslib/synced_bcf_reader.h"
}

using namespace std;

class marker {
public:
    marker(bcf_hdr_t *header, bcf1_t *line);
    ~marker();
    string ref, alt;        // Reference and alternate alleles
    uint32_t AN;            // number of alleles present (homref included)
    uint32_t AC;            // Allele count
    uint32_t npass;         // Number of alleles passed qc
    float passrate;         // proportion of ALT genotypes passing filters
    uint32_t pos;           // Position
    float qual;             // Quality
    vector<unsigned int> counts;        // counts of 0/0,0/1,1/1
};

class region {
public:
    region(const string& c,uint32_t a,uint32_t b);
    string chrom;
    uint32_t start, stop;
    string toString();
};


class homBlock {
public:
    homBlock(uint32_t a, uint32_t b, bool p, bool m);
    uint32_t start, stop;
    bool pass, missing;
};


int addPosition(bcf_hdr_t *header ,bcf1_t *line, map < string , map < unsigned int, map < pair< string, string > , marker > > >  &d, pair< byte, byte > g);

vector<region> getRegions(char *fname,int buf);

int isBlock(bcf_hdr_t *header,bcf1_t *line);

int collate1(int argc,char **argv);

int genotype1(int argc,char **argv);

int count1(int argc,char **argv);

int update1(int argc,char **argv);

double hwe(vector<unsigned int> count);

int buildSiteList(const string& fname, const string& region,  vector<marker> & out);

int writeSiteList(const string& fname, map<string,map<unsigned int,map< pair<string,string> ,marker> > >  & d);

int writeSiteList(const string& chromosome, 
		  vector<marker> & sites,
		  const string& fname, 
		  int output_type,
		  int nsample,
		  const string& sites_file);

htsFile *vcf_wopen(const string& out_filename, int output_type);

int fillHeader(bcf_hdr_t *hdr,int nsample,string sites_file);//fills in the standard stuff for an agg header.
int copyContigs(bcf_hdr_t *src,bcf_hdr_t *dst);
