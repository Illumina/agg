#include "agg.h"

//iterates through a sites file and builds our dictionary.
int buildSiteList(const string& fname, map<string,map<unsigned int,map< pair<string,string> ,marker> > >  & d) {
    cerr << "Getting variant database from "<< fname << endl;
    bcf_srs_t *rdr = bcf_sr_init(); 

    if(!(bcf_sr_add_reader (rdr, fname.c_str()) ))
        die("problem opening "+fname);

    bcf1_t * line;
    int nvar = 0;
    while (bcf_sr_next_line (rdr)) {
        line = bcf_sr_get_line(rdr, 0);
        if(line->n_allele != 2) {
            cerr << endl << line->pos + 1 << endl;
            die("sites had a line with n_allele != 2");
        }
        string chrom = bcf_seqname(rdr->readers[0].header, line) ;
        unsigned int pos = line->pos;
        marker m(rdr->readers[0].header, line);
        pair<string, string> key(line->d.allele[0], line->d.allele[1]);
        d[chrom][pos].insert( pair<pair<string, string>, marker> (key, m));
        nvar++;
    }
    cerr << nvar << " variants currently in database..." << endl;
    return(nvar);
}

int updateSiteList(const string& file_list, map<string,map<unsigned int,map< pair<string,string> ,marker> > >  & d) {
    cerr << "Updating DB from gvcfs in "<< file_list <<"."<< endl;

    string sampleid, filename, blockfile;                               // Are these ever used?
    bcf1_t *line;
    ifstream is(file_list.c_str());
    if(!is) die("problem opening "+((string)file_list));
    bcf_hdr_t * header;
    int ngt,mgt_arr = 0, *gt_arr = NULL;
    pair<byte, byte> g;
    int nfile = 0;
    string hdr_copy;
    bool missing;
    int nvar = 0;
    while(is >> sampleid&&is >> filename&&is >> blockfile) {
        if(!is) continue; 
        if(nfile == 0) hdr_copy = filename;
        cerr << "Processing " << filename << endl;
        bcf_srs_t *reader= bcf_sr_init();

        // From here is the same as agg_collate
        if(!bcf_sr_add_reader (reader, filename.c_str())) 
            die("problem opening " + filename);

        header = reader->readers[0].header;

        while(bcf_sr_next_line (reader)) {
            line = bcf_sr_get_line(reader, 0);
            ngt = bcf_get_genotypes(header, line, &gt_arr, &mgt_arr);

            if(ngt == 2) {
                g.first = bcf_gt_allele(gt_arr[0]);         // ref
                g.second = bcf_gt_allele(gt_arr[1]);        // alt
                missing = (gt_arr[0] == bcf_gt_missing && gt_arr[1] == bcf_gt_missing);
            }
            else if(ngt==1) {
                g.first = bcf_gt_allele(gt_arr[0]);
                missing = gt_arr[0] == bcf_gt_missing;
            }
            else  die("Not a gvcf");

            if((g.first > 0 || g.second > 0) && !missing)
                nvar +=	addPosition(header, line, d, g);	
        }
        nfile++;
    }
    cerr << nfile << " files processed"<<endl;
    cerr << nvar << " novel variants found (not necessarily passing variants)" << endl;  
    return(nvar);
}


int updateSamples(const string& infile, const string& outfile) {
    cerr << "Updating sample database..." << endl;
    map< string, pair<string, string> > sampledb;
    ifstream is(outfile.c_str());

    string id, var, hr;
    while(is >> id && is >> var && is >> hr) 
        sampledb[id] = pair<string, string>(var, hr);    
    is.close();
    is.open(infile.c_str());
    int n = 0;
    while(is >> id && is >> var && is >> hr)  {
        if(sampledb.count(id)) {
            die("ERROR: "+ id +" was already in database!\nExiting without changes...");
        }
        else {
            sampledb[id] = pair<string, string>(var, hr); 
            n++;
        }
    }  
    cerr << n << " new samples found" << endl;
    return(0);
}

int update1(int argc,char **argv) {

    if(argc != 4) {
        cerr << "Usage: agg update <path to database dir> <file_list.txt>" << endl;
        return(1);
    }
    cerr << "agg update is currently disabled"<<endl;
    return(1);
    string db = argv[2];
    string file_list = argv[3];
    updateSamples(file_list, db + "/samples.txt");
    map<string, map<unsigned int, map< pair<string, string>, marker> > > d;
    buildSiteList(db + "/sites.bcf", d);
    updateSiteList(file_list, d);
    cerr << "Writing new variant database...";
    writeSiteList(db + "/sites.bcf", d);

    cerr << "Finished." << endl;
    return(0);
}

