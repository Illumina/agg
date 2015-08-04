#include "agg.h"

int collate1(int argc, char **argv) {

    // Check the correct arguments were supplied
    if(argc!=5) {
        cerr << "Usage: agg collate database_dir/ output.bcf region" << endl;
        return(1);
    }

    string dbdir = argv[2];
    string filelist = (dbdir+"/samples.txt").c_str();       // list of gvcfs
    char *outf = argv[3];           // output vcf file

    // Check the output file doesn't already exist
    if(fileexists(outf))
        die(((string)outf)+" exists! Will not overwrite");

    cerr << "Reading gvcfs from " << filelist << "." << endl;
    
    char *reg;          // region specified in arguments
    bool useRegion;     // specifies if we are looking at a specific region or not
    if(argc == 5) {
        reg = argv[4];      // Region of interest
        cerr << "Region: "<< reg << endl;
        useRegion = argc == 5;
    }    
    else {
        cerr << "No region set (this might be slow and use lots of memory)" << endl;
    }
    cerr << "Writing output to "<<outf<<" (bgzipped vcf)"<<endl;

    // Open file list
    ifstream is(filelist.c_str());
    if(!is) die("problem opening "+((string)filelist));

    string filename;
    bcf1_t * line;
    map< string, map< unsigned int, map< pair<string, string>, marker> > > d;       // map of variants
    int ngt, mgt_arr = 0, * gt_arr = NULL;                                          // number of genotypes? some array, gentoype array?
    pair<byte, byte> g;                                                             // genotype
    bcf_hdr_t * header;                                                             // header
    int nfile = 0;                                                                  // number of vcf files
    string hdr_copy;                                                                // copy of the header
    bool missing;                                                                   // are their missing genotypes
    bcf_hdr_t * ohdr;                                                               // our header
    int nvar = 0;                                                                   // number of variants
    int missint = 0;                                                                // number of missing genotypes
    string f2,f3;
    while(is) {
        is >> filename;
	is >> f2;is>>f3;
	filename = dbdir+"/variants/"+filename+".bcf";
        if(!is) continue; 
        if(nfile == 0) hdr_copy = filename;
        cerr << "Processing " << filename << endl;

        // Initialise a bcf stream reader
        bcf_srs_t *reader = bcf_sr_init();

        // Exit gracefully if we can't set the region
        if(useRegion) {
            reader->require_index = 1;
            if(bcf_sr_set_regions(reader, reg, 0) != 0)
                die("Problem setting regions");
        }

        // Exit gracefully if we can't open the gVcf with the stream reader
        if(!bcf_sr_add_reader (reader, filename.c_str())) 
            die("problem opening " + filename);

        header = reader->readers[0].header;  

        // Read the gVcf
        while(bcf_sr_next_line(reader)) {
            line = bcf_sr_get_line(reader, 0);
            ngt = bcf_get_genotypes(header, line, &gt_arr, &mgt_arr);
            // Get the genotype and count missing ones
            if(ngt == 2) {
                g.first = bcf_gt_allele(gt_arr[0]);
                g.second = bcf_gt_allele(gt_arr[1]);  
                missing = (gt_arr[0] == bcf_gt_missing && gt_arr[1] == bcf_gt_missing);
            }
            else if(ngt == 1) {
                g.first = bcf_gt_allele(gt_arr[0]);
                missing = gt_arr[0] == bcf_gt_missing;
            }
            else die("Not a gvcf");

            // If its a variant add it to our database
	    // note this allows variants where guys are homref. this can happen when there is evidence for a variant but 
	    // not enough to genotype it as het 
            if((g.first >= 0 || g.second >= 0) && !missing) {
                nvar +=	addPosition(header, line, d, g);
            }
            else if (missing) {
                missint++;
            }
        }

        if(nfile == 0){
            char *const* s = NULL;
            ohdr = bcf_hdr_subset(header, 0, s, NULL);
        }
        bcf_sr_destroy(reader);
        nfile++;
    }
    cerr << nfile << " files processed"<< endl;
    cerr << nvar << " variants found"<< endl;
    cerr << missint << " missing" << endl;
    htsFile * fp = hts_open(outf, "wb");

    //  bcf_hdr_set(ohdr,hdr_copy.c_str());
    //  bcf_hdr_t *ohdr = bcf_hdr_dup(header);//bcf_hdr_init("w");
    bcf1_t *orec    = bcf_init1();

    bcf_hdr_append(ohdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
    bcf_hdr_append(ohdr, "##INFO=<ID=NPASS,Number=1,Type=Integer,Description=\"Number of ALT genotypes passing\">");

    bcf_hdr_add_sample(ohdr, NULL);      // to update internal structures
    bcf_hdr_write(fp, ohdr);

    //  cout << "header written" << endl;
    map<string, map< unsigned int  ,map< pair<string, string> , marker> > >::iterator it1;
    map<unsigned int, map< pair<string, string> , marker> >::iterator it2;
    map< pair<string, string> , marker>::iterator it3;
    
    const char *allele[2];

    // Write the map out to the vcf
    for (it1 = d.begin(); it1 !=d.end(); ++it1) {
        for (it2 = it1->second.begin(); it2 != it1->second.end(); ++it2) {
            for (it3 = it2->second.begin(); it3!=it2->second.end(); ++it3) {
                bcf_clear1(orec);
                orec->rid = bcf_hdr_name2id(ohdr, it1->first.c_str());
                orec->pos = it2->first;
                allele[0] = it3->first.first.c_str();
                allele[1] = it3->first.second.c_str();
                //	cout << it1->first << "\t"<<it2->first<< "\t"<<it3->first.first<<"\t"<<it3->first.second<<"\t"<<it3->second.npass<<" "<<allele <<endl;
                bcf_update_alleles(ohdr, orec, allele,2);
                bcf_update_info_int32(ohdr, orec, "AC", &(it3->second.AC), 1);	
                orec->qual = it3->second.qual;
                //	bcf_update_info_int32(ohdr, orec, "AN", &(it3->second.AN), 1);	
                bcf_update_info_int32(ohdr, orec, "NPASS", &(it3->second.npass), 1);	
                bcf_write1(fp, ohdr, orec); 
            }
        }
    }

    bcf_hdr_destroy(ohdr);
    bcf_destroy1(orec);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n", argv[2], ret);
        exit(ret);
    }

    cerr << "Finished." << endl;
    return(0);
}
