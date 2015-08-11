#include "agg.h"

double hwe(vector<unsigned int> count) {
  double na = count[0]*2 + count[1];
  double nb = count[2]*2 + count[1];
  double n = count[0]+count[1]+count[2];
  if(n==0) return(0.0);
  double fa = na/(2*n);
  double fb = nb/(2*n);
  /*
    double *E[3];
    E[0] = fa*fa*n;
    E[1] = 
    E[2] = fb*fb*n;
  */
  int max_het;
  if(na<nb) max_het=(int)na;
  else max_het=(int)nb;

  int e_idx = (int) (2 * fa * fb * n);
  if(!(e_idx%2==max_het%2)) e_idx++;
  //  cout << e_idx << " " << max_het << endl;
  assert(e_idx<=max_het);
  vector<double> p(max_het+1,0);
  p[e_idx] = 1.0;
  double den = 1.0;
  for(int i=e_idx;i>1;i-=2) {
    //    cout <<"i="<< i << " " << den << endl;
    double nab = i;
    double naa = (na-i)/2;
    double nbb = (nb-i)/2;
    //    assert( (naa+nab+nbb)==n);
    //    assert( (i-2) > -1 );
    p[i-2] = p[i] * ( nab * (nab-1) ) / ( 4 * (naa+1) * (nbb+1) );
    den +=  p[i-2] ;
  }

  for(int i=e_idx;i<max_het;i+=2) {
    //    cout <<"i="<< i << " " << den << endl;
    double nab = i;
    double naa = (na-i)/2;
    double nbb = (nb-i)/2;
    //    assert( (naa+nab+nbb)==n);
    //    assert( (i+2) < (1+max_het) );
    p[i+2] = p[i] * 4*naa*nbb/( (nab+2)*(nab+1) );
    den +=  p[i+2] ;
  }
  double o=p[count[1]];
  double pval=0;
  for(int i=(max_het%2);i<=max_het;i+=2) {
    //    cout << i << " " << o << " " << p[i]<<"/"<<den<<endl;
    if(o>=p[i])
      pval += p[i]/den;
  }  
  //  return(pval);
  pval = -log10(pval);
  if(pval<0.0)
    pval=0.0;
  return(pval);
}

//returns 0 if  a compress block, END (which is >0) if is a block
int isBlock(bcf_hdr_t *header,bcf1_t *line) {
  int end=-1;
  int *endptr = &end;
  int ndst=1;
  int val=bcf_get_info_int32(header,line,"END",&endptr,&ndst);
  if(val==1) 
    return(end);
  else
    return(0);
}

int addPosition(bcf_hdr_t *header,bcf1_t *line, map < string , map < unsigned int, map < pair<string,string> , marker > > >  & d, pair<byte,byte> g) {
  // Adds a new position to the map----------map < string , map < unsigned int, map < pair<string,string> , marker > > >
  //------------------------------------------------Chr20      ->       pos       ->        ref,   alt ->   marker
  //  cout << isNonVar(header,line) <<endl;

  if(isBlock(header, line) > 0)    {
    //cout << "Skipping "<< chrom << " " << pos+1 << endl;    
    return(0);//block ref/ref segment. skip.  
  }

  string chrom = bcf_seqname(header, line) ;
  unsigned int pos = line->pos ;                  //bcf coordinates are 0-based
  marker m(header, line);
  m.AN = 0;
  m.npass = bcf_has_filter(header, line, ".");
  int  count = 0;
  // Find the position and add the detail if it exists
  for(unsigned int i = 1; i < line->n_allele; i++) { 
    m.AC = 0;
    bool seen = false;
    pair<string,string> key(line->d.allele[0], line->d.allele[i]);          // ref, alt
    if(g.first == i) m.AC++;                                                // Allele count
    if(g.second == i) m.AC++;
    if(d.count(chrom)) {
      if(d[chrom].count(pos)) {
	if(d[chrom][pos].count(key)) {
	  d[chrom][pos].at(key).AC += m.AC;
	  d[chrom][pos].at(key).npass += m.npass;
	  if(m.qual > d[chrom][pos].at(key).qual)  d[chrom][pos].at(key).qual = m.qual;       //take max qual at this site.
	  seen = true;
	}
      }
    }

    // Create a new entry if it doesn't exist
    if(!seen)  {
      // cout << chrom << " " << pos+1 << " " << line->d.allele[0] << " " << line->d.allele[i] << endl;
      d[chrom][pos].insert(pair< pair<string,string>, marker > (key, m));
      count++;
    }
  }
  return(count);
}

vector<region> getRegions(char *fname, int buf) {

  bcf_srs_t *sr = bcf_sr_init();
  if(!bcf_sr_add_reader(sr,fname)) 
    die("problem opening" + (string)fname);
  vector<region> rgs;
  int count = 0;
  while ( bcf_sr_next_line(sr) ) {
    bcf1_t *line = bcf_sr_get_line(sr,0); 
    string chrom = bcf_seqname(sr->readers[0].header,line);   
    uint32_t pos = line->pos+1;//bcf coordinates are 0-based
    if(count==0||rgs[count-1].chrom!=chrom) {
      region r(chrom,pos-buf,pos);
      rgs.push_back(r);
      count++;
    }  else {
      rgs[count-1].stop=pos+buf;
    }
  }
  bcf_sr_destroy(sr);  
  return(rgs);
}

// marker::marker() {
//   ref="";alt="";qual=-1.;
// }

marker::marker(bcf_hdr_t *header,bcf1_t *line) {
  qual = line->qual;
  npass=0;            //bcf_has_filter(header, line, ".") ;
  AN=0;
  AC=0;
  counts.assign(3,0);
}


marker::~marker() {

}

string region::toString() {
  stringstream sstm;
  sstm<<chrom<<":"<<start<<"-"<<stop;
  return(sstm.str().c_str());  
}

region::region(const string& c,uint32_t a,uint32_t b) {
  chrom = c;
  start=a;
  stop=b;
}

homBlock::homBlock(uint32_t a,uint32_t b,bool p,bool m) {
  start=a;
  stop=b;
  pass=p;
  missing=m;
}

//iterates through a sites file and builds our site list.
int buildSiteList(const string& fname, const string& region,  vector<marker> & out) {
  //  cerr << "Getting variants from "<<fname<<" "<<region<<endl;
  bcf_srs_t *rdr=bcf_sr_init() ; 

  if(bcf_sr_set_regions(rdr,region.c_str(),0)!=0)
    die("problem setting regions");

  if(!(bcf_sr_add_reader (rdr, fname.c_str()) ))
    die("problem opening "+fname);

  bcf1_t *line;
  int nvar =0;
  while (bcf_sr_next_line (rdr)) {
    line = bcf_sr_get_line(rdr,0);
    if(line->n_allele!=2) {
      cerr << endl << line->pos+1 <<endl;
      die("sites had a line with n_allele!=2");
    }

    string chrom = bcf_seqname(rdr->readers[0].header,line) ;
    marker m(rdr->readers[0].header,line);
    m.pos = line->pos ;
    m.ref = line->d.allele[0];
    m.alt = line->d.allele[1];
    out.push_back(m);
    nvar++;
  }
  bcf_sr_destroy(rdr);
  cerr << "Found "<< nvar << " variants in region." <<endl;

  return(nvar);
}

int writeSiteList(const string& fname, map<string,map<unsigned int,map< pair<string,string> ,marker> > >  & d) {
  htsFile *fp = hts_open(fname.c_str(),"wz");

  bcf_srs_t *reader = bcf_sr_init();

  bcf1_t *orec = bcf_init1();
  bcf_hdr_t *hdr;
  if(!bcf_sr_add_reader(reader, fname.c_str()))
    die("problem opening "+ fname);
  else
    hdr = bcf_hdr_dup(reader->readers[0].header);
  bcf_sr_destroy(reader);
  
    bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
  bcf_hdr_write(fp, hdr);
  //  cout << "header written"<<endl;
  map<string,map<unsigned int,map< pair<string,string> ,marker> > >::iterator it1;
  map<unsigned int,map< pair<string,string> ,marker> >::iterator it2;
  map< pair<string,string> ,marker>::iterator it3;
  const char *allele[2];
  for (it1=d.begin();it1!=d.end();++it1) {
    for (it2=it1->second.begin();it2!=it1->second.end();++it2) {
      for (it3=it2->second.begin();it3!=it2->second.end();++it3) {
	bcf_clear1(orec);
	orec->rid = bcf_hdr_name2id(hdr, it1->first.c_str());
	orec->pos = it2->first;
	allele[0] = it3->first.first.c_str();
	allele[1] = it3->first.second.c_str();
	//cout << it1->first << "\t"<<it2->first<< "\t"<<it3->first.first<<"\t"<<it3->first.second<<"\t"<<it3->second.npass<<" "<<allele <<endl;
	bcf_update_alleles(hdr, orec, allele, 2);
	bcf_update_info_int32(hdr, orec, "AC", &(it3->second.AC), 1);	
	orec->qual = it3->second.qual;
	//	bcf_update_info_int32(hdr, orec, "AN", &(it3->second.AN), 1);	
	bcf_update_info_int32(hdr, orec, "NPASS", &(it3->second.npass), 1);	
	bcf_write1(fp, hdr, orec);
      }
    }
  }

  bcf_hdr_destroy(hdr);
  bcf_destroy1(orec);
  int ret;
  if ( (ret=hts_close(fp)) )
    {
      fprintf(stderr,"hts_close(%s): non-zero status %d\n", fname.c_str(),ret);
      exit(ret);
    }
  cerr << "done" << endl;
  return(0);
}
