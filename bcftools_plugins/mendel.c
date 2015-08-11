/*  mendel.c counts mendel errors and phases where possible (trios and bi-allelic sites only)

    Copyright (C) 2015 Illumina

    Author: Jared O'Connell <joconnell@illumina.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <math.h>
#include <getopt.h>

typedef struct {
  char **fid,**id,**dad,**mum;
  int *dadidx,*mumidx,n;
} pedigree;

bcf_hdr_t *in_hdr, *out_hdr;
int *gt = NULL, ngt = 0;
int count[3];
pedigree ped;
int *pmap;

int read_pedigree(char *fname,pedigree *p) {
  int maxl = 10000;
  char line[maxl],*pch;
  FILE *fp = fopen(fname,"r"); 
  int i,j,count1,count2;
  p->n=0;
  if( fp == NULL )
    {
      fprintf(stderr,"Error while opening %s.\n",fname);
      exit(-1);
    }
  while( fgets(line,maxl, fp) !=NULL) {
    //    fprintf(stderr,"%d %s",p->n,line);
    p->n++;
  }
  p->id = (char **)malloc(p->n * sizeof(char *));
  p->mum = (char **)malloc(p->n * sizeof(char *));
  p->dad = (char **)malloc(p->n * sizeof(char *));
  p->dadidx = (int *)malloc(p->n * sizeof(int));
  p->mumidx = (int *)malloc(p->n * sizeof(int));

  fclose(fp);
  fp = fopen(fname,"r"); 

  for(i=0;i<p->n;i++) {
    assert(fgets(line,maxl, fp) !=NULL);
    strtok(line,"\t ");
    pch = strtok(NULL,"\t ");
    p->id[i] = (char *)malloc(strlen(pch)*sizeof(char *));
    strcpy(p->id[i],pch);
    pch = strtok(NULL,"\t ");
    p->dad[i] = (char *)malloc(strlen(pch)*sizeof(char *));
    strcpy(p->dad[i],pch);
    pch = strtok(NULL,"\t ");
    p->mum[i] = (char *)malloc(strlen(pch)*sizeof(char *));
    strcpy(p->mum[i],pch);
  }
  //  for(i=0;i<p->n;i++)   fprintf(stderr,"%d %s %s %s\n",i,p->id[i],p->mum[i],p->dad[i]);

  //index the pedigree. 
  for(i=0;i<p->n;i++) {
    p->mumidx[i]=-1;
    p->dadidx[i]=-1;
    
    if(strcmp(p->mum[i],"0")!=0) 
      for(j=0;j<p->n;j++) 
	if(strcmp(p->mum[i],p->id[j])==0) 
	  p->mumidx[i] = j;

    if(strcmp(p->dad[i],"0")!=0)
      for(j=0;j<p->n;j++)
	if(strcmp(p->dad[i],p->id[j])==0)
	  p->dadidx[i] = j;    
  }  

  count1=0;count2=0;
  for(i=0;i<p->n;i++) {
    if(p->dadidx[i]!=-1&&p->mumidx[i]!=-1)
      count1++;
    else if(p->dadidx[i]!=-1||p->mumidx[i]!=-1)
      count2++;
  }

  fprintf(stderr,"%d trios and %d duos found in %s\n",count1,count2,fname);
  return(p->n);
}

int map_pedigree() {
  int i,j,count=0;
  pmap = (int *)malloc(ped.n * sizeof(int));
  fprintf(stderr,"%d %d\n",ped.n,bcf_hdr_nsamples(in_hdr));
  for(i=0;i<ped.n;i++) {
    pmap[i]=-1;
    for(j=0;j<bcf_hdr_nsamples(in_hdr);j++) {
      //      fprintf(stderr,"%s %s\n",in_hdr->samples[j],ped.id[i]);
      if(strcmp(in_hdr->samples[j],ped.id[i])==0) {
	pmap[i]=j;
	count++;
      }
    }
  }
  fprintf(stderr,"%d samples in pedigree file were in vcf\n",count);
  return(count);
}

int mendel_main(bcf1_t *rec, bcf_hdr_t *hdr) {
  int i;
  int32_t nmendel=0;
  int k0,k1,m0,m1,d0,d1;//haploid gneotypes
  int k,m,d;//diploid genotypes

  for(i=0;i<ped.n;i++) {
    int is_mendel_inconsistent=1;
    if(pmap[i]==-1 || ped.dadidx[i]==-1 || ped.mumidx[i]==-1 )
      continue;

    k0=gt[pmap[i]*2]; k1=gt[pmap[i]*2+1];
    d0=gt[pmap[ped.dadidx[i]]*2]; d1=gt[pmap[ped.dadidx[i]]*2+1];
    m0=gt[pmap[ped.mumidx[i]]*2]; m1=gt[pmap[ped.mumidx[i]]*2+1];

    if(k0>bcf_gt_missing&&k1>bcf_gt_missing&&d0>bcf_gt_missing&&d1>bcf_gt_missing&&m0>bcf_gt_missing&&m1>bcf_gt_missing) {
      k=bcf_gt_allele(k0)+bcf_gt_allele(k1);
      m=bcf_gt_allele(m0)+bcf_gt_allele(m1);
      d=bcf_gt_allele(d0)+bcf_gt_allele(d1);

      assert(k>=0&&k<3);
      assert(d>=0&&d<3);
      assert(m>=0&&m<3);
      //check mendel consistency.
      if(m==0&&d==0)
	if(k==0) is_mendel_inconsistent=0;
      if((m==0&&d==1)||(m==1&&d==0))
	if(k!=2) is_mendel_inconsistent=0;
      if((m==0&&d==2)||(m==2&&d==0))
	if(k==1) is_mendel_inconsistent=0;
      if(m==1&&d==1)
	is_mendel_inconsistent=0;
      if((m==1&&d==2)||(m==2&&d==1))
	if(k!=0) is_mendel_inconsistent=0;      
      if(m==2&&d==2)
	if(k==2) is_mendel_inconsistent=0;
      nmendel+=is_mendel_inconsistent;
      //      if(is_mendel_inconsistent) fprintf (stderr,"%d\t%d x %d -> %d\tERROR\n",rec->pos+1,m,d,k);
      //      else fprintf (stderr,"%d\t%d x %d -> %d\n",rec->pos+1,m,d,k);
      //phase.
      if(!is_mendel_inconsistent) {
	if(!(k==1&&m==1&&d==1) && (k==1||m==1||d==1)) {// can/should we phase this site?
	  int du=-10,dt=-10,mu=-10,mt=-10,kd=-10,km=-10;
	  //phase rules
	  if(m!=1) {
	    mt=m/2;
	    mu=m/2;
	    km=mt;
	    kd=k-km;
	    dt=kd;
	    du=d-dt;
	  }
	  if(d!=1) {
	    dt=d/2;
	    du=d/2;
	    kd=dt;
	    km=k-kd;
	    mt=km;
	    mu=m-mt;
	  }
	  if(k!=1) {
	    kd=k/2;
	    km=k/2;
	    dt=kd;
	    du=d-dt;
	    mt=km;
	    mu=m-mt;	    
	  }
	  assert(k==(kd+km));
	  assert(m==(mu+mt));
	  assert(d==(du+dt));
	  gt[pmap[i]*2] = bcf_gt_phased(kd);
	  gt[pmap[i]*2+1]= bcf_gt_phased(km);
	  gt[pmap[ped.dadidx[i]]*2] = bcf_gt_phased(dt);
	  gt[pmap[ped.dadidx[i]]*2+1] = bcf_gt_phased(du);
	  gt[pmap[ped.mumidx[i]]*2]  = bcf_gt_phased(mt);
	  gt[pmap[ped.mumidx[i]]*2+1] = bcf_gt_phased(mu);
	}
      }
    }
  }
  bcf_update_genotypes(hdr,rec,(void *)gt,bcf_hdr_nsamples(hdr)*2);
  //  fprintf(stderr,"NMENDEL=%d\n",nmendel);
  bcf_update_info_int32(hdr, rec, "NMENDEL", &nmendel, 1);
  return(0);
}


const char *about(void)
{
    return "Rudimentary trio analysis. Counts mendel inconsistencies and phase (where possible). Bi-allelic sites only.\n";
}

char *usage(void)
{
  return "guess\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  in_hdr  = in;
  out_hdr = out;
  bcf_hdr_append(out_hdr, "##INFO=<ID=NMENDEL,Number=1,Type=Integer,Description=\"number of mendel inconsistences observed at this site\">");

  int c;
  char *fname = NULL;

  static struct option loptions[] =
    {
      {"pedigree",1,0,'p'},
      {0,0,0,0}
    };

  while ((c = getopt_long(argc, argv, "p:?h",loptions,NULL)) >= 0)    {
    switch (c) 
      {
      case 'p': fname = optarg; break;
      case 'h':
      case '?':
      default: fprintf(stderr,"%s", usage()); exit(1); break;
      }
  }
  if ( !fname )    {
    fprintf(stderr,"Missing the -p option.\n");
    return -1;
  }
  fprintf(stderr,"n=%d\n", bcf_hdr_nsamples(out_hdr));
  read_pedigree(fname,&ped);
  map_pedigree();
  return 0;
}

bcf1_t *process(bcf1_t *rec)
{
  int ret;
  if(rec->n_allele==2) {
    ret = bcf_get_genotypes(in_hdr, rec, &gt, &ngt);
    assert(ret>0);
    ret = mendel_main(rec,out_hdr);
  }
  
  return rec;
}

void destroy(void)
{
  free(gt);
}


