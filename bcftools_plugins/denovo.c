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
//GT:GQ:DP:DPF:AD:PF
int *gt = NULL, ngt = 0;
int32_t *gq,ngq=0;
int mingq=10;//"denovos" where the minimum pedigree gq is below this will not be considered.
int32_t *dp,ndp=0;
int32_t *dpf,ndpf=0;
int32_t *ad,nad=0;
int32_t *ft,nft=0;



int count[3];
pedigree ped;
int *pmap;

int read_pedigree(char *fname,pedigree *p) {
  int maxl = 10000;
  char line[maxl],*pch;
  FILE *fp = fopen(fname,"r"); 
  int i,j,count1,count2;
  p->n=0;
  if( fp == NULL )    {
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
    if(pch==NULL) {
      p->n=i;
      break;
    }
      
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


int mendel_main(bcf1_t *rec, bcf_hdr_t *hdr) 
{
    int i,j;
    int k0,k1,m0,m1,d0,d1;//haploid gneotypes
    int k,m,d;//diploid genotypes

  for(i=0;i<ped.n;i++) {
    if(pmap[i]==-1 || ped.dadidx[i]==-1 || ped.mumidx[i]==-1 )
      continue;

    k0=gt[pmap[i]*2]; k1=gt[pmap[i]*2+1];
    d0=gt[pmap[ped.dadidx[i]]*2]; d1=gt[pmap[ped.dadidx[i]]*2+1];
    m0=gt[pmap[ped.mumidx[i]]*2]; m1=gt[pmap[ped.mumidx[i]]*2+1];

    k=0;
    m=0;
    d=0;
    /* fprintf(stdout,"%d/%d\n",bcf_gt_allele(k0),bcf_gt_allele(k1)); */
    /* fprintf(stdout,"%d %d %d\n",bcf_gt_is_missing(k0),bcf_gt_is_missing(k1),!bcf_gt_is_missing(k0)&&!bcf_gt_is_missing(k1)); */
    /* fprintf(stdout,"%d %d %d\n",k,d,m); */
    if(!bcf_gt_is_missing(k0)&&!bcf_gt_is_missing(k1))
	k=bcf_gt_allele(k0)+bcf_gt_allele(k1);
    if(!bcf_gt_is_missing(m0)&&!bcf_gt_is_missing(m1))
	m=bcf_gt_allele(m0)+bcf_gt_allele(m1);
    if(!bcf_gt_is_missing(d0)&&!bcf_gt_is_missing(d1))
	d=bcf_gt_allele(d0)+bcf_gt_allele(d1);
    //fprintf(stdout,"%d %d %d\n",k,d,m);
    assert(k>=0&&k<3);
    assert(d>=0&&d<3);
    assert(m>=0&&m<3);
    //check mendel consistency.

    if(m==0&&d==0&&k>0)
    {
	if(k!=0)
	{ //denovo
//GT:GQ:DP:DPF:AD:PF
//	    fprintf(stdout,"%d %d %d\n",k,d,m);
	    bcf_unpack(rec, BCF_UN_FMT);
	    assert(bcf_get_format_int32(hdr,rec,"GQ",&gq,&ngq)>0);
	    if(gq[pmap[i]]>=mingq && gq[pmap[ped.dadidx[i]]]>=mingq &&  gq[pmap[ped.mumidx[i]]]>=mingq)
	    {
		assert(bcf_get_format_int32(hdr,rec,"DP",&dp,&ndp)>0);
		assert(bcf_get_format_int32(hdr,rec,"DPF",&dpf,&ndpf)>0);
		assert(bcf_get_format_int32(hdr,rec,"AD",&ad,&nad)>0);
		assert(bcf_get_format_int32(hdr,rec,"PF",&ft,&nft)>0);
		fprintf(stdout,"%s\t%d\t%s\t%s\t%s\t%s\t%s",bcf_hdr_id2name(hdr,rec->rid),rec->pos+1,rec->d.allele[0],rec->d.allele[1],hdr->samples[pmap[i]],hdr->samples[pmap[ped.dadidx[i]]],hdr->samples[pmap[ped.mumidx[i]]]);
		int pedlook[3];
		pedlook[0] = pmap[i];
		pedlook[1] = pmap[ped.dadidx[i]];
		pedlook[2] = pmap[ped.mumidx[i]];
		kstring_t str = {0,0,0};
		for(j=0;j<3;j++)
		{
		    int idx = pedlook[j];
//		fprintf(stderr,"idx=%d\n",idx);
		    kputc('\t',&str);
		    if(!bcf_gt_is_missing(gt[idx*2]))	 kputw(bcf_gt_allele(gt[idx*2]),&str);else kputc('.',&str);
		    kputc('/',&str);
		    if(!bcf_gt_is_missing(gt[idx*2+1]))	 kputw(bcf_gt_allele(gt[idx*2+1]),&str);else kputc('.',&str);
		    kputc(':',&str);
		    if(gq[idx]>=0)	 kputw(gq[idx],&str);else kputc('.',&str);
		    kputc(':',&str);
		    if(dp[idx]>=0)	 kputw(dp[idx],&str);else kputc('.',&str);
		    kputc(':',&str);
		    if(dpf[idx]>=0)	 kputw(dpf[idx],&str);else kputc('.',&str);
		    kputc(':',&str);
		    if(ad[2*idx]>=0)	 kputw(ad[2*idx],&str);else kputc('.',&str);
		    kputc(',',&str);
		    if(ad[2*idx+1]>=0)	 kputw(ad[2*idx+1],&str);else kputc('.',&str);
		    kputc(':',&str);
		    if(ft[idx]>=0)	 kputw(ft[idx],&str);else kputc('.',&str);
		}
		fprintf(stdout,"%s",str.s);
		free(str.s);
		fprintf(stdout,"\n");
	    }
	}
    }
  }
  return(0);
}



char *usage(void)
{
  return "prints out nominal denovo mutations\n";
}

char *about(void)
{
  return "prints out nominal denovo mutations\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  in_hdr  = in;

  int c;
  char *fname = NULL;

  static struct option loptions[] =
    {
      {"pedigree",1,0,'p'},
      {"min-gq",1,0,'g'},
      {0,0,0,0}
    };

  while ((c = getopt_long(argc, argv, "g:p:?h",loptions,NULL)) >= 0)    {
    switch (c) 
      {
      case 'p': fname = optarg; break;
      case 'g': mingq = atoi(optarg); break;
      case 'h':
      case '?':
      default: fprintf(stderr,"%s", usage()); exit(1); break;
      }
  }
  if ( !fname )    {
    fprintf(stderr,"Missing the -p option.\n");
    return -1;
  }
  fprintf(stderr,"n=%d\n", bcf_hdr_nsamples(in_hdr));
  read_pedigree(fname,&ped);
  map_pedigree();
  return 1;
}

bcf1_t *process(bcf1_t *rec)
{
  int ret;
  if(rec->n_allele==2) {
    ret = bcf_get_genotypes(in_hdr, rec, &gt, &ngt);

    assert(ret>0);
    ret = mendel_main(rec,in_hdr);
  }
  
  return(NULL);
}

void destroy(void)
{
  free(gt);
}


