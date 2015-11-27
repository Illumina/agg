/*  plugins/fill-AN-AC.c -- fills AN and AC INFO fields per group (eg. ethnicity, pipeline, etc).

    This code is derivative of the original bcftools plugins by Petr Danecek.

    Copyright (C) 2015 Illumina
    Copyright (C) 2014 Genome Research Ltd.

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
#include <htslib/khash.h>
#include <getopt.h>


bcf_hdr_t *in_hdr, *out_hdr;
int *arr = NULL, marr = 0;
int nsample;
void *group_dict;
int *group_int;
int ngroup=0;
char **group_names;
int *group_an,*group_ac;
khiter_t k;
int *gt = NULL, ngt = 0;
char **an_tag,**ac_tag;

const int khStrInt = 33;
KHASH_MAP_INIT_STR(khStrInt, int) // setup khash to handle string key with int payload
khash_t(khStrInt) *h;
	     
		  const char *about(void) {
		    return "bcftools +fill-AN-AC-bygroup input.bcf -- -g groups.txt\nFill INFO fields AN and AC by group listed in --g group.txt.\n";
		  }

int read_groups(char *fname) {
  int maxl = 10000;
  char line[maxl];
  FILE *fp = fopen(fname,"r"); 

  if( fp == NULL )    {
    fprintf(stderr,"Error while opening %s.\n",fname);
    exit(-1);
  }

  int nline=0;
  int  absent;
  h = kh_init(khStrInt); // create a hashtable
  while( fgets(line,maxl, fp) !=NULL) {   
    assert(strlen(line)>1);
    line[strlen(line)-1]='\0';
    k = kh_put(khStrInt, h, line, &absent);
    if (absent) {
      kh_value(h, k) = ngroup++;
      kh_key(h, k) = strdup(line);
    }
    k=kh_get(khStrInt,h,line);
    group_int[nline] = kh_value(h,k);
    //    fprintf(stderr,"%d %s %d\n",nline,line,kh_val(h,k));
    nline++;
  } 

  assert(nline==nsample);
  return(0);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  char *fname = NULL;
  int c;
  static struct option loptions[] =
    {
      {"groups",1,0,'g'},
      {0,0,0,0}
    };
  
  while ((c = getopt_long(argc, argv, "p:?h",loptions,NULL)) >= 0)    {
    switch (c) 
      {
      case 'g': fname = optarg; break;
      case 'h':
      case '?':
      default: fprintf(stderr,"%s", about()); exit(1); break;
      }
  }
  if ( !fname )    {
    fprintf(stderr,"Missing the -p option.\n");
    return -1;
  }
 
  in_hdr  = in;
  out_hdr = out;
  nsample =  bcf_hdr_nsamples(in_hdr);
  gt = (int *)malloc(2*nsample*sizeof(int));
  fprintf(stderr,"nsample=%d\n", bcf_hdr_nsamples(out_hdr));
  group_int = (int *)malloc(nsample*sizeof(int));
  read_groups(fname);   


  bcf_hdr_append(out_hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">");
  bcf_hdr_append(out_hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">");


  fprintf(stderr,"# of distinct groups: %d (", kh_size(h));
  assert(kh_size(h)==ngroup);
  group_ac = (int *)malloc(ngroup*sizeof(int));
  group_an = (int *)malloc(ngroup*sizeof(int));
  group_names = (char **)malloc(ngroup*sizeof(char *));

  int i=0;
  an_tag = (char **)malloc(sizeof(char *)*ngroup);
  ac_tag = (char **)malloc(sizeof(char *)*ngroup);
  for (k = 0; k < kh_end(h); ++k) {
    if (kh_exist(h, k))  {
      int idx=kh_val(h,k);
      group_names[idx] = kh_key(h,k);
      fprintf(stderr,"%s,",       group_names[idx]  );
      an_tag[idx]=(char *)malloc(1000*sizeof(char));
      ac_tag[idx]=(char *)malloc(1000*sizeof(char));
      sprintf(an_tag[idx], "##INFO=<ID=%s_AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes for group %s\">",group_names[idx],group_names[idx]);
      bcf_hdr_append(out_hdr, an_tag[idx]);
      sprintf(ac_tag[idx], "##INFO=<ID=%s_AC,Number=A,Type=Integer,Description=\"Alternate allele count for group %s\">",group_names[idx],group_names[idx]);
      bcf_hdr_append(out_hdr, ac_tag[idx]);
      sprintf(ac_tag[idx],"%s_AC",group_names[idx]);
      sprintf(an_tag[idx],"%s_AN",group_names[idx]);
      bcf_hdr_append(out_hdr, ac_tag[idx]);
      i++;
    }
  }
  fprintf(stderr,")\n");

  return 0;
}

bcf1_t *process(bcf1_t *rec)
{
  int i;
  assert(rec->n_allele==2);
  hts_expand(int,rec->n_allele,marr,arr);
  int ret = bcf_calc_ac(in_hdr,rec,arr,BCF_UN_FMT);
  for(i=0;i<ngroup;i++) {
    group_ac[i]=0;
    group_an[i]=0;
  }
  if (ret>0)  {
    int  an = 0;
    for (i=0; i<rec->n_allele; i++) an += arr[i];
    bcf_update_info_int32(out_hdr, rec, "AN", &an, 1);
    bcf_update_info_int32(out_hdr, rec, "AC", arr+1, rec->n_allele-1);
  }
  ret = bcf_get_genotypes(in_hdr, rec, &gt, &ngt);
  assert(ret==2*nsample);
  for(i=0;i<nsample;i++){
    if(!bcf_gt_is_missing(gt[2*i])&&!bcf_gt_is_missing(gt[2*i+1])){
      group_an[group_int[i]]+=2;
      group_ac[group_int[i]]+=(bcf_gt_allele(gt[2*i])+bcf_gt_allele(gt[2*i+1]));
    }
  } 

  for(i=0;i<ngroup;i++) {
    bcf_update_info_int32(out_hdr, rec, an_tag[i], group_an+i, 1);
    bcf_update_info_int32(out_hdr, rec, ac_tag[i], group_ac+i, 1);    
  }

  return rec;
}

void destroy(void)
{
  khiter_t k;
  for (k = 0; k < kh_end(h); ++k)
    if (kh_exist(h, k))
      free((char*)kh_key(h, k));
  kh_destroy(khStrInt, h);
  free(group_int);
  free(arr);
}
