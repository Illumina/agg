/*  update-qual.c gives a crude re-estimate of QUAL by summing GQ across all alternate genotypes.

    Copyright (C) 2016 Illumina

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

bcf_hdr_t *in_hdr, *out_hdr;
int *gt = NULL, ngt = 0,n;
int *f_gq=NULL,ngq=0;
float sum_gq;


const char *about(void)
{
  return "re-estimates QUAL as sum(GQ) across all ALT genotypes\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  in_hdr  = in;
  out_hdr = out;
  n =  bcf_hdr_nsamples(in_hdr);
  gt = (int *)malloc(2*n*sizeof(int));
  f_gq = (int *)malloc(n*sizeof(int));
  return 0;
}

bcf1_t *process(bcf1_t *rec)
{
  int i;
  if(rec->n_allele!=2) {
    fprintf(stderr,"ERROR: multi-allelic variant found. this is not agg output\n");
    exit(1);
  }
  
  assert(bcf_get_format_int32(in_hdr, rec, "GQ", &f_gq, &ngq)>0);


  bcf_get_genotypes(in_hdr, rec, &gt, &ngt);
  sum_gq=0;
  for(i=0;i<n;i++) {
    if(bcf_gt_allele(gt[2*i])>0 || bcf_gt_allele(gt[2*i+1])>0 ) {
      if(f_gq[i]!=bcf_int32_missing){
	sum_gq+=(float)f_gq[i];
      }
    }        
  }
  rec->qual = sum_gq;
  return rec;
}

void destroy(void)
{
  free(f_gq);
  free(gt);
}
