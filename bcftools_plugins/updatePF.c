/*  updatePF.c updates INFO/PF which is the average number of ALT genotypes that PASSed at this site in the original gvcf

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

bcf_hdr_t *in_hdr, *out_hdr;
int *gt = NULL, ngt = 0,n;
int *f_pf=NULL,npf=0;

const char *about(void)
{
  return "Adds -log10(p-value) for HWE exact test to INFO columns. Bi-allelic sites only at the moment.\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  in_hdr  = in;
  out_hdr = out;
  n =  bcf_hdr_nsamples(in_hdr);
  gt = (int *)malloc(2*n*sizeof(int));
  f_pf = (int *)malloc(n*sizeof(int));
  bcf_hdr_append(out_hdr, "##INFO=<ID=PF,Number=A,Type=Float,Description=\"proportion of genotypes containing an ALT that passed the original single sample gvcf filter\">");

  return 0;
}

bcf1_t *process(bcf1_t *rec)
{
  int i;
  bcf_get_genotypes(in_hdr, rec, &gt, &ngt);
  assert(bcf_get_format_int32(in_hdr, rec, "PF", &f_pf, &npf)>0);
  float i_pf=0.;
  float nalt=0;
  for(i=0;i<n;i++) {
    if(bcf_gt_allele(gt[2*i])>0 || bcf_gt_allele(gt[2*i+1])>0 ) {
      if(f_pf[i]) i_pf++;
      nalt++;
    }      
  }
  //  fprintf(stderr,"%d %f %f\n",rec->pos+1,i_pf,nalt);
  if(i_pf>0.0) i_pf/=nalt;
  bcf_update_info_float(out_hdr, rec, "PF", &i_pf, 1);
  return rec;
}

void destroy(void)
{
  free(f_pf);
  free(gt);
}


