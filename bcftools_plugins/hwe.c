/*  hwe.c exact test for Hardy-Weinberg Equilibrium.  

    See this paper for details:

    Wigginton, Janis E., David J. Cutler, and Gonçalo R. Abecasis. "A note on exact tests of Hardy-Weinberg equilibrium." The American Journal of Human Genetics 76.5 (2005): 887-893.

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
int32_t count[3];
double *p;

float hwe(bcf1_t *rec) {
  assert(rec->n_allele==2);
  float fpval,af;//alternate allele frequency
  float F;//inbreeding coefficeint.
  int an=0,i;
  int ret = bcf_get_genotypes(in_hdr, rec, &gt, &ngt);
  assert(ret==2*n);
  for(i=0;i<3;i++) count[i]=0;

  if(ret!=2*n)
    return(0.);

  for(i=0;i<n;i++) {
    if(gt[i*2]!=bcf_gt_missing && gt[i*2+1]!=bcf_gt_missing) {
      int g = bcf_gt_allele(gt[i*2]);
      if(gt[i*2+1]>0)
	g+=bcf_gt_allele(gt[i*2+1]);
      else
	g*=2;
      //      fprintf(stderr,"g=%d/%d\n",bcf_gt_allele(gt[i*2]),bcf_gt_allele(gt[i*2+1]));
      if(!(g>=0&&g<3)) {
	fprintf(stderr,"bad genotype at pos %d sample %d g=%d/%d\n",rec->pos+1,i,bcf_gt_allele(gt[i*2]),bcf_gt_allele(gt[i*2+1]));
	exit(1);
      }
      count[g]++;
      an++;
    }
  }

  if(an==count[0]||an==count[2])//monomorphic.stop.
    return(0.0);

  af = (float)( count[1] + 2*count[2] ) / (float)(2. * an);
  F = 1. -  ((float)count[1]) / (an * 2. * af * (1-af) );
  bcf_update_info_float(out_hdr, rec, "F", &F, 1);

  /* fprintf(stderr,"pos=%d\n",rec->pos+1); */
  /* fprintf(stderr,"%d %d %d\n",count[0],count[1],count[2]); */
  /* fprintf(stderr,"F=%f\n",F); */
    
  double na = count[0]*2 + count[1];
  double nb = count[2]*2 + count[1];
  double n = count[0]+count[1]+count[2];
  double fa = na/(2*n);
  double fb = nb/(2*n);
  int max_het;
  if(na<nb) max_het=na;
  else max_het=nb;
	      
  int e_idx = 2 * fa * fb * n;
  if(!(e_idx%2==max_het%2)) e_idx++;
  assert(e_idx<=max_het);
  p[e_idx] = 1.0;
  double den = 1.0;
  for(i=e_idx;i>1;i-=2) {
    double nab = i;
    double naa = (na-i)/2;
    double nbb = (nb-i)/2;
    p[i-2] = p[i] * ( nab * (nab-1) ) / ( 4 * (naa+1) * (nbb+1) );
    den +=  p[i-2] ;
  }

  for(i=e_idx;i<max_het;i+=2) {
    double nab = i;
    double naa = (na-i)/2;
    double nbb = (nb-i)/2;
    p[i+2] = p[i] * 4*naa*nbb/( (nab+2)*(nab+1) );
    den +=  p[i+2] ;
  }
  double o=p[count[1]];
  double pval=0;
  for(i=(max_het%2);i<=max_het;i+=2) {
    if(o>=p[i])
      pval += p[i]/den;
  }  

  pval = -log10(pval);
  if(pval<=0.0)    
    pval=0.0;//avoids HWE=-0 formatting
  if(pval>1000.)
    pval=1000.; //stops Inf values.
  //  fprintf(stderr,"pval=%f\n",pval);
  fpval=pval;
  bcf_update_info_float(out_hdr, rec, "HWE", &fpval, 1);
  bcf_update_info_int32(out_hdr, rec, "GN", count, 3);
  return(pval);
}

const char *about(void)
{
  return "Adds -log10(p-value) for HWE exact test to INFO columns. Bi-allelic sites only at the moment.\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
  in_hdr  = in;
  out_hdr = out;
  bcf_hdr_append(out_hdr, "##INFO=<ID=GN,Number=G,Type=Integer,Description=\"genotype counts\">");
  bcf_hdr_append(out_hdr, "##INFO=<ID=HWE,Number=1,Type=Float,Description=\"-log10(pvalue) for HWE exact test\">");
  bcf_hdr_append(out_hdr, "##INFO=<ID=F,Number=1,Type=Float,Description=\"inbreeding coefficient\">");
  n =  bcf_hdr_nsamples(in_hdr);
  ngt=2*n;
  p = (double *)malloc(ngt*sizeof(double));
  gt = (int *)malloc(2*n*sizeof(int));
  return 0;
}

bcf1_t *process(bcf1_t *rec)
{
  if(rec->n_allele==2) {
    hwe(rec);
  }
  
  return rec;
}

void destroy(void)
{
  free(p);
  free(gt);
}


