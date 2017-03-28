/*  update-info.c updates INFO fields generated by agg (DP,AD,PF,AN,AC). This is useful when you subset some agg output and want to refresh annotations accordingly.

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
int *f_pf=NULL,npf=0,*f_dp,*f_dpf,ndp=0,*f_ad,nad=0,*f_gq=NULL,ngq=0;


const char *about(void)
{
    return "fills FORMAT/AD from FORMAT/DP\n";
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    in_hdr  = in;
    out_hdr = out;
    n =  bcf_hdr_nsamples(in_hdr);
    gt = (int *)malloc(2*n*sizeof(int));
    f_pf = (int *)malloc(n*sizeof(int));
    f_dp = (int *)malloc(n*sizeof(int));
    f_dpf = (int *)malloc(n*sizeof(int));
    f_ad = (int *)malloc(2*n*sizeof(int));
    return 0;
}


bcf1_t *update_info_fields(bcf1_t *rec,bcf_hdr_t *_in_hdr, bcf_hdr_t *_out_hdr)
{
    int i;
    if(rec->n_allele!=2) 
    {
	fprintf(stderr,"ERROR: multi-allelic variant found. this is not agg output\n");
	exit(1);
    }
    int nval=bcf_get_format_int32(_in_hdr, rec, "DP", &f_dp, &ndp);
    assert(nval==n);
    assert(bcf_get_format_int32(_in_hdr, rec, "AD", &f_ad, &nad)>0);
    assert(bcf_get_format_int32(_in_hdr, rec, "GQ", &f_gq, &ngq)==n);
  
    bcf_get_genotypes(_in_hdr, rec, &gt, &ngt);

    for(i=0;i<n;i++) 
    {
	if(f_gq[i]==bcf_int32_missing)
	{
	    f_gq[i]=0; 
	}
	if(f_dp[i]==bcf_int32_missing)
	{
	    f_dp[i]=0;
	}
	if( (bcf_gt_allele(gt[2*i])==0&&bcf_gt_allele(gt[2*i+1])==0) || (gt[i*2+1]==bcf_gt_missing && gt[i*2]==bcf_gt_missing) )
	{
	    if(f_ad[i*2]==bcf_int32_missing || f_ad[i*2+1]==bcf_int32_missing)
	    {
		f_ad[i*2]=f_dp[i];
		f_ad[i*2+1]=0;
	    }
		
	}
    }
    bcf_update_format_int32(out_hdr,rec,"GQ",f_gq,n);
    bcf_update_format_int32(out_hdr,rec,"AD",f_ad,2*n);
    return rec;
}

void destroy(void)
{
    free(f_dp);
    free(f_ad);
    free(f_pf);
    free(f_gq);
    free(gt);
}

bcf1_t *process(bcf1_t *rec)
{
    return(update_info_fields(rec,in_hdr,out_hdr)) ;
}
