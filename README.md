agg: a utility for aggregating gvcfs. This software is not commercially supported.

Copyright (c) 2016, Illumina, Inc. All rights reserved. 

The agg source code is provided under the [GPLv3 license] (LICENSE).

###Summary

This tool implements a basic pipeline to merge Illumina gvcfs in a dynamic fashion. That is, not all gvcfs need to merged at once, new groups of samples can be added periodically.  It achieves this by storing variants in a standard bcf and storing depth/GQ information in an auxilliary file (.dpt - depth track). The variants and depth files make up an "agg chunk" and these chunks can then be merged and genotyped, with the depth file allowing the union of variants to be genotyped across all samples.  Note, here "genotyping" simply means we can gauge the coverage/GQ at a variant position for samples that do not have that variant, allowing us to examine the evidence that this sample was homozygous reference at that location.

For example, if you had 3000 gvcfs, your pipeline might be:

1. Create 6 agg chunks of samples 500 (ingestion)
2. Create 3000 sample vcf/bcf from these 6 chunks (genotyping)

The advantage of this approach is that when another 500 samples come along, you only have to build one new chunk (step 1) and then re-genotype from the 6+1=7 chunks to create a 3500 sample bcf (step 2). This is faster than merging from all gvcfs every time some new samples arrive.

**Note:** if you are working with relatively few gvcfs, then [gvcftools](https://github.com/sequencing/gvcftools) merge_variants is more appropriate.  The agg pipeline is more complex to run and is geared towards 1000s of gvcfs.

###Installation
The only compilation dependency is [htslib](http://www.htslib.org/) which is included with the software.  

```
git clone git@github.com:BEETL/agg.git
cd agg/
make
./agg
```

It should be noted that parts of the agg source code were taken directly from the excellent [bcftools](https://github.com/samtools/bcftools) which is licensed permissively under BSD.

###Building an agg chunk
The input to agg's genotyping routine is one or more agg "chunks".  As new batches of samples arrive they can be rolled into a new chunk, without the need to modify previous chunks containing older samples. 

The easiest way to make a chunk is with the provided python script which wraps several agg commands. Say we have a plain text list of a few thousand gvcfs in `gvcfs.txt`. We will build chunks of size 500 (watch out for file handle limits), so first split your gvcf list via:
```
split -d -l 500 gvcfs.txt chunk_
```
then run the script on each chunk
```
mkdir chunks/
for i in chunk_*;
do
        python ~/agg/make_chunk.py $i -o chunks/${i} -nproc 16 -tmp /scratch/
done        
```
Note the above command will:
* use 16 processes
* take a while (~24 hours)
* put a lot of temporary files on /scratch/ (~1.5GB per sample)

In practice, a user would want to submit each of these commands to a cluster node with multiple cores and sufficient local scratch.

The `make_chunk.py` script simply wraps some agg commands for ease-of-use. Users may be able to design more efficient bespoke pipelines for their respective systems. For a description of how to do this manual see [doc/ingest.md](doc/ingest.md)

###Genotyping and merging agg chunks
Once you have your chunks, life is easy.  Simply call `agg genotype` on any number of chunks to produce a typical multi-sample bcf/vcf that contains all the samples in all the chunks genotyped at all variants seen across the chunks. 
```
$ agg genotype -r chr1 chunk_*.bcf -Ob -o merged.chr1.bcf
$ bcftools index merged.chr1.bcf
```
Note you can (optionally) use the `-r` argument to specify chromosomes or smaller regions for easy parallelism, again `xargs` is your friend:
```
#genotype autosomes separately
$ for i in {1..22};
  do 
     echo genotype -r chr${i} chunk_*.bcf -Ob -o merged.chr${i}.bcf;
  done | xargs agg -l -P 16

#concatenate autosomes into one big file if desired
$ for i in {1..22};do echo merged.chr${i}.bcf;done > files_to_concat.txt
$ bcftools concat -f files_to_concat.txt -Ob -o merged.bcf
bcftools index merged.bcf
```
**Note:** Variants are not normalised by agg, meaning occasionally you could see the the same indel with different representations.  Using a normalisation tool such as [vt](https://github.com/atks/vt) on agg output is a good idea.

####Filtering
The output from `agg` is very raw, containing all variants called in any sample, whether they passed filter in the single sample gvcfs or not. How exactly to filter these down to a high quality list of variants is a research topic in itself.  A simplistic first pass may involve:

* set genotypes where GQ<10 to missing
* remove sites where < 90% of genotypes are called (after GQ<10 removal)
* remove sites where QUAL<30

this can be achieved with something like:
```
bcftools filter -e 'FMT/GQ<10' -S . -O u | bcftools annotate -x FILTER -Ou| bcftools view -i 'QUAL>=30 & AN>900' -Ob -o merged.flt.bcf
```
This is very crude, typically one may also filter on extreme depth, allelic imbalance, divergence from HWE etc etc.

####Creating a site list
For applications such as annotating variants in an individual with a rare disease.  Often all that is needed is a site-only vcf with summary statistics of interest (such as allele frequency) stored in the INFO field.  This is straightforward to generate from the multi-sample bcf that was created in the previous section.
```
bcftools view -G merged.flt.bcf -Ou |  bcftools +fill-AN-AC | bcftools view -Oz -o merged.sites.vcf.gz
tabix merged.sites.vcf.gz
```
We may also wish to add some custom stuff to the INFO field. For example, the hwe.c plugin included with this package will add the -log10(p-value) for Fisher's Exact test for divergence from Hardy-Weinberg Equilibrium, as well as the inbreeding coefficient.
```
bcftools view merged.flt.bcf -Ou | bcftools +fill-AN-AC | bcftools +hwe | bcftools view -G -Oz -o merged.sites.vcf.gz
tabix merged.sites.vcf.gz
```

###A note on genotyping homref positions
Genotyping an individual (from their gvcf) who does not have an ALT allele called at SNP location is relatively easy (at least I think so). We simply take the DP and the homref GQ at that base. For indels, things are not so straightforward, I have implemented what I think is a reasonable scheme.  For deletions, agg reports the average depth across the length of the deletion and the minimum homref GQ observed. For insertions, agg reports the average depth of the two bases flanking the insertion and the minimum homref GQ of these two bases. This is of course inferior to proper joint calling where reads are aligned to candidate haplotypes to generate a likelihood for each possible genotype. I would be happy to hear about better alternatives to these rules.

###Known issues
Multi-allelic annotations for individuals that are ALT in one (or more) of the alleles might be slightly off. Why? Currently DP/GQ are filled from the .dpt file for an individual if they do not have the variant being processed *regardless* of whether they have a different variant at the same position.  This could be solved by buffering all co-occurring variants which is something we are looking at for future versions.

Currently, one should "trust" the annotations where the variant is typed as ALT over ones where they are typed as HOMREF at the same position.


