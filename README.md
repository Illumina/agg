agg: a utility for aggregating gvcfs

Copyright (c) 2015, Illumina, Inc. All rights reserved. See LICENSE.pdf for further details.

####Summary

This tool implements a basic pipeline to merge Illumina gvcfs in a dynamic fashion. That is, not all gvcfs need to merged at once, new groups of samples can be added periodically.  It achieves this by storing variants in a standard bcf and storing depth/GQ information in an auxilliary file (.dpt - depth track). The variants and depth files make up an "agg chunk" and these chunks can then be merged on the fly, with the depth file allowing the union of variants to be genotyped across all samples.  Note, here "genotyping" simply means we can gauge the coverage/GQ at a variant position for samples that do not have that variant, allowing us to examine the evidence that this sample was homozygous reference at that location.

If you are working with relatively few gvcfs, then [gvcftools](https://github.com/sequencing/gvcftools) merge_variants is more appropriate.  The agg pipeline is more complex to run and is geared towards 1000s of gvcfs.

####Installation
The only compilation dependency is [htslib](http://www.htslib.org/) which is included with the software.  

```
git clone git@github.com:BEETL/agg.git
cd agg/
make
```

It should be noted that parts of the agg source code were taken directly from the excellent [bcftools](https://github.com/samtools/bcftools) which is licensed permissively under BSD.

####Usage

```
./agg
Program:        agg 75c2488 (aggregation tool for multiple samples)
Contact:        joconnell@illumina.com

Copyright (c) 2015, Illumina, Inc. All rights reserved. See LICENSE.pdf for further details.

Usage:  agg <command> [options]

Commands:
        ingest1         converts gvcfs to input suitable for agg ingest2
        ingest2         uses output files from ingest1 to build an agg chunk
        genotype        genotypes and merges agg chunks into a multi-sample bcf/vcf
```

#####Building an agg chunk
The input to agg's genotyping routine is one or more agg "chunks".  As new batches of samples arrive they can be rolled into a new chunk, without the need to modify previous chunks containing older samples. 

Creating a chunk is currently a two stage process using the `agg ingest1` and `agg ingest2` commands.  The individual gvcfs are first pre-processed with `ingest1` and then merged into a chunk with `agg ingest2`.  In the future, we plain to wrap these these two steps into a single (faster) `ingest` command.

######ingest1: pre-process gvcfs
You can do this to one gvcf like so:
```
$ mkdir ingest1/
$ agg ingest1 sample1.genome.vcf.gz -o ingest1/sample1
Input: sample1.genome.vcf.gz    Output: ingest1/sample1
depth: ingest1/sample1.dpt
variants: ingest1/sample1.bcf
Indexing ingest1/sample1.bcf
Done.
```
this takes ~10 minutes.

In practice, we have 100s to 1000s of gvcfs and want to automate things. Say you have a 16 core node and want to build an agg chunk from 500 gvcfs listed in `gvcfs.txt`.  We can leverage the multiple CPUs using unix `xargs`.
```
$ for i in `cat gvcfs.txt`;do out=$(basename ${i%.genome.vcf.gz});echo ingest1 $i -o ingest1/${out};done | xargs -p 16 -n 4 agg
```
note you can replace `cat gvcfs.txt` with `find . -name '*.genome.vcf.gz` or similar.

These files are the input for `agg ingest2`, which builds a chunk and is explained next. These intermediate files are rather large (1GB-2GB per sample), but after building a chunk they can be disposed of. You may wish to store them on local scratch space if it is available.

######ingest2: merge temporary files into a chunk
Let's roll how intermediate files into five separate agg chunks. We will use `xargs` to leverage multiple cores again.
```
$ find ingest1/ -name '*.bcf' > ingest1.txt ##get bcfs from ingest1 step
$ split -d -l 100 ingest1.txt chunk_ #splits the file into 5 groups of 100
$ ls chunk_*
chunk_00  chunk_01  chunk_02  chunk_03  chunk_04
$ for i in chunk_*;do echo ingest2 -l $i -@4 -o $i;done | xargs -n 7 -P 16
...
$ ls chunk_*.*
```

#####Genotyping and merging agg chunks
Once you have your chunks, life is easy.  Simply call `agg genotype` on any number of chunks to produce a typical multi-sample bcf/vcf that contains all the samples in all the chunks fully genotyped at all variants seen in all the chunks. 
```
$ agg genotype -r chr1 chunk_00.bcf chunk_01.bcf chunk_02.bcf chunk_04.bcf -Ob -o merged.chr1.bcf
$ bcftools index merged.chr1.bcf
```
Note you can (optionally) use the `-r` argument to specify chromosomes or smaller regions for easy parallelism, again `xargs` is your friend:
```
#genotype autosomes separately
$ for i in {1..22};do echo genotype -r chr${i} chunk_00.bcf chunk_01.bcf chunk_02.bcf chunk_04.bcf -Ob -o merged.chr${i}.bcf;done | xargs agg -n 10 -P 16

#concatenate autosomes into one big file if desired
$ for i in {1..22};do echo merged.chr${i}.bcf;done > files_to_concat.txt
$ bcftools concat -f files_to_concat.txt -Ob -o merged.bcf
bcftools index merged.bcf
```

######Filtering
The output from `agg` is very raw, containing all variants called in any sample, filtered or not. How exactly to filter this down to a high quality list of variants is a research topic in itself.  A simplistic first pass may involve:

* set genotypes where GQ<10 to missing
* remove sites where only 50% of genotypes are called (after GQ<10 removal)
* remove sites where QUAL<30

this can be achieved by:
```
bcftools filter -e 'FMT/GQ<10' -S . -O u | bcftools view -i 'QUAL>=30 & AN>500' -Ob -o merged.flt.bcf
```
This is very crude, typically one may also filter on extreme depth, allelic imbalance, divergence from HWE etc etc.

#####Creating a site list
For applications such as annotating variants in a rare disease study.  Often all that is needed is a site-only vcf with summary statistics of interest (such as allele frequency). stored in the INFO field.  This is straightforward to generate from the multi-sample bcf that was created in the previous section.
```
bcftools view -G merged.flt.bcf -Oz -o merged.sites.vcf.gz
tabix merged.sites.vcf.gz
```
We may also wish to add some custom stuff to the INFO field. For example, the hwe.c plugin included with this package will add the -log10(p-value) for Fisher's Exact test for divergence from Hardy-Weinberg Equilibrium, as well as the inbreeding coefficient.
```
bcftools view merged.flt.bcf -Ou | bcftools +hwe | bcftools view -G -Oz -o merged.sites.vcf.gz
tabix merged.sites.vcf.gz
```

