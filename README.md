agg: a utility for aggregating gvcfs

Copyright (c) 2015, Illumina, Inc. All rights reserved. See LICENSE.pdf for further details.

####Summary

This tool implements a simple pipeline to merge gvcfs in a dynamic fashion. That is, not all gvcfs need to merged at once, new groups of samples can be added periodically.  It achieves this by storing sample variants in a standard bcf and storing depth information in an auxilliary file (.dpt - depth track). The variants and depth files make up an "agg chunk" and these chunks can then be merged on the fly, with the depth file allowing the union of variants to be genotyped across all samples.  Note, here "genotyping" simply means we can gauge the coverage at a variant position for samples that do not have that variant, allowing us to examine the evidence that this sample was homozygous reference at that location.

If you are working with relatively few gvcfs, then [gvcftools](https://github.com/sequencing/gvcftools) is probably more appropriate.  The agg pipeline is more complex to run and really geared to 1000s of gvcfs.

####Installation
The only compilation dependency is [htslib](http://www.htslib.org/) which is included with the software.  

```
git clone git@github.com:BEETL/agg.git
cd agg/htslib-1.2.1/
make
cd ../
make
```

It should be noted parts of the agg source code were taken directly from the excellent [bcftools](https://github.com/samtools/bcftools) which is licensed permissively under BSD.

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

