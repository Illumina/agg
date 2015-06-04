agg2: a utility for aggregating gvcfs

Copyright (c) 2015, Illumina, Inc. All rights reserved. See LICENSE.pdf for further details.

####Installation

The only compilation dependency is [htslib](http://www.htslib.org/) which is included with the software.  These tools are needed during ingestion:

* [gvcftools 0.16](https://sites.google.com/site/gvcftools/)
* [bcftools 1.2](http://samtools.github.io/bcftools/bcftools.html)

the Makefile attempts to download and compile both of these.

```
cd kimura/agg/
make
```

####Usage

```

$ ./agg

Program:        agg 059d9a7 (aggregation tool for multiple samples)
Contact:        joconnell@illumina.com

Copyright (c) 2015, Illumina, Inc. All rights reserved. See LICENSE.pdf for further details.

Usage:  agg <command> [options]

Commands:
        collate         initialises the database's variant list
        update          adds new samples to the database
        count           summary statistics (genotype counts, passrate, etc)
        genotype        produce multisample vcf from the database
```

####Initial ingestion

Assume gvcfs.txt contains a list of all our gvcfs we wish to initialise our database with.  The database is stored in a directory database/

Steps:

1. "Ingest" our data into separate files of homref beds and variant BCFs
2. Crawl our data to build a list of all variants present in the samples

This is all peformed by the `initdb.py` script and utilises either Sun Grid Engine (qsub) or Slurm (sbatch)

```
#qsubs our ingestion jobs. gvcfs.txt is simply a list of gvcfs. agg is just the agg install path.
python agg/scripts/initdb.py -db database/ -files gvcfs.txt 
       			     	       \ -agg ~/agg/
				       \ -ref /illumina/scratch/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
				       \ -chrom_prefix chr
				       \ -nprocess 8
```

Wait until all cluster jobs finished and then the database is ready to query.

####Updating the databse

The whole point of agg is that it allows batches of samples to be added to the database (otherwise why not just use a multisample vcf). Adding samples to an existing data is similar to initialising a database, we use the script `updatedb.py`.  This script can also be used to check/correct an existing database in the case where some ingestion jobs failed.
```
python /home/joconnell/kimura/agg/scripts/updatedb.py -db database/ -agg /home/joconnell/kimura/agg/ -files gvcf.list -ref /illumina/scratch/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
```
This assumes `database` already exists.  If gvcf.list contains samples that are either:
1. Not in the database
2. In the database, but incorrectly ingested (invalid bed or bcf)

these samples will be (re-)ingested.  If a sample is already correctly in the database, it will not be re-ingested.  That is, you can have redundant samples in gvcf.list and not worry about wasting computation re-ingesting them.

####Genotyping

```
$ ./agg genotype

About:   genotype samples from <database> at a region in the genome.
Usage:   agg genotype <database>

Required options:
    -r, --regions <region>              region to genotype eg. chr1 or chr20:5000000-6000000

Output options:
    -o,   --output-file <file>          output file name [stdout]
    -O,   --output-type <b|u|z|v>       b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]

Subset options:
    -s, --samples [^]<list>       comma separated list of samples to include (or exclude with "^" prefix)
    -S, --samples-file [^]<file>  file of samples to include (or exclude with "^" prefix)

Filter options:
    TBA
```

Example:

```
[joconnell@ussd-prd-lndt-b-8-1 test]$ ~/agg/agg genotype database/ -r chr20 -O b -o test.chr20.bcf
Using database database/
Writing output to test.chr20.bcf
Genotyping region chr20
221798 variants to genotype...
Including all 10 samples in database.
Writing output...
0 LP6005600-DNA_H12
1 LP6005599-DNA_C04
2 LP6005600-DNA_C03
3 LP6005600-DNA_B08
4 LP6005600-DNA_D10
5 LP6005600-DNA_A09
6 LP6005600-DNA_A12
7 LP6005600-DNA_A08
8 LP6005616-DNA_C01
9 LP6005616-DNA_F01

$ bcftools view -H test.chr20.bcf | head
chr20   61098   .       C       T       176     .       AN=20;AC=4;NPASS=4      GT      0/0     0/1     0/0     0/1     0/0     0/0     0/1     0/0     0/1     0/0
chr20   61138   .       C       CT      151     .       AN=20;AC=2;NPASS=0      GT      0/0     0/1     0/0     0/0     0/0     0/0     0/0     0/0     0/1     0/0
chr20   61270   .       A       C       62      .       AN=14;AC=2;NPASS=1      GT      0/0     ./.     ./.     0/1     0/0     0/0     0/1     0/0     0/0     ./.
chr20   61795   .       G       T       208     .       AN=20;AC=6;NPASS=6      GT      0/1     0/1     0/1     0/1     0/0     0/0     0/1     0/0     0/1     0/0
chr20   62731   .       C       A       179     .       AN=20;AC=2;NPASS=2      GT      0/1     0/0     0/1     0/0     0/0     0/0     0/0     0/0     0/0     0/0
chr20   62975   .       T       TA      133     .       AN=20;AC=1;NPASS=0      GT      0/0     0/0     0/0     0/0     0/1     0/0     0/0     0/0     0/0     0/0
chr20   63244   .       A       C       189     .       AN=20;AC=3;NPASS=3      GT      0/0     0/1     0/0     0/0     0/0     0/0     0/1     0/0     0/1     0/0
chr20   63452   .       C       G       195     .       AN=20;AC=1;NPASS=1      GT      0/0     0/0     0/0     0/1     0/0     0/0     0/0     0/0     0/0     0/0
chr20   63799   .       C       T       264     .       AN=20;AC=6;NPASS=6      GT      0/1     0/1     0/1     0/1     0/0     0/0     0/1     0/0     0/1     0/0
chr20   63967   .       A       G       244     .       AN=20;AC=1;NPASS=1      GT      0/0     0/0     0/0     0/1     0/0     0/0     0/0     0/0     0/0     0/0
```

You might notice this process is rather slow.  We are currently working on better performing data structures.  In the interim, you can simply perform queries in genomics chunks and concatenate them. We provide a tool to specify sensibly size chunks.  See below:

```
$ ~/agg/chunker -n 100000 testdb2/sites.bcf > regions.txt
$ for r in `cat regions.txt`;do echo genotype testdb2/ -r ${r} -O b -o ${r}.bcf ;done | xargs -P 8 -n 8 ~/agg/agg
$ for r in `cat regions.txt`;do echo ${r}.bcf;done > file_list.txt
$ bcftools concat -f file_list.txt -O b -o final.bcf
```
The above example generates a list of windows with ~100,000 variants in them using the `chunker` tool, with a large database you would want to set this to perhaps ~5,000,000.  We then run one `agg genotype` process per region, we use the unix `xargs` command to run eight jobs in parallel. Finally, we concatenate the data with `bcftools`.

####Counting genotypes
Rather than creating multisample vcfs, you may prefer just to calculate some summary statistics. 

```
 ./agg count ~/scratch/agg/test -r chr20:10000000-10010000 -O z -o ~/test.vcf.gz

$ bcftools view -H ~/test.vcf.gz  | head
chr20	10000117	.	C	T	125	.	AC=2;AN=4;AF=0.5;PASSRATE=1,0;CALLRATE=1
chr20	10000211	.	C	T	199	.	AC=2;AN=4;AF=0.5;PASSRATE=1,0;CALLRATE=1
chr20	10000439	.	T	G	499	.	AC=3;AN=4;AF=0.75;PASSRATE=1,0.5;CALLRATE=1
chr20	10000598	.	T	A	491	.	AC=3;AN=4;AF=0.75;PASSRATE=1,0.5;CALLRATE=1
chr20	10000694	.	G	A	177	.	AC=2;AN=4;AF=0.5;PASSRATE=1,0;CALLRATE=1
chr20	10000758	.	T	A	521	.	AC=3;AN=4;AF=0.75;PASSRATE=1,0.5;CALLRATE=1
chr20	10001019	.	T	G	122	.	AC=2;AN=4;AF=0.5;PASSRATE=1,0;CALLRATE=1
chr20	10001298	.	T	A	470	.	AC=3;AN=4;AF=0.75;PASSRATE=1,0.5;CALLRATE=1
chr20	10001436	.	A	AAGGCT	1154	.	AC=3;AN=4;AF=0.75;PASSRATE=1,0.5;CALLRATE=1
chr20	10001474	.	C	T	456	.	AC=3;AN=4;AF=0.75;PASSRATE=1,0.5;CALLRATE=1

```
Again, this can be sped up by chunking and concatenating using the same procedure described in the genotyping section.
