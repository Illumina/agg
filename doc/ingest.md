#####Making an agg chunk
This describes how to build an agg chunk "manually", that is, not using the `make_chunk.py` script.

Creating a chunk is currently a two stage process using the `agg ingest1` and `agg ingest2` commands.  The individual gvcfs are first pre-processed with `ingest1` and then merged into a chunk with `agg ingest2`.  In the future, we plan to wrap these these two steps into a single (faster) `ingest` command.

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

In practice, we have 100s to 1000s of gvcfs and want to automate things. Say you have a 16 core node and want to build an agg chunk from 500 gvcfs listed in `gvcfs.txt`.  We can leverage the multiple CPUs using the `xargs` command.
```
$ for i in `cat gvcfs.txt`;
  do 
     out=$(basename ${i%.genome.vcf.gz});
     echo ingest1 $i -o ingest1/${out};
  done | xargs -P 16 -l agg
```
note you can replace `cat gvcfs.txt` with `find . -name '*.genome.vcf.gz'` or similar.

These files are the input for `agg ingest2`, which builds a chunk and is explained next. These intermediate files are rather large (1GB-2GB per sample), but after building a chunk they can be disposed of. You may wish to store them on local scratch space if it is available.

######ingest2: merge temporary files into a chunk
Let's roll these intermediate files into five separate agg chunks (100 samples per chunk). We will use `xargs` to leverage multiple cores again.
```
##get bcfs from ingest1 step
$ find ingest1/ -name '*.bcf' > ingest1.txt 
#splits the file into 5 groups of 100
$ split -d -l 100 ingest1.txt chunk_ 
$ ls chunk_*
chunk_00  chunk_01  chunk_02  chunk_03  chunk_04
$ for i in chunk_*;do echo ingest2 -l $i -@4 -o $i;done | xargs -l -P 16 agg
...
$ ls chunk_*.*
```
