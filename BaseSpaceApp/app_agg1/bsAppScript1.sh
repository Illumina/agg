#!/bin/bash -x
# Script v3 goes through all AppResults in a designated project to locate genome vcf files

echo "Starting BaseSpace App: $*"
echo "Project1=$Project1"
echo "OUTDIR=$OUTDIR"

export AGG_PATH=/agg/
export PATH=$PATH:${AGG_PATH}

env
df -h
#find /genomes | grep hg19
find /data/input

cd /data/scratch


mkdir headers
find /data/input/appresults -name '*.genome.vcf.gz' -print > chunk1

echo "CHUNK1:"
cat chunk1


# Extract references
#find /data/input/ -name '*.genome.vcf.gz' -exec zgrep -H '##reference=' $i; done 2>/dev/null > greppedReferences_vcf
##for i in headers/*vcf.gz ; do zgrep -H '##reference=' $i; done 2>/dev/null > greppedReferences_vcf_gz
#for i in headers/*vcf.gz.uncompressed; do  grep -H '##reference=' $i; done 2>/dev/null > greppedReferences_vcf_gz
#
#
## Show all possible reference names
#cat greppedReferences_vcf greppedReferences_vcf_gz | tr '\\' '/' | tr -d '\r' | sed 's|^.*/[gG]enomes/||' | sort -f | uniq -ic
#
##
#cat greppedReferences_vcf greppedReferences_vcf_gz | tr '\\' '/' | tr -d '\r' | grep -i "Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" | cut -d ':' -f 1 | sed 's|^headers/||' | sed 's|.uncompressed$||' > hg19.files
#sed 's/\([^_]*\)_/\1\t0\t/' hg19.files > hg19.files.tsv
#show total size: for i in `cat hg19.files`; do cat metadata/$i | jq .Response.Size; done | ~/Sum



# Run agg
mkdir chunks
time python ${AGG_PATH}/make_chunk.py chunk1 -o ${OUTDIR}/aggChunk -ref /genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -tmp `pwd` -nproc 32 | tee ${OUTDIR}/mylog


# Temporary, for debugging
cd /data/scratch
mkdir ${OUTDIR}/scratch
ls -laR > ${OUTDIR}/scratch/find.txt
ls -laR /data > ${OUTDIR}/scratch/find2.txt
find /data/input >  ${OUTDIR}/scratch/find3.txt
#mv * ${OUTDIR}/scratch/

