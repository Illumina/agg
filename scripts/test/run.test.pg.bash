#!/usr/bin/env bash

# Simple test of agg pipeline

# set this to the location of some test gvcf files
TESTDATA=/illumina/build/platinumgenomes/Platinum_Genomes_Short_Insert/isaac_02.14/*/Data/Intensities/BaseCalls/Alignment/
DB=gvcf.pg.db

PYTHON=/illumina/thirdparty/Python-2.7.3/bin/python

AGG=../..
REF=/illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFASTA/genome.fa

GT_DEF_OUT=gt.chr20.5000000-6000000.default.vcf
CNT_DEF_OUT=cnt.chr20.5000000-6000000.default.vcf

# clean up
rm -rf ${DB}
rm -f ingest_gvcf.q.*
rm -f *canon.bcf
rm -rf *tmp.vcf.gz

# initialize database
ls ${TESTDATA}/*genome.vcf.gz > gvcf.list
${PYTHON} ../initdb.py -db ${DB} -files gvcf.list -agg ${AGG} -ref ${REF} -chrom_prefix chr

# wait until sge jobs are done
RET=`qstat -u ${USER} | grep -c concat`
while [ ${RET} -ne 0 ] ; do
    echo "Waiting for initdb.py"
    sleep 20
    RET=`qstat -u ${USER} | grep -c concat`
done
echo "Done"

#echo "Testing genotyping output"
#GT_OUT=gt.chr20.5000000-6000000.vcf
#rm -f ${GT_OUT}
#../../agg genotype ${DB} -r chr20:5000000-6000000 -o ${GT_OUT} > /dev/null 2>&1
#if ! diff --brief ${GT_OUT} ${GT_DEF_OUT}; then
#    echo "Difference in output files"
#else
#    echo "All good"
#fi

#echo "Testing agg count"
#CNT_OUT=cnt.chr20.5000000-6000000.vcf
#rm -f ${CNT_OUT}
#../../agg count gvcf.db -r chr20:5000000-6000000 -o ${CNT_OUT} > /dev/null 2>&1
#if ! diff --brief ${CNT_OUT} ${CNT_DEF_OUT}; then
#    echo "Difference in output files"
#else
#    echo "All good"
#fi

