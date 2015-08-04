#!/usr/bin/env bash

# Simple test of agg pipeline

# set this to the location of some test gvcf files
TESTDATA=/illumina/scratch/kimura/small.agg.test
DB=gvcf.small.db

PYTHON=/illumina/thirdparty/python/Python-2.7.3/bin/python

AGG=../..
REF=/illumina/development/Isis/Genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFASTA/genome.fa

GT_DEF_OUT=gt.chr20.5000000-6000000.default.vcf
CNT_DEF_OUT=cnt.chr20.5000000-6000000.default.vcf

# clean up
rm -rf ${DB}
#rm -f ingest_gvcf.q.*
#rm -f *canon.bcf
#rm -rf *tmp.vcf.gz

# initialize database
ls ${TESTDATA}/*.vcf.gz > gvcf.list
${PYTHON} ../initdb.py -db ${DB} -files gvcf.list -agg ${AGG} -ref ${REF} -chrom_prefix chr

# wait until sge jobs are done
RET=`qstat -u ${USER} | grep -c collate`
while [ ${RET} -ne 0 ] ; do
    echo "Waiting for initdb.py"
    sleep 20
    RET=`qstat -u ${USER} | grep -c collate`
done
echo "Done"

echo "Testing genotyping output"
GT_OUT=gt.chr20.5000000-6000000.vcf
rm -f ${GT_OUT}
../../agg genotype ${DB} -r chr20:5000000-6000000 -o ${GT_OUT}
if ! diff --brief ${GT_OUT} ${GT_DEF_OUT}; then
    echo "Difference in output files"
    #exit 1
else
    echo "All good"
fi

echo "Testing agg count"
CNT_OUT=cnt.chr20.5000000-6000000.vcf
rm -f ${CNT_OUT}
../../agg count ${DB} -r chr20:5000000-6000000 -o ${CNT_OUT}  
if ! diff --brief ${CNT_OUT} ${CNT_DEF_OUT}; then
    echo "Difference in output files"
    #exit 1
else
    echo "All good"
fi

exit $?
