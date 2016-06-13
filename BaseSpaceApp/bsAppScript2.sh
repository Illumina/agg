#!/bin/bash -x
# Agg app 2
#	 //texlive-latex-base
echo "Starting BaseSpace App: $*"

InputAppResultDir=$1
echo "InputAppResultDir=${InputAppResultDir}"
echo "OUTDIR=$OUTDIR"

env

#ServerUri=`cat /data/input/AppSession.json | jq --raw-output '.OriginatingUri' | sed 's/^https:\/\///'`

cd /data/scratch
export PATH=/agg2:$PATH
AKT_DATA=/agg2/data
AKT_SCRIPTS=/agg2/scripts
S3_DATA=/data/scratch/s3


# setup
#export PATH=$PATH:/illumina/thirdparty/bcftools/bcftools-1.2/:~/workspace/git/agg:~/workspace/git/akt
#mkdir /tmp/agg
#ln -s /illumina/build/1000GenomesReference/20130502/${S3_DATA}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz .
#ln -s /illumina/build/1000GenomesReference/20130502/${S3_DATA}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi .
#ln -s /illumina/build/1000GenomesReference/20130502/cgi_variant_calls/${S3_DATA}/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf .
#ln -s /illumina/build/1000GenomesReference/20130502/cgi_variant_calls/${S3_DATA}/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf.csi .

#cd "/basespace/ljanin./hoth_basespaceuser2/Projects/AutomaticallyGeneratedAppTest/AppResults/kmers (2)/Files"
# execute the following in parallel with the db download

# Discover chunks
CHUNKS=`find /data/input/appresults -name aggChunk.bcf -print`
echo CHUNKS=${CHUNKS}

for i in `seq 1 22` X Y ; do
  time agg genotype -r chr${i} ${CHUNKS} -Ob -o merged.chr${i}.bcf & # 1min
done


# Download db files
mkdir -p ${S3_DATA}
for i in \
  ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf.csi \
  ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf \
  ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz \
  ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz.tbi \
    ; do time curl https://s3.amazonaws.com/illumina-ukch-compbio/1000GenomesReference/20130502/$i -o ${S3_DATA}/$i & done
wait


for i in `seq 1 22` X Y ; do
 (
  input=merged.chr${i}.bcf 
  bcftools index ${input}

  # basic summary of number of sites, intersection with 1000G, etc.
  time bcftools stats $input ${S3_DATA}/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz > stats_chr${i}.txt # 4.5min
  plot-vcfstats stats_chr${i}.txt -p stats_chr${i}_dir/ # fails quickly
 ) &
done
wait

#Rudy's tool for PCA
#NO_OR_WITH=no
#time akt pca -w ${AKT_DATA}/1000G.snps.${NO_OR_WITH}chr.vcf.gz  ${S3_DATA}/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf > pca.txt # 1min
#Rscript ${AKT_SCRIPTS}/1000G_pca.R pca.txt > 1000G_pca.txt # null device



# 
##akt ibd ${S3_DATA}/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -T ${AKT_DATA}/1000G.snps.nochr.vcf.gz -n 32 > kin.txt
#time akt kin ${S3_DATA}/ALL.cgi_multi_sample.20130725.snps_indels.high_coverage_cgi.normalized.uniq.genotypes.gtonly.cr90.ic10.bcf -T ${AKT_DATA}/1000G.snps.nochr.vcf.gz -n 32 > kin.txt # 1.5min
#akt relatives kin.txt > fam.txt
#(
#grep '^Dup' fam.txt | awk '{print $1}' | uniq | awk 'END {printf "%s Duplicates\n", NR}'
#grep '^Unrel' fam.txt | awk '{print $1}' | uniq | awk 'END {printf "%s Unrelated\n", NR}'
#grep '^Fam' fam.txt | awk '{print $1}' | uniq | awk 'END {printf "%s Families\n", NR}'
#grep 'Parent' fam.txt | awk 'END {printf "%s Parent/Child pairs\n", NR}'
#grep 'Sibling' fam.txt | awk 'END {printf "%s Sibling pairs\n", NR}'
#grep 'Second-order' fam.txt | awk 'END {printf "%s Second-order pairs\n", NR}'
#) > fam_out.txt


ls -lartR
rm -rf ${S3_DATA}

#find . -maxdepth 1 -type f -size '-1000k' -print -exec mv {} ${OUTDIR}/ \;
mv * ${OUTDIR}/

ls -lartR

