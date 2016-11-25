wget --continue https://s3-eu-west-1.amazonaws.com/agg-examples/agg_test.tar.gz
tar -xzvf agg_test.tar.gz

for i in 1 2;do echo ../make_chunk.py gvcfs${i}.txt -output chunk${i} -ref ref.fa -tmp /tmp/ -nproc 4 -agg ../agg;done | xargs -P2 -l python

#for i in chrX chr21 chr22; do echo genotype -r $i chunk1.bcf chunk2.bcf -Ob -o ${i}.bcf -@ 4; done | xargs -P5 -l ../agg

../agg genotype -o - -r chr21 chunk1.bcf chunk2.bcf -Ou -@4 | ./bcftools filter -Ou -i 'FMT/DP>=10 & FMT/GQ>=20' -S. | ./bcftools annotate -Ob -x FILTER -o chr22.bcf


../agg anno chr21.bcf -i 'QUAL>=20'  -Ob -o chr21.anno.bcf

./vcftools --bcf chr21.bcf --diff-bcf omni.bcf  --diff-discordance-matrix

awk '{if(NR>1) {sum1+=$NR;sum2+=($2+$3+$4)} }END{print 100*(sum1/sum2)" % concordance"}' out.diff.discordance_matrix 
