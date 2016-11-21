wget --continue https://s3-eu-west-1.amazonaws.com/agg-examples/agg_test.tar.gz
tar -xzvf agg_test.tar.gz

for i in 1 2;do echo ../make_chunk.py gvcfs${i}.txt -output chunk${i} -ref ref.fa -tmp /tmp/ -nproc 24 -agg ../agg;done | xargs -P2 -l python

for i in chrX chr21 chr22; do echo genotype -r $i chunk1.bcf chunk2.bcf -Ob -o ${i}.bcf -@ 4; done | xargs -P5 -l ../agg


bcftools index chr21.bcf
bcftools index chr22.bcf

bcftools concat chr21.bcf chr22.bcf -Ou  | bcftools view -Ob -o autosomes.bcf --threads 8
bcftools index autosomes.bcf

../agg anno autosomes.bcf  -Ob -o autosomes.anno.bcf

