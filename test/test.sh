wget --continue https://s3-eu-west-1.amazonaws.com/agg-examples/agg_test.tar.gz
tar -xzvf agg_test.tar.gz

for i in 1 2;do echo ../make_chunk.py gvcfs${i}.txt -output chunk${i} -ref ref.fa -tmp /tmp/ -nproc 24 -agg ../agg;done | xargs -P2 -l python

for i in chrX chr21 chr22; do echo genotype -r $i chunk1.bcf chunk2.bcf -Ob -o ${i}.bcf -@ 4; done | xargs -P5 -l ../agg


../agg anno chr21.bcf  -Ob -o chr21.anno.bcf


