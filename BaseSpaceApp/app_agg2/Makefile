docker.remote.done: app_agg2.json docker.local
	~/automatedNativeAppCreation/bs-create-native-app $<

docker.local: FilesForDockerImage.tgz
	sudo docker build -t agg-part-2-local .

FilesForDockerImage.tgz: ../../agg bsAppScript2.sh ../../../akt/akt ../../../akt/scripts ../../../akt/data /illumina/thirdparty/bcftools/bcftools-1.2/bcftools /illumina/thirdparty/bcftools/bcftools-1.2/plot-vcfstats
	rm -rf agg2
	mkdir agg2
	cp -r $^ agg2/
	tar czf FilesForDockerImage.tgz agg2
	rm -rf agg2
