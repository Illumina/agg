docker.remote.done: app_agg1.json docker.local
	~/automatedNativeAppCreation/bs-create-native-app $<

docker.local: FilesForDockerImage.tgz
	sudo docker build -t agg-part-1-local .

FilesForDockerImage.tgz: ../../agg ../../make_chunk.py ../curlWithRetry bsAppScript1.sh
	rm -rf agg
	mkdir agg
	cp $^ agg/
	tar czf FilesForDockerImage.tgz agg
	rm -rf agg
