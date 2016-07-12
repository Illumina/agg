docker: FilesForDockerImage.tgz
	sudo docker build -t agg-part-1-local .

FilesForDockerImage.tgz: ../../agg ../../make_chunk.py ../curlWithRetry bsAppScript1.sh
	rm -rf agg
	mkdir agg
	cp $^ agg/
	tar czf FilesForDockerImage.tgz agg
	rm -rf agg
