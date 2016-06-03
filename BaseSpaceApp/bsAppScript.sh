#!/bin/bash -x
# Script v3 goes through all AppResults in a designated project to locate genome vcf files

echo "Starting BaseSpace App: $*"
echo "Project1=$Project1"
echo "OUTDIR=$OUTDIR"

export PATH=$PATH:/agg

env
df -h
#find /genomes | grep hg19

ServerUri=`cat /data/input/AppSession.json | jq --raw-output '.OriginatingUri' | sed 's/^https:\/\///'`

cd /data/scratch

# Get better access token
#appId=`cat /data/input/AppSession.json |jq '.Application.Id | tonumber'`
#AccessToken=`curl -H x-access-token:${AccessToken} https://api.cloud-hoth.illumina.com/v1pre3/oauthv2/privatetoken?ApplicationId=${appId} | jq --raw-output .Response.AccessToken`


#optional input project
#compulsory destination project

#set default min timestamp = 0
#Go through destination project's appresults (it has to be all of them, we can't restrict it to the current app id, as versioning would change the app id), and try to detect the last one that contains an aggregate result (with its cutoff timestamp and id of its optional input project, which will override the defaults form above)
outputProjectId=`cat /data/input/AppSession.json | jq --raw-output '.Properties.Items[] | select(.Name=="Input.project-id") | .Content.Id'`
#set default input project
inputProjectId=`cat /data/input/AppSession.json | jq --raw-output '.Properties.Items[] | select(.Name=="Input.Project1") | .Content.Id'`

#todo loop until find what we want
#curl -H x-access-token:${AccessToken} "https://api.${ServerUri}/v1pre3/projects/${inputProjectId}/appresults?SortBy=DateCreated&SortDir=desc&limit=1000"


#query all appsessions
#	(not appresults, as we wouldn't be able to do a whole account query, and doing multiple calls would lead to "race conditions" when figuring out the cutoff timestamp to save)
#	optionally with project filter
#	only Completed status
#	sorted by [creation or modification?] date
#	remember total number of results for this page 1
#ask for next pages until we reach min timestamp
#Do the same query again for 1 item, and check that the total number of results is still the same as before. If it's not, it means that some appsessions have completed just in the last few seconds => we should restart the whole query just to make sure our cutoff point is precise enough so that we don't miss anything next time we run the tool.
userId=`cat /data/input/AppSession.json | jq -r .UserCreatedBy.Id`
pageSize=1024
offset=0
echo > fileInfos
#time curl -H x-access-token:${AccessToken} "https://api.${ServerUri}/v1pre3/users/current/appsessions?limit=1024&sortby=DateCreated&sortdir=desc&userCreatedBy=${userId}&status=Complete" > appSessions.json
while [ 1 ] ; do
time curl --silent -H x-access-token:${AccessToken} "https://api.${ServerUri}/v1pre3/users/current/appsessions?output.projects=${inputProjectId}&offset=${offset}&limit=${pageSize}&sortby=DateCreated&sortdir=desc&userCreatedBy=${userId}&status=Complete&include=properties&propertyFilters=Output.AppResults" > appSessions.json
if [ ${offset} == 0 ] ; then cp appSessions.json appSessions0.json ; fi
#cat appSessions.json | jq -r .Response.Items[].Id | sort -g | uniq > appSessionIds
count=`cat appSessions.json | jq -r .Response.DisplayedCount`
if [ "${count}" == "null" ] ; then exit 1; fi
if [ "x${count}" == "x" ] ; then exit 2; fi
if [ ${count} == 0 ] ; then break; fi
cat appSessions.json | jq -r .Response.Items[].Properties.Items[0].Items[].Id > appResultIds 2> /dev/null


#save first app sessions's timestamp (+1?) as cutoff date for next run
#go through all appsessions' appresults, and make a list of all vcf file ids
for id in `cat appResultIds`; do 
  curl --silent -H x-access-token:${AccessToken} "https://api.${ServerUri}/v1pre3/appresults/${id}/files?limit=1024" > files
  cat files | jq -r '.Response.Items[] | select(.Name | contains("genome.vcf")) | "\(.Id)\t\(.Size)\t\(.Name)"' >> fileInfos
done

let offset=offset+pageSize
done



# Filter
grep -P ".genome.vcf$" fileInfos | sort -gk2 > fileInfos.genome_vcf
grep -P ".genome.vcf.gz$" fileInfos  | sort -gk2 > fileInfos.genome_vcf_gz



#if 0 vcf: ?

#add list of vcf to report

#generate restricted list of files that we want to use
cat fileInfos.genome_vcf > fileInfos.genome_vcf.restricted
cat fileInfos.genome_vcf_gz > fileInfos.genome_vcf_gz.restricted

#generate Makefile to download metadata
mkdir metadata
echo -n "all: " > metadata/Makefile
cut -f 1,3 fileInfos.genome_vcf.restricted fileInfos.genome_vcf_gz.restricted | tr '\t \n' '__ ' >> metadata/Makefile
echo >> metadata/Makefile
echo >> metadata/Makefile
cat fileInfos.genome_vcf.restricted fileInfos.genome_vcf_gz.restricted | awk -v AccessToken=${AccessToken} -v ServerUri=${ServerUri} '{ print $1 "_" $3 ":\n\tcurl --silent -o $@ -H x-access-token:" AccessToken " \"https://api." ServerUri "/v1pre3/files/" $1 "&FileHrefContentResolution=true\"\n" }' >> metadata/Makefile

#download all metadata
cd metadata
make -j 32
cd ..


#generate Makefile to download headers (and decompress the .vcf.gz ones while ignoring the "file too short" errors)
mkdir headers
echo -n "all: " > headers/Makefile
cut -f 1,3 fileInfos.genome_vcf.restricted fileInfos.genome_vcf_gz.restricted | tr '\t \n' '__ ' >> headers/Makefile
echo >> headers/Makefile
echo >> headers/Makefile
cat fileInfos.genome_vcf.restricted | awk -v AccessToken=${AccessToken} -v ServerUri=${ServerUri} '{ print $1 "_" $3 ":\n\tcurlWithRetry --silent --retry 3 -o $@ -LH x-access-token:" AccessToken " \"https://api." ServerUri "/v1pre3/files/" $1 "/content\" --range 0-8192\n" }' >> headers/Makefile
cat fileInfos.genome_vcf_gz.restricted | awk -v AccessToken=${AccessToken} -v ServerUri=${ServerUri} '{ print $1 "_" $3 ":\n\tcurlWithRetry --silent --retry 3 -o $@ -LH x-access-token:" AccessToken " \"https://api." ServerUri "/v1pre3/files/" $1 "/content\" --range 0-8192 && ( zcat $@ > $@.uncompressed 2> /dev/null ; true )\n" }' >> headers/Makefile

#download all headers
cd headers
make -j 16
cd ..


# Extract references
for i in headers/*vcf    ; do  grep -H '##reference=' $i; done 2>/dev/null > greppedReferences_vcf
#for i in headers/*vcf.gz ; do zgrep -H '##reference=' $i; done 2>/dev/null > greppedReferences_vcf_gz
for i in headers/*vcf.gz.uncompressed; do  grep -H '##reference=' $i; done 2>/dev/null > greppedReferences_vcf_gz


# Show all possible reference names
cat greppedReferences_vcf greppedReferences_vcf_gz | tr '\\' '/' | tr -d '\r' | sed 's|^.*/[gG]enomes/||' | sort -f | uniq -ic

#
cat greppedReferences_vcf greppedReferences_vcf_gz | tr '\\' '/' | tr -d '\r' | grep -i "Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa" | cut -d ':' -f 1 | sed 's|^headers/||' | sed 's|.uncompressed$||' > hg19.files
sed 's/\([^_]*\)_/\1\t0\t/' hg19.files > hg19.files.tsv
#show total size: for i in `cat hg19.files`; do cat metadata/$i | jq .Response.Size; done | ~/Sum



#generate Makefile to download files
LIST_OF_FILES_TO_DOWNLOAD=hg19.files.tsv
#LIST_OF_FILES_TO_DOWNLOAD=fileInfos.genome_vcf.restricted fileInfos.genome_vcf_gz.restricted

mkdir downloadedFiles
echo -n "all: " > downloadedFiles/Makefile
cut -f 1,3 ${LIST_OF_FILES_TO_DOWNLOAD} | tr '\t \n' '__ ' >> downloadedFiles/Makefile
echo >> downloadedFiles/Makefile
echo >> downloadedFiles/Makefile
cat ${LIST_OF_FILES_TO_DOWNLOAD} | awk -v AccessToken=${AccessToken} -v ServerUri=${ServerUri} '{ print $1 "_" $3 ":\n\tcurlWithRetry --silent --retry 3 -o $@ -LH x-access-token:" AccessToken " \"https://api." ServerUri "/v1pre3/files/" $1 "/content\"\n" }' >> downloadedFiles/Makefile


#download all vcfs
cd downloadedFiles
make -j 16
#for i in `cut -f 1,3 fileInfos.genome_vcf fileInfos.genome_vcf_gz | tr '\t ' '__'`; do id=${i%_*}; name=${i#*_}; echo "$i $id $name" ; curl -o downloadedFiles/${i} -LH x-access-token:${AccessToken} "https://api.${ServerUri}/v1pre3/files/${id}/content" ; done

# decompress vcf.gz (needed?)
#gunzip *vcf.gz

# Run agg
ls *vcf *.vcf.gz > chunk1
mkdir ../chunks
time python ~/agg/make_chunk.py chunk1 -o ${OUTDIR}/aggChunk -ref /genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -tmp `pwd` --allow-duplicates -nproc 32 | tee ${OUTDIR}/mylog


# Temporary, for debugging
#cd /data/scratch
#mkdir ${OUTDIR}/scratch
#mv * ${OUTDIR}/scratch/

