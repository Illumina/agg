CC = gcc
CXX = g++

HTSLIB = htslib-1.2.1/
IFLAGS = -I$(HTSLIB)
LFLAGS = -lz -lpthread -lm -ldl 

//CFLAGS = -pg -g -Wall $(IFLAGS)  -Wno-write-strings 
CFLAGS = -O3  $(IFLAGS)  -Wno-write-strings  

release: CFLAGS = -O3  $(IFLAGS)   -Wno-write-strings  -fopenmp
release: all

static: CFLAGS = -O3 -static $(IFLAGS)  -Wno-write-strings
static: all

debug: CFLAGS = -pg -g -Wall  $(IFLAGS) -Wno-write-strings
debug: all 

ALL=bedmerge chunker canon agg version.h bcftools-1.2 gvcftools-0.16 vt-0.57/vt $(HTSLIB) $(HTSLIB)/libhts.a htslib-1.2.1/bgzip htslib-1.2.1/tabix # sqltest

$(HTSLIB) $(HTSLIB)/libhts.a htslib-1.2.1/bgzip htslib-1.2.1/tabix:
	tar -xjf htslib-1.2.1.tar.bz2 && \
	cd htslib-1.2.1/ && \
	make all
version.h: 
	git log --pretty=format:'#define VERSION "%h"' -n 1 > version.h
	echo >> version.h
sqlite3/sqlite3.o:
	$(CXX) $(CFLAGS) sqlite3/sqlite3.c -c -o sqlite3/sqlite3.o
utils.o: utils.cpp utils.h
	$(CXX) $(CFLAGS) -c utils.cpp 
#sqltest: sqlite3/sqlite3.o utils.o aggutils.cpp agg.h  sqltest.cpp
#	$(CXX) $(CFLAGS)  -o $@ sqltest.cpp aggutils.cpp utils.o sqlite3/sqlite3.o $(HTSLIB) $(LFLAGS) -Isqlite3/
sampleReader.o: sampleReader.cpp sampleReader.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c sampleReader.cpp  
sampleMerger.o: sampleMerger.cpp sampleMerger.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c sampleMerger.cpp  
canon: canon.cpp  utils.o $(HTSLIB)/libhts.a
	$(CXX) $(CFLAGS) utils.o canon.cpp -o canon $(HTSLIB)/libhts.a $(LFLAGS)
tabReader.o: tabReader.cpp tabReader.h
	$(CXX) $(CFLAGS) -c tabReader.cpp 
bedmerge: bedmerge.cpp tabReader.o utils.o  $(HTSLIB)/libhts.a
	$(CXX) $(CFLAGS) tabReader.o utils.o bedmerge.cpp -o bedmerge $(HTSLIB)/libhts.a $(LFLAGS)
chunker: chunker.cpp $(HTSLIB)/libhts.a
	$(CXX) $(CFLAGS) chunker.cpp -o chunker $(HTSLIB)/libhts.a $(LFLAGS)
agg_utils.o: agg_utils.cpp agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $< 
agg_collate.o: agg_collate.cpp agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
agg_genotype.o: agg_genotype.cpp agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
agg_count.o: agg_count.cpp agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
agg_update.o: agg_update.cpp agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
agg: agg.cpp sampleReader.o sampleMerger.o utils.o agg_utils.o agg_collate.o  agg_genotype.o  agg_update.o agg_count.o version.h $(HTSLIB)/libhts.a
	$(CXX) $(CFLAGS)  -o agg agg.cpp sampleReader.o sampleMerger.o utils.o agg_utils.o agg.h agg_collate.o agg_genotype.o agg_update.o $(HTSLIB)/libhts.a agg_count.o $(LFLAGS) 
mergeBlocks: mergeBlocks.cpp
	$(CXX) $(CFLAGS) mergeBlocks.cpp -o mergeBlocks
gvcftools-0.16:
	bash get_gvcftools.sh && \
	$(MAKE) -C "$@" 
bcftools-1.2:
	tar -xjf bcftools-1.2.tar.bz2 && \
	cd bcftools-1.2/ && \
	make
vt-0.57/vt:
	tar -xzvf vt-0.57.tar.gz && \
	cd vt-0.57 && \
	make
all:  $(ALL)
clean:
	rm -rf *.o $(ALL)
