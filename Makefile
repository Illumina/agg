CC=gcc
CXX=g++

all: agg test

CFLAGS = -O2  $(ALLFLAGS) 

debug: CFLAGS =  -pg -g -Wall $(ALLFLAGS) 
debug: all 

#linker stuff
HTSDIR = htslib-1.3
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR)
LFLAGS = -lz -lpthread -lm -ldl 
ALLFLAGS =  $(IFLAGS) 
CXX_FLAGS = -std=c++0x -Wno-write-strings


GIT_HASH := $(shell git describe --abbrev=4 --always )
BCFTOOLS_VERSION=1.3
VERSION = 0.3.4
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
version.h:
	echo '#define VERSION "$(VERSION)"' > $@
	echo '#define BCFTOOLS_VERSION "$(BCFTOOLS_VERSION)"' >> $@

##agg source code
utils.o: utils.cpp utils.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c utils.cpp 
agg_genotyper.o: agg_genotyper.cpp agg_genotyper.h version.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c agg_genotyper.cpp  
agg_utils.o: agg_utils.cpp agg.h 
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c $< 
depthMerger.o:  depthMerger.cpp depthMerger.h version.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c depthMerger.cpp 
agg_ingest2.o: agg_ingest2.cpp agg.h 
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c $<  
agg_ingest1.o: agg_ingest1.cpp agg_ingest1.h agg.h version.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c $<  
##these files were taken from bcftools hence compiled as straight C
version.o: version.c 
	$(CC) $(CFLAGS) -c $<  
vcfmerge.o: vcfmerge.c 
	$(CC) $(CFLAGS) -c $<  
vcfnorm.o: vcfnorm.c
	$(CC) $(CFLAGS) -c $<  
vcmp.o: vcmp.c
	$(CC) $(CFLAGS) -c $<  
##binary
agg: agg.cpp depthMerger.o vcfnorm.o vcmp.o vcfmerge.o  agg_ingest2.o utils.o agg_utils.o agg_genotyper.o  agg_ingest1.o version.o version.h  $(HTSLIB)
	$(CXX) $(CFLAGS)  -o agg agg.cpp vcfnorm.o vcmp.o vcfmerge.o agg_ingest2.o agg_genotyper.o utils.o agg_utils.o  agg_ingest1.o depthMerger.o version.o  $(HTSLIB) $(LFLAGS) 
agg_ingest: agg_ingest.cpp vcfnorm.o vcmp.o vcfmerge.o  agg_ingest2.o utils.o agg_utils.o agg_genotyper.o  agg_ingest1.o version.h version.o $(HTSLIB)
	$(CXX) $(CFLAGS)  -o agg agg.cpp vcfnorm.o vcmp.o vcfmerge.o agg_ingest2.o agg_genotyper.o utils.o agg_utils.o  agg_ingest1.o $(HTSLIB) $(LFLAGS) 
test:	agg
	echo > test

#housekeeping
all:  $(ALL)
clean:
	rm -rf *.o $(ALL) version.h

