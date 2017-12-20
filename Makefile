CC=gcc
CXX=g++

all: agg

CFLAGS = -O2  $(ALLFLAGS) 
CXXFLAGS = -std=c++11 -Wno-write-strings -O2

debug: CFLAGS =  -O0 -pg -g -Wall $(IFLAGS)
debug: all 

#linker stuff
HTSDIR = htslib-1.5
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR) -I./
LFLAGS = -lz -lpthread -lm -ldl 

GIT_HASH := $(shell git describe --abbrev=4 --always )
BCFTOOLS_VERSION=1.4
VERSION = 0.3.7
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
version.h:
	echo '#define VERSION "$(VERSION)"' > $@
	echo '#define BCFTOOLS_VERSION "$(BCFTOOLS_VERSION)"' >> $@

# Plugin rules
PLUGINC = $(foreach dir, bcftools_plugins, $(wildcard $(dir)/*.c))
PLUGINS = $(PLUGINC:.c=.so)

%.so: %.c version.h version.c
	$(CC) -fPIC -shared -g -Wall -Wc++-compat -O2 -I. -I$(HTSDIR) -o $@ version.c $<

plugins: $(PLUGINS)

OBJS=agg_anno.o depthMerger.o vcfnorm.o vcmp.o vcfmerge.o  agg_ingest2.o utils.o agg_utils.o agg_genotyper.o  agg_ingest1.o version.o filter.o regidx.o needle.o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o $@ $<
.c.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c -o $@ $<

##agg source code
utils.o: utils.cpp utils.h
agg_genotyper.o: agg_genotyper.cpp agg_genotyper.h version.h
agg_anno.o: agg_anno.cpp agg_anno.h  version.h filter.h
agg_utils.o: agg_utils.cpp agg.h 
depthMerger.o:  depthMerger.cpp depthMerger.h version.h
agg_ingest2.o: agg_ingest2.cpp agg.h 
agg_ingest1.o: agg_ingest1.cpp agg_ingest1.h agg.h version.h
needle.o: needle.cpp needle.h
##these files were taken from bcftools hence compiled as straight C
filter.o: filter.c 
version.o: version.c 
vcfmerge.o: vcfmerge.c 
vcfnorm.o: vcfnorm.c
vcmp.o: vcmp.c
regidx.o: regidx.c

##binary
agg: agg.cpp $(OBJS) version.h  $(HTSLIB)
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o agg agg.cpp $(OBJS) $(HTSLIB) $(LFLAGS) 
test: agg
	cd test/;bash -e test.sh 

#housekeeping
all:  $(ALL)
clean:
	rm -rf *.o $(ALL) version.h bcftools_plugins/*.so
