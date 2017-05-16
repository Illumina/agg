CC=gcc
CXX=g++

all: agg plugins

CFLAGS = -O2  $(ALLFLAGS) 

debug: CFLAGS =  -O0 -pg -g -Wall $(ALLFLAGS) 
debug: all 

#linker stuff
HTSDIR = htslib-1.4
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR)
LFLAGS = -lz -lpthread -lm -ldl 
ALLFLAGS =  $(IFLAGS) 
CXX_FLAGS = -std=c++0x -Wno-write-strings


GIT_HASH := $(shell git describe --abbrev=4 --always )
BCFTOOLS_VERSION=1.4
VERSION = 0.3.5
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
#	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(LIBS)

plugins: $(PLUGINS)

##agg source code
utils.o: utils.cpp utils.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c utils.cpp 
agg_genotyper.o: agg_genotyper.cpp agg_genotyper.h version.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c agg_genotyper.cpp  
agg_anno.o: agg_anno.cpp agg_anno.h  version.h filter.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c agg_anno.cpp  
agg_utils.o: agg_utils.cpp agg.h 
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c $< 
depthMerger.o:  depthMerger.cpp depthMerger.h version.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c depthMerger.cpp 
agg_ingest2.o: agg_ingest2.cpp agg.h 
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c $<  
agg_ingest1.o: agg_ingest1.cpp agg_ingest1.h agg.h version.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c $<
needle.o: needle.cpp needle.h
	$(CXX) $(CFLAGS) $(CXX_FLAGS) -c $<
##these files were taken from bcftools hence compiled as straight C
filter.o: filter.c 
	$(CC) $(CFLAGS) -c $<  
version.o: version.c 
	$(CC) $(CFLAGS) -c $<  
vcfmerge.o: vcfmerge.c 
	$(CC) $(CFLAGS) -c $<  
vcfnorm.o: vcfnorm.c
	$(CC) $(CFLAGS) -c $<  
vcmp.o: vcmp.c
	$(CC) $(CFLAGS) -c $<
regidx.o: regidx.c
	$(CC) $(CFLAGS) -c $<  
##binary
agg: agg.cpp  agg_anno.o depthMerger.o vcfnorm.o vcmp.o vcfmerge.o  agg_ingest2.o utils.o agg_utils.o agg_genotyper.o  agg_ingest1.o version.o filter.o regidx.o needle.o version.h  $(HTSLIB)
	$(CXX) $(CFLAGS)  -o agg agg.cpp regidx.o vcfnorm.o vcmp.o vcfmerge.o agg_ingest2.o agg_genotyper.o utils.o agg_utils.o  agg_ingest1.o depthMerger.o version.o agg_anno.o filter.o needle.o $(HTSLIB) $(LFLAGS) 
test: agg
	cd test/;bash -e test.sh 

#housekeeping
all:  $(ALL)
clean:
	rm -rf *.o $(ALL) version.h bcftools_plugins/*.so
