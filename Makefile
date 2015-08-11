all: agg

CFLAGS = -O2  $(ALLFLAGS) 

debug: CFLAGS =  -pg -g -Wall $(ALLFLAGS) 
debug: all 

#linker stuff
HTSDIR = htslib-1.2.1
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR)
LFLAGS = -lz -lpthread -lm -ldl 
ALLFLAGS = -std=c++0x  -Wno-write-strings  $(IFLAGS) 

##objects libs etc
version.h: 
	git log --pretty=format:'#define VERSION "%h"' -n 1 > version.h
	echo >> version.h
utils.o: utils.cpp utils.h
	$(CXX) $(CFLAGS) -c utils.cpp 
agg_genotyper.o: agg_genotyper.cpp agg_genotyper.h
	$(CXX) $(CFLAGS) -c agg_genotyper.cpp  
agg_utils.o: agg_utils.cpp agg.h 
	$(CXX) $(CFLAGS) -c $< 
agg_ingest2.o: agg_ingest2.cpp agg.h 
	$(CXX) $(CFLAGS) -c $<  
vcfmerge.o: vcfmerge.c 
	$(CXX) $(CFLAGS) -c $<  
vcfnorm.o: vcfnorm.c
	$(CXX) $(CFLAGS) -c $<  
vcmp.o: vcmp.c
	$(CXX) $(CFLAGS) -c $<  
agg_ingest1.o: agg_ingest1.cpp agg_ingest1.h agg.h
	$(CXX) $(CFLAGS) -c $<  

##binaries
agg: agg.cpp vcfnorm.o vcmp.o vcfmerge.o agg_ingest2.o utils.o agg_utils.o agg_genotyper.o  agg_ingest1.o version.h $(HTSLIB)
	$(CXX) $(CFLAGS)  -o agg agg.cpp vcfnorm.o vcmp.o vcfmerge.o agg_ingest2.o agg_genotyper.o utils.o agg_utils.o  agg_ingest1.o $(HTSLIB) $(LFLAGS) 

#housekeeping
all:  $(ALL)
clean:
	rm -rf *.o $(ALL)
