#linker stuff
HTSLIB = htslib-1.2.1/
IFLAGS = -I$(HTSLIB)
LFLAGS = -lz -lpthread -lm -ldl 

ALLFLAGS = -std=c++0x  -Wno-write-strings  $(IFLAGS) 

release: CFLAGS = -O3  $(ALLFLAGS) 
release: all

debug: CFLAGS =  -pg -g -Wall $(ALLFLAGS) 
debug: all 

ALL=agg

##objects libs etc
version.h: 
	git log --pretty=format:'#define VERSION "%h"' -n 1 > version.h
	echo >> version.h
utils.o: utils.cpp utils.h
	$(CXX) $(CFLAGS) -c utils.cpp 
agg_genotyper.o: agg_genotyper.cpp agg_genotyper.h $(HTSLIB)
	$(CXX) $(CFLAGS) -c agg_genotyper.cpp  
agg_utils.o: agg_utils.cpp agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $< 
agg_ingest2.o: agg_ingest2.cpp agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
vcfmerge.o: vcfmerge.c htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
vcfnorm.o: vcfnorm.c htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
vcmp.o: vcmp.c  htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  
agg_ingest1.o: agg_ingest1.cpp agg_ingest1.h agg.h htslib-1.2.1/
	$(CXX) $(CFLAGS) -c $<  

##binaries
agg: agg.cpp vcfnorm.o vcmp.o vcfmerge.o agg_ingest2.o utils.o agg_utils.o agg_genotyper.o  agg_ingest1.o version.h $(HTSLIB)/libhts.a
	$(CXX) $(CFLAGS)  -o agg agg.cpp vcfnorm.o vcmp.o vcfmerge.o agg_ingest2.o agg_genotyper.o utils.o agg_utils.o  agg_ingest1.o $(HTSLIB)/libhts.a  $(LFLAGS) 

#housekeeping
all:  $(ALL)
clean:
	rm -rf *.o $(ALL)
