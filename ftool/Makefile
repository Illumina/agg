CC=gcc
CXX=g++

# ifndef RDIR
#     $(error RDIR is undefined)
# endif

all: ftool

CFLAGS = -O2  $(ALLFLAGS) 

debug: CFLAGS =  -pg -g -Wall $(ALLFLAGS) 
debug: all 

#linker stuff
HTSDIR = htslib-1.3
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR) -I./ 
LFLAGS = -lz -lpthread
ALLFLAGS =  $(IFLAGS) 
CXX_FLAGS = -std=c++0x -Wno-write-strings
RLIB=$(RDIR)/src/nmath/standalone/libRmath.a

GIT_HASH := $(shell git describe --abbrev=4 --always )

VERSION = 0.x.x
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
version.h:
	echo '#define VERSION "$(VERSION)"' > $@
##bcftools module
filter.o: filter.c 
	$(CC) $(CFLAGS) -c $<  
##agg source code
power.o: power.cpp 
	$(CXX) $(CFLAGS) -c power.cpp $(IFLAGS)
annotate.o: annotate.cpp 
	$(CXX) $(CFLAGS) -c $<   $(IFLAGS)
evaluate.o: evaluate.cpp 
	$(CXX) $(CFLAGS)  -c $<  $(IFLAGS) 
utils.o: utils.cpp utils.h
	g++ -c utils.cpp $(IFLAGS)
ftool: ftool.cpp  annotate.o evaluate.o power.o utils.o filter.o $(HTSLIB)
	$(CXX) $(CXX_FLAGS) -o $@ ftool.cpp *.o $(HTSLIB) $(LFLAGS) $(IFLAGS) 
#housekeeping
all:  $(ALL)
clean:
	rm -rf *.o $(ALL) version.h ftool 

