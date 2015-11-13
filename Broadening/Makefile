SHELL = /bin/bash

.DELETE_ON_ERROR:

.PHONY: all clean


ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS    := $(shell $(ROOTCONFIG) --libs)
ROOTINCDIR  := $(shell $(ROOTCONFIG) --incdir)

CXX       := g++
CXXFLAGS  += -std=c++11 -O2 -Wall -fPIC $(ROOTCFLAGS)
LD        = g++
LDFLAGS   = -O2 $(ROOTLDFLAGS)

INCLUDES  := -I/$(ROOTINCDIR)
LIBS      := $(ROOTLIBS)

##############################################################################
all: AcceptPt2 Pt2Broad

Pt2Broad:Pt2Broad.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@
	
Pt2Broad.o:Pt2Broad.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@	
	
AcceptPt2:AcceptPt2.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

AcceptPt2.o:AcceptPt2.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	@rm -rf AcceptPt2 AcceptPt2.o Pt2Broad Pt2Broad.o
