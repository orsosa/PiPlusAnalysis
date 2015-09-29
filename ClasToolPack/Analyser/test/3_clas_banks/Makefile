.DELETE_ON_ERROR:

ifndef CLASTOOL
    $(error "Please set the variable CLASTOOL")
endif

ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS   := $(shell $(ROOTCONFIG) --glibs)

CXX         := g++
CXXFLAGS    := -O2 -Wall -fPIC $(ROOTCFLAGS)
LD          := g++
LDFLAGS     := -O2 $(ROOTLDFLAGS)

INCLUDES    := -I../../analysis_lib/include               \
               -I$(CLASTOOL)/include
LIBS        := $(ROOTGLIBS)                               \
               -L$(CLASTOOL)/slib/Linux -lClasTool        \
               -L../../analysis_lib/slib -lTIdentificator

SRC_FILES   := $(wildcard *.cxx)
EXE_FILES   := $(SRC_FILES:%.cxx=%) 


.PHONY: all clean distclean

all: $(EXE_FILES)


%: %.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


clean:
	rm -f $(EXE_FILES:%=%.o) $(EXE_FILES)

distclean: clean
	rm -f particle*.root
