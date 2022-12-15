CXX=clang++

ROOTCONFIG := root-config

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += -I${FAIRROOTPATH}/include
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(VMCWORKDIR)/r3bdata
CFLAGS += -I$(VMCWORKDIR)/tracking
CFLAGS += -I$(VMCWORKDIR)/r3bdata/fibData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/califaData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/trackerData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/neulandData
CFLAGS += -I$(VMCWORKDIR)/r3bdata/landData
CFLAGS += -I$(SIMPATH)/include
CFLAGS += -I$(ROOT_INCLUDE_PATH)
CFLAGS += -I$(ROOT_INCLUDE_DIR)
CFLAGS += --std=c++11 -g -O0 -fexceptions
#CFLAGS += -g -Wall -W -Wconversion -Wshadow -Wcast-qual -Wwrite-strings 

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += -lEG $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(ROOT_LIBRARY_DIR) -L$(FAIRROOTPATH)/lib




invariant_mass : invariant_mass.cpp 
	$(CXX) $(CFLAGS) $(LDFLAGS) -v -O0 $^ -o $@ 


clean:
	rm -f *.o
	rm -f invariant_mass
	rm -rf *.dSYM
