#for pitt smp cluster
export GCC_ROOT := /ihome/crc/install/gcc/5.4.0
export PATH := $(GCC_ROOT)/bin:$(PATH)
export CPATH := $(GCC_ROOT)/include:$(CPATH)
export LIBRARY_PATH := $(GCC_ROOT)/lib64:$(LIBRARY_PATH)
export LD_LIBRARY_PATH := $(GCC_ROOT)/lib64:$(LD_LIBRARY_PATH)
export LD_RUN_PATH := $(GCC_ROOT)/lib64:$(LD_RUN_PATH)
export LIBPATH := $(GCC_ROOT)/lib/gcc/x86_64-unknown-linux-gnu/5.4.0:$(LIBPATH)
CXXFLAGS = -I/ihome/crc/install/intel-2017.1.132/eigen/3.3.1/include/eigen3 -std=c++11
#
#targets=p1
#all: $(targets)
#
#clean:
#	rm -vf *~ $(targets) *.o

SOURCES=p1.cpp molecule.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=p1

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $< -o $@

clean:
	rm -vf $(OBJECTS) $(EXECUTABLE)
