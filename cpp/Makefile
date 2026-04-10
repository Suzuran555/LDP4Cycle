CXXFLAGS    = -march=native -O3
CXXFLAGS    += -std=c++11
CXXFLAGS    += -isystem ./include
CXXFLAGS    += -MMD
CXXFLAGS    += -fopenmp    # Added for OpenMP compilation

LDFLAGS     += -fopenmp    # Added for OpenMP linking

# Note: Make uses $(CC) for C and $(CXX) for C++. 
# Since you are compiling C++, it's technically more standard to use CXX, 
# but keeping CC=g++ works perfectly fine with your implicit rule below.
CC  = g++
CXX = g++

TARGET  = SubgraphShuffle

all: $(TARGET)

SubgraphShuffle: SubgraphShuffle.o

clean:
	rm -f *~ $(TARGET) *.d *.o

-include *.d