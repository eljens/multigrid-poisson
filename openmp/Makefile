APP?=multigrid

ifndef GPU_ARCH
GPU_ARCH = sm_80
$(info GPU_ARCH is not set and defaults to $(GPU_ARCH))
endif

ARCH_PATH=/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/libomptarget/libomptarget-nvptx-$(GPU_ARCH).bc

TARGETFLAGS = -fopenmp --offload-arch=$(GPU_ARCH) --libomptarget-nvptx-bc-path=$(ARCH_PATH)

CXXFLAGS = -Wall -g -O3 -std=c++17 $(TARGETFLAGS)

CXX = /work3/s174515/LLVM/bin/clang++

#### Standard settings which are identical for all compilers
# Add flags start start with -I here (include path)
CPPFLAGS = -Iinclude -I/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/runtime/src/
# Add flags that start with -L here (link path)
LDFLAGS = -Llib -L/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/runtime/src/ -lomp -L/appl/gcc/11.3.0-binutils-2.38/lib64/ -Wl,-rpath,/appl/gcc/11.3.0-binutils-2.38/lib64/
# Add flags stat start with -l here (libbrary)
LDLIBS = -lm -lstdc++ -lpoisson

# Source files
SRC = $(wildcard src/*.cpp)

# Object files
OBJ = $(patsubst src/%.cpp,lib/%.o,$(SRC))

TARGET?=bin/$(APP)

bin/$(APP): lib/$(APP).o lib/libpoisson.a
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

lib/$(APP).o: drivers/$(APP).cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

lib/libpoisson.a: $(OBJ)
	ar rcs lib/libpoisson.a $(OBJ)

lib/%.o: src/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ -c $<

archive:
	make TARGET=lib/libpoisson.a

clean:
	-rm -f top3d *.o lib/*.o lib/*.a bin/*

info:
	$(info "COMPILER=$(COMPILER)")
	$(info "CXX=$(CXX)")
	$(info "CXXFLAGS=$(CXXFLAGS)")
	$(info "LDFLAGS=$(LDFLAGS)")
	$(info "LDLIBS=$(LDLIBS)")
	$(info "DEFS=$(DEFS)")

