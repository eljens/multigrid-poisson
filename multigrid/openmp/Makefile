APP?=driver

ifndef GPU_ARCH
GPU_ARCH = sm_80
$(info GPU_ARCH is not set and defaults to $(GPU_ARCH))
endif

ARCH_PATH=/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/libomptarget/libomptarget-nvptx-$(GPU_ARCH).bc

TARGETFLAGS = -fopenmp --offload-arch=$(GPU_ARCH) --libomptarget-nvptx-bc-path=$(ARCH_PATH)

CFLAGS += -Wall -g -std=c++17 -O3 $(TARGETFLAGS)
CFLAGS += -I/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/runtime/src/

LIBS = -lm -L/appl/gcc/11.3.0-binutils-2.38/lib64/ -Wl,-rpath,/appl/gcc/11.3.0-binutils-2.38/lib64/
LIBS += -L/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/runtime/src/ -lomp


$(info Compiler flags are: $(CFLAGS))
$(info Linker flags are: $(LIBS))

CC = /work3/s174515/LLVM/bin/clang
CXX = /work3/s174515/LLVM/bin/clang++

OBJ = parser.o

all: $(APP)

$(APP): $(APP).cpp $(OBJ)
	$(CXX) $(CFLAGS) -o $@ $^ $(LIBS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

clean:
	-rm -f $(APP) *.core *.o 
