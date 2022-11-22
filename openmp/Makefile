APP?=driver

ifndef GPU_ARCH
GPU_ARCH = sm_80
$(info GPU_ARCH is not set and defaults to $(GPU_ARCH))
endif

ARCH_PATH=/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/libomptarget/libomptarget-nvptx-$(GPU_ARCH).bc

TARGETFLAGS = -fopenmp --offload-arch=$(GPU_ARCH) --libomptarget-nvptx-bc-path=$(ARCH_PATH)

CFLAGS += -Wall -g -O3 $(TARGETFLAGS)
CFLAGS += -I/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/runtime/src/

LIBS = -lm -L/appl/gcc/11.3.0-binutils-2.38/lib64/ -Wl,-rpath,/appl/gcc/11.3.0-binutils-2.38/lib64/
LIBS += -L/work3/s174515/LLVM/runtimes/runtimes-bins/openmp/runtime/src/ -lomp

#$(info Compiler flags are: $(CFLAGS))
#$(info Linker flags are: $(LIBS))

CC = /work3/s174515/LLVM/bin/clang
CXX = /work3/s174515/LLVM/bin/clang++

SRCDIR   = src
OBJDIR   = lib
BINDIR   = bin
INCLDIR  = include

OBJ = parser.o problem_definition.o settings.o domainsettings.o

all: $(BINDIR)/$(APP)

$(BINDIR)/$(APP): $(APP).cpp $(OBJDIR)/parser.o $(OBJDIR)/problem_definition.o $(OBJDIR)/settings.o $(OBJDIR)/domainsettings.o
	$(CXX) $(CFLAGS) -std=c++17 -o $@ $^ $(LIBS)

$(OBJDIR)/parser.o: $(SRCDIR)/parser.cpp
	$(CXX) $(CFLAGS) -std=c++17 -o $@ -c $<

$(OBJDIR)/settings.o: $(SRCDIR)/settings.cpp
	$(CXX) $(CFLAGS) -std=c++17 -o $@ -c $<

$(OBJDIR)/domainsettings.o: $(SRCDIR)/domainsettings.cpp
	$(CXX) $(CFLAGS) -std=c++17 -o $@ -c $<

$(OBJDIR)/problem_definition.o: $(SRCDIR)/problem_definition.cpp
	$(CXX) $(CFLAGS) -std=c++17 -o $@ -c $<

clean:
	-rm -f $(BINDIR)/$(APP) *.core *.o 
	-rm -f $(OBJDIR)/*.o

realclean:
	make clean
	-rm -rf results/*.vtk

