APP?=multi_gpu

#### User must provide a desired compiler
ifndef COMPILER
CONFIGFILES = $(wildcard *config/*.mkf)
COMPILERS = $(patsubst config/%.mkf,%,$(CONFIGFILES))
$(info The following compilers are supported: $(COMPILERS))
$(error Please set COMPILER before running make. You can use export COMPILER=<compiler>)
else
ifeq ("$(wildcard config/$(COMPILER).mkf)","")
$(error "config/$(COMPILER).mkf does not exist")
endif
include config/$(COMPILER).mkf
endif

ifndef GPU
GPUFILES = $(wildcard *config/gpu/*.mkf)
GPUS = $(patsubst config/gpu/%.mkf,%,$(GPUFILES))
$(info The following GPUs are supported: $(GPUS))
$(error Please set GPU before running make. You can use export GPU=<gpu_name>)
else
ifeq ("$(wildcard config/gpu/$(GPU).mkf)","")
$(error "config/gpu/$(GPU).mkf does not exist")
endif
include config/gpu/$(GPU).mkf
endif

ifndef PROBLEM
PROBLEM=1
endif

CXXFLAGS+= -DPROBLEM$(PROBLEM)

ifdef SCHED
CXXFLAGS+= -DSCHED=$(SCHED)
endif

ifdef BLOCK_SIZE
ifeq ($(shell test $(BLOCK_SIZE) -gt 0; echo $$?),0)
CXXFLAGS+= -DBLOCK_SIZE=$(BLOCK_SIZE)
endif
endif

ifdef DIST
CXXFLAGS+= -DDIST=$(DIST)
endif

#### Standard settings which are identical for all compilers
# Add flags start start with -I here (include path)
CPPFLAGS += -Iinclude
# Add flags that start with -L here (link path)
LDFLAGS += -Llib 
# Add flags stat start with -l here (libbrary)
LDLIBS += -lpoisson -lm -lstdc++

# Source files
SRC = $(wildcard src/*.cpp)

# Object files
OBJ = $(patsubst src/%.cpp,lib/%.o,$(SRC))

# Drivers
APPS = $(wildcard drivers/*.cpp)
APPSCLEAN = $(patsubst drivers/%.cpp,lib/%.o,$(APPS))

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
	-rm -f bin/* $(APPSCLEAN)

realclean:
	-rm -f lib/*.o lib/*.a bin/*

info:
	$(info "COMPILER=$(COMPILER)")
	$(info "CXX=$(CXX)")
	$(info "CXXFLAGS=$(CXXFLAGS)")
	$(info "LDFLAGS=$(LDFLAGS)")
	$(info "LDLIBS=$(LDLIBS)")
	$(info "DEFS=$(DEFS)")

