git clone https://github.com/llvm/llvm-project.git
cd llvm-project
git checkout 710a834c4c822c5c444fc9715785d23959f5c645
cd ..

module load cmake/3.23.2
module load gcc/11.3.0-binutils-2.38
export CC=`which gcc`
export CXX=`which g++`
cmake llvm-project/llvm -DLLVM_ENABLE_PROJECTS='clang;lld' -DLLVM_ENABLE_RUNTIMES='openmp' -DCMAKE_BUILD_TYPE=Release -DLLVM_TARGETS_TO_BUILD='host;NVPTX;AMDGPU' -DCLANG_OPENMP_NVPTX_DEFAULT_ARCH=sm_80 -DLIBOMPTARGET_NVPTX_COMPUTE_CAPABILITIES=80,70 -DGCC_INSTALL_PREFIX=/appl/gcc/11.3.0-binutils-2.38/
make -j 24
