export LUMI=1
module load rocm
export LD_LIBRARY_PATH=~/GCC/13/gcc-offload/install/lib64:$LD_LIBRARY_PATH
export PATH=/pfs/lustrep2/users/rydahlan/GCC/13/gcc-offload/install/bin:$PATH
export ROCR_VISIBLE_DEVICES=0
export OMP_TARGET_OFFLOAD=MANDATORY
echo $PATH
echo $LD_LIBRARY_PATH
