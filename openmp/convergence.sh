#!/bin/sh -x

source modules.sh
source paths.sh
source exports.sh

arch=sm_80
MAXITER=30

export CUDA_VISIBLE_DEVICES=1

for PROB in 1 2 3
do
make clean
make GPU_ARCH=$arch PROBLEM=$PROB

FILENAME=vcycle_problem$PROB.txt
rm -rf results/$FILENAME

LEVELS=3
for N in 9 17 33 65 129 257 513
do
    M=$((2*($N-1)+1))
    M=$((2*($M-1)+1))
    M=$((2*($M-1)+1))
    ./bin/multigrid -x $M -y $N -z $N -maxiter $MAXITER -levels $LEVELS -length 0.5 -stats $FILENAME
    LEVELS=$(($LEVELS + 1))
done
done


