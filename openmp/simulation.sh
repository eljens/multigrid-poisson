#!/bin/sh -x

source modules.sh
source paths.sh
source exports.sh

arch=sm_80
MAXITER=30

export CUDA_VISIBLE_DEVICES=1

if 0
then
PROB=1
make clean
make GPU_ARCH=$arch PROBLEM=$PROB
./bin/multigrid -x 513 -y 513 -z 257 -maxiter $MAXITER -levels 9 -length 2 -save
fi

PROB=2
make clean
make GPU_ARCH=$arch PROBLEM=$PROB
./bin/multigrid -x 513 -y 513 -z 257 -maxiter $MAXITER -levels 9 -length 5 -save

if 0
then
PROB=3
make clean
make GPU_ARCH=$arch PROBLEM=$PROB
./bin/multigrid -x 513 -y 513 -z 513 -maxiter $MAXITER -levels 9 -length 2.5 -save
fi