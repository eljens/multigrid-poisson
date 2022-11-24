#!/bin/sh -x

source modules.sh
source paths.sh
source exports.sh

export CUDA_VISIBLE_DEVICES=1

make clean
make

FILENAME=vcycle.txt

MAXITER=10
LEVELS=3

for N in 9 17 33 65 129 257 513 1025
do
    ./bin/driver -x $N -y $N -z $N -maxiter $MAXITER -levels $LEVELS -save 0 -length 5 -stats $FILENAME
    LEVELS=$(($LEVELS + 1))
done


