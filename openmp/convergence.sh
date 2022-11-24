#!/bin/sh -x

source modules.sh
source paths.sh
source exports.sh

make clean
make

FILENAME=vcycle.txt

MAXITER=30
LEVELS=3

for N in 33 65 129 257 513
do
    ./bin/driver -x $N -y $N -z $N -maxiter $MAXITER -levels $LEVELS -save 0 -length 5 -stats $FILENAME
    LEVELS=$(($LEVELS + 1))
done


