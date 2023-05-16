#!/bin/bash
#SBATCH -p small-g
#SBATCH -A project_465000434
#SBATCH --time=0:44:59
#SBATCH --nodes=1
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=460GB
#SBATCH --job-name=transfer_clang
#SBATCH --output=%x.%j.out

source paths/lumi/clang.sh
export COMPILER=clang
export GPU=MI250X
export OMP_NUM_THREADS=63
export OMP_WAIT_POLICY=ACTIVE

export ROCR_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

for size in 9 17 33 65 129 257 513 1025;
do
	filename=results/$GPU/clang/transfer_$size.txt	
	rm -rf $filename
	make clean; make SAVEP2P='"'${filename}'"' APP=multi_gpu_weak	
	./bin/multi_gpu_weak -x $size -y $size -z $size -levels 1 -maxiter 100 -tol 1e-10

	filename1=results/$GPU/clang/transfer_p2p_$size.txt
	rm -rf $filename1
	make clean; make PEER2PEER=1 SAVEP2P='"'${filename1}'"' APP=multi_gpu_weak	
	./bin/multi_gpu_weak -x $size -y $size -z $size -levels 1 -maxiter 100 -tol 1e-10
	
done
