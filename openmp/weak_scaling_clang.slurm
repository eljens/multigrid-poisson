#!/bin/bash
#SBATCH -p small-g
#SBATCH -A project_465000434
#SBATCH --time=0:44:59
#SBATCH --nodes=1
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=400GB
#SBATCH --job-name=weak_scaling_clang
#SBATCH --output=%x.%j.out

source paths/lumi/clang.sh
export COMPILER=clang
export GPU=MI250X
export OMP_NUM_THREADS=63
export OMP_WAIT_POLICY=ACTIVE

make realclean; make APP=multi_gpu_weak PEER2PEER=1
mv ./bin/multi_gpu_weak ./bin/weak_scaling_p2p
make clean;make APP=multi_gpu_weak
mv ./bin/multi_gpu_weak ./bin/weak_scaling
for level in 1 2 3 4 5 6 7;
do
	filename=results/$GPU/clang/weak_scaling_p2p_$level.txt
	filename1=results/$GPU/clang/weak_scaling_$level.txt
	rm -rf $filename $filename1
	command="./bin/weak_scaling_p2p -x 641 -y 1025 -z 641 -levels $level -length 4 -maxiter 30 -tol 1e-6 -stats $filename"
	command1="./bin/weak_scaling -x 641 -y 1025 -z 641 -levels $level -length 4 -maxiter 30 -tol 1e-6 -stats $filename1"
	export ROCR_VISIBLE_DEVICES=0
	$command
	$command1
	export ROCR_VISIBLE_DEVICES=0,1
	$command
	$command1
	export ROCR_VISIBLE_DEVICES=0,1,2,3
	$command
	$command1
        export ROCR_VISIBLE_DEVICES=0,1,2,3,4,5,6,7
        $command
        $command1
done
