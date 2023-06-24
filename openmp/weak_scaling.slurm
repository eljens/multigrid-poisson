#!/bin/bash
#SBATCH -p small-g
#SBATCH -A project_465000434 
#SBATCH --time=0:44:59 
#SBATCH --nodes=1  
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=400GB
#SBATCH --job-name=weak_scaling
#SBATCH --output=%x.%j.out

source paths/lumi/clang.sh
export COMPILER=clang
export GPU=MI250X
export OMP_NUM_THREADS=63 #$SLURM_JOB_CPUS_PER_NODE
export OMP_WAIT_POLICY=ACTIVE
echo $SLURM_JOB_CPUS_PER_NODE

make clean; make APP=multi_gpu_weak
cp ./bin/multi_gpu_weak ./bin/weak_scaling
for level in 1 2 3 4 5 6 7;
do
	filename=results/$GPU/weak_scaling_$level.txt
	rm -rf $filename
	command="./bin/weak_scaling -x 1025 -y 1281 -z 1025 -levels $level -length 4 -maxiter 30 -tol 1e-6 -stats $filename"
	export ROCR_VISIBLE_DEVICES=0
	$command
	export ROCR_VISIBLE_DEVICES=0,1
	$command
	export ROCR_VISIBLE_DEVICES=0,1,2,3
	$command
	export ROCR_VISIBLE_DEVICES=0,1,2,3,4,5,6,7
	$command
done
