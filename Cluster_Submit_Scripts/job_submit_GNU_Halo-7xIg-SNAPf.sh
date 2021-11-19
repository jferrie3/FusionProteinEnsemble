#!/bin/bash
# Job name:
#SBATCH --job-name=Halo7xIgSnap
#
# Account:
#SBATCH --account=fc_3dgenome
#
# Partition:
#SBATCH --partition=savio
#
# Nodes
#SBATCH --nodes=1
#
# Tasks per node
#SBATCH --cpus-per-task=1
#
# Wall clock limit:
#SBATCH --time=05:00:00
#
## Command(s) to run:
module load gnu-parallel/2019.03.22
module load python/2.7
## Default is to run in HOME directory
export WDIR=/global/home/users/jferrie/TG_HaloSnap/Halo-7xIg-SNAPf/
cd $WDIR
# set the number of jobs
export JOBS_PER_NODE=$(( $SLURM_CPUS_ON_NODE / $SLURM_CPUS_PER_TASK ))
## For a single node
parallel --jobs $SLURM_CPUS_ON_NODE --joblog task.log --resume --progress sh run_cmd_7.sh {} ::: `seq 100`
## For parallelizing across multiple nodes
## echo $SLURM_JOB_NODELIST |sed s/\,/\\n/g > hostfile
## parallel --jobs $JOBS_PER_NODE --slf hostfile --wd $WDIR --joblog task.log --resume --progress sh run_cmd.sh {} ::: `seq 100`