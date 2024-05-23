#!/bin/bash
#SBATCH -J pyED
#SBATCH -o output.%a
#SBATCH -e error.%a
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --array=1-48

echo $SLURM_ARRAY_TASK_ID
subfolder=$(awk "NR==$SLURM_ARRAY_TASK_ID" job.list)
echo $subfolder
cd ../calc/$subfolder

pyED -s result.p
