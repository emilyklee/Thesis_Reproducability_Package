#!/bin/bash
#SBATCH -J pythonCode                        # name of job
#SBATCH -A niemeyek                          # name of my sponsored account, e.g. class or research group
#SBATCH -p mime4                             # name of partition or queue
#SBATCH -F ./nrg-nodes
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH -o bfmTimerRuns.out                  # name of output file for this submission script
#SBATCH -e bfmTimerRuns.err                  # name of error file for this submission script
#SBATCH --mail-type=BEGIN,END,FAIL           # send email when job begins, ends or aborts
#SBATCH --mail-user=kleeem@oregonstate.edu   # send email to this address
echo $SLURM_JOB_ID
python timer_simulations.py species_removed.txt
