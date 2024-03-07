#!/bin/bash -l
#
#SBATCH -J NEMO_file_check           # the name of your job   
#SBATCH -p normal           # request normal partition, job takes > 1 hour (this line can also be left out because 'normal' is the default)  
#SBATCH -t 120:00:00         # time in hh:mm:ss you want to reserve for the job
#SBATCH -n 1               # the number of cores you want to use for the job, SLURM automatically determines how many nodes are needed
#SBATCH -o NEMO_file_check.%j.o  # the name of the file where the standard output will be written to. %j will be the jobid determined by SLURM
#SBATCH -e NEMO_file_check.%j.e  # the name of the file where potential errors will be written to. %j will be the jobid determined by SLURM
#SBATCH --mail-user=c.m.pierard@uu.nl
#SBATCH --mail-type=ALL

echo 'Running corrupted_files_check_NEMO_Ensemble.py'
# cd ${HOME}/3DModelling_SouthAtlantic/

python3 corrupted_files_check_NEMO_Ensemble.py

echo 'Finished computation.'