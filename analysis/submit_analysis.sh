#!/bin/bash -l
#
#SBATCH -J Analizando           # the name of your job   
#SBATCH -p normal           # request normal partition, job takes > 1 hour (this line can also be left out because 'normal' is the default)  
#SBATCH -t 20:00:00         # time in hh:mm:ss you want to reserve for the job
#SBATCH -n 1               # the number of cores you want to use for the job, SLURM automatically determines how many nodes are needed
#SBATCH -o logs/analysis.%j.o  # the name of the file where the standard output will be written to. %j will be the jobid determined by SLURM
#SBATCH -e logs/analysis.%j.e  # the name of the file where potential errors will be written to. %j will be the jobid determined by SLURM
#SBATCH --mail-user=c.m.pierard@uu.nl
#SBATCH --mail-type=ALL

conda activate nemo-ensemble

python connectivity_mixture_stats_spatial.py
# python connectivity_single_member_stats.py
# python connectivity_mixture_stats_temporal.py
# python mixture_distributions_temp.py
# python mixture_distributions_spatial.py

echo 'Finished computation.'
