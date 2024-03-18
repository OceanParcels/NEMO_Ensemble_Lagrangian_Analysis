#!/bin/bash -l
#
#SBATCH -J P-Ensemble           # the name of your job   
#SBATCH -p normal           # request normal partition, job takes > 1 hour (this line can also be left out because 'normal' is the default)  
#SBATCH -t 120:00:00         # time in hh:mm:ss you want to reserve for the job
#SBATCH -n 1               # the number of cores you want to use for the job, SLURM automatically determines how many nodes are needed
#SBATCH -o P-Ensemble.%j.o  # the name of the file where the standard output will be written to. %j will be the jobid determined by SLURM
#SBATCH -e P-Ensemble.%j.e  # the name of the file where potential errors will be written to. %j will be the jobid determined by SLURM
#SBATCH --mail-user=c.m.pierard@uu.nl
#SBATCH --mail-type=ALL

conda activate nemo-ensemble

for i in {1..50}
do
   echo "Member: $i"
   python ensemble_Member_simulation.py -m $i

done

echo 'Finished computation.'