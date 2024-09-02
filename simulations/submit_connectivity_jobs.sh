#!/bin/bash -l
#
#SBATCH -J GSnSPG           # the name of your job   
#SBATCH -p normal           # request normal partition, job takes > 1 hour (this line can also be left out because 'normal' is the default)  
#SBATCH -t 120:00:00         # time in hh:mm:ss you want to reserve for the job
#SBATCH -n 1               # the number of cores you want to use for the job, SLURM automatically determines how many nodes are needed
#SBATCH -o logs/connectivity.%j.o  # the name of the file where the standard output will be written to. %j will be the jobid determined by SLURM
#SBATCH -e logs/connectivity.%j.e  # the name of the file where potential errors will be written to. %j will be the jobid determined by SLURM
#SBATCH --mail-user=c.m.pierard@uu.nl
#SBATCH --mail-type=ALL

conda activate nemo-ensemble

for j in {41..50}
do 
   echo "Member $j -- Week span: 12, depth: 1, delta R: 1.75"
   python ensemble_Member_connectivity.py -m $j -d 1 -dr 1.75 -w 12
done

echo "Finished computation."
