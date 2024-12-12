#!/bin/bash -l
#
#SBATCH -J espacio           # the name of your job   
#SBATCH -p normal           # request normal partition, job takes > 1 hour (this line can also be left out because 'normal' is the default)  
#SBATCH -t 120:00:00         # time in hh:mm:ss you want to reserve for the job
#SBATCH -n 1               # the number of cores you want to use for the job, SLURM automatically determines how many nodes are needed
#SBATCH -o logs/hatteras.%j.o  # the name of the file where the standard output will be written to. %j will be the jobid determined by SLURM
#SBATCH -e logs/hatteras.%j.e  # the name of the file where potential errors will be written to. %j will be the jobid determined by SLURM
#SBATCH --mail-user=c.m.pierard@uu.nl
#SBATCH --mail-type=ALL

conda activate nemo-ensemble

for j in {1..50}
do
   for i in 0.01 1.0 2.0
   do
      echo "MEMEMber $j -- Delta r: $i"
      python ensemble_Member_spatial.py -m $j -dr $i
   done
done
echo "Finished computation."

# echo "MEMEMber 39 -- Delta r: 2.0"
# python ensemble_Member_spatial.py -m 39 -dr 2.0
   

# for j in {1..50}
# do
#    echo "MEMEMber $j -- Delta r: 2.0"
#    python ensemble_Member_spatial.py -m $j -dr 2.0
   
# done

# for j in {26..50}
# do
#    for i in 0.1 1.0 2.0
#    do
#       echo "MEMEMber $j -- Delta r: $i"
#       python ensemble_Member_spatial.py -m $j -dr $i
#    done
# done