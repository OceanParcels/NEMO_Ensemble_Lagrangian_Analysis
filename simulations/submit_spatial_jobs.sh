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
   echo "MEMEMber $j -- Delta r: 2.0"
   python ensemble_Member_spatial.py -m $j -dr 2.0
   
done

# for j in {26..50}
# do
#    for i in 0.1 1.0 2.0
#    do
#       echo "MEMEMber $j -- Delta r: $i"
#       python ensemble_Member_spatial.py -m $j -dr $i
#    done
# done

# for j in {43..43}
# do
#    for i in $(seq 0.01 1.0 0.1)
#    do
#       # # Skip the iterations where i equals 0.05 or 0.15
#       # if [ $(echo "$i == 0.01" | bc) -eq 1 ] || [ $(echo "$i == 0.1" | bc) -eq 1 ]; then
#       #    continue
#       # fi
#       # echo "MEMEMber $j -- STD: $i"
#       # python ensemble_Member_simulation.py -m $j -s $i
#    done
# done
echo "Finished computation."
