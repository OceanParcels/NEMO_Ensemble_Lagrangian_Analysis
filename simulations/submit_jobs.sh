#!/bin/bash -l
#
#SBATCH -J NELS           # the name of your job   
#SBATCH -p normal           # request normal partition, job takes > 1 hour (this line can also be left out because 'normal' is the default)  
#SBATCH -t 120:00:00         # time in hh:mm:ss you want to reserve for the job
#SBATCH -n 1               # the number of cores you want to use for the job, SLURM automatically determines how many nodes are needed
#SBATCH -o logs/hatteras.%j.o  # the name of the file where the standard output will be written to. %j will be the jobid determined by SLURM
#SBATCH -e logs/hatteras.%j.e  # the name of the file where potential errors will be written to. %j will be the jobid determined by SLURM
#SBATCH --mail-user=c.m.pierard@uu.nl
#SBATCH --mail-type=ALL

conda activate nemo-ensemble

# for i in {41..48}
# # for i in $(seq 0.01 0.01 0.2)
# do
#    echo "Member: $i"
#    python ensemble_Member_simulation.py -m $i -s 0.1
#    # echo "MEMEMber 50 -- STD: $i"
#    # python ensemble_Member_simulation.py -m 49 -s $i

# done

for j in {43..43}
do
   for i in $(seq 0.01 0.01 0.2)
   do
      # Skip the iterations where i equals 0.05 or 0.15
      if [ $(echo "$i == 0.01" | bc) -eq 1 ] || [ $(echo "$i == 0.1" | bc) -eq 1 ]; then
         continue
      fi
      echo "MEMEMber $j -- STD: $i"
      python ensemble_Member_simulation.py -m $j -s $i
   done
done
echo 'Finished computation.'


# for i in $(seq 0.09 0.01 0.2)
# do
#    # Skip the iterations where i equals 0.05 or 0.15
#    if [ $(echo "$i == 0.01" | bc) -eq 1 ] || [ $(echo "$i == 0.1" | bc) -eq 1 ]; then
#       continue
#    fi
#    echo "MEMEMber 48 -- STD: $i"
#    python ensemble_Member_simulation.py -m 48 -s $i
# done

# echo "MEMEMber 48 -- STD: 0.2"
# python ensemble_Member_simulation.py -m 2 -s 0.2

echo 'Finished computation.'
