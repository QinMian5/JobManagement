#!/bin/bash

#SBATCH --partition=p_pamish
#SBATCH --job-name="op_${QBAR.CENTER}"
#SBATCH --output=job.log
#SBATCH --error=job.error
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

if [[ -d "/scratch/pamish1/mian" ]]; then
    module purge
    module load gcc-9.2.0/openmpi/4.1.0
fi

gmx grompp -p ../../../../topol.top -f grompp.mdp -c ../../conf.gro -r ../../conf.gro -maxwarn 10

mpirun -np 8 gmx_mpi mdrun -s topol.tpr -ntomp 1 -cpt 10 -dd 2 2 2 -op
echo -e "0\n0" | gmx trjconv -f traj_comp.xtc -s ../../topol.tpr -pbc mol -center

(
    cd post_processing || exit
    OrderParameters post_processing.dat
)

(
    cd post_processing_with_PI || exit
    OrderParameters post_processing_with_PI.dat
)

(
    cd post_processing_chillplus || exit
    OrderParameters post_processing_chillplus.dat
)
