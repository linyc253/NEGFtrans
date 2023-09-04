#!/bin/sh
#SBATCH -A MST111107
#SBATCH -J NEGFtrans
#SBATCH -p ct560            ## Name of queue
#SBATCH -n 500              ## Total number of cores
#SBATCH -c 1       
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH -t 01:00:00         ## Wall time limit (days-hrs:min:sec)

module load compiler/gcc/10.2.0 OpenMPI/4.0.5
export UCX_TLS="ud,dc,rc,sm,self"
export OMPI_MCA_btl="^vader,tcp,openib,uct"
export OMPI_MCA_pml="ucx"
export KMP_AFFINITY="compact,noverbose"
export UCX_NET_DEVICES=mlx5_0:1

mpirun -n $SLURM_NTASKS /home/linyc253/public/NEGFtrans/NEGFtrans_1.0.0
