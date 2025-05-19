# Introduction
The program solve the electron density, local density of state, and transmission coefficient of an open system using Non-Equilibrium Green's Function (NEGF) method. In particular, the program computes the quantities based on plane wave basis in transverse direction (xy), and real space grid in transport direction (z).
# Installation on Taiwania 3
Get the source code by
```
git clone https://github.com/linyc253/NEGFtrans.git
```
Then, go to `build` directory, and compile with the following command
```
module load gcc/10.5.0 intelmpi/2021.11
make
```
# Execution on Taiwania 3
Use the following job script
```
#!/bin/sh
#SBATCH -A XXXXXXXX
#SBATCH -J NEGFtrans
#SBATCH -p ctest
#SBATCH -n 300          ## Total cores
#SBATCH -c 1            ## Without this option, the controller will just try to allocate one processor per task.
##SBATCH -N 1            ## Number of nodes
#SBATCH -o %j.out
#SBATCH -e %j.err
##SBATCH -t 05:00:00     ## Wall time limit (days-hrs:min:sec)

module load gcc/10.5.0 intelmpi/2021.11

mpirun /[your_path]/NEGFtrans/build/NEGFtrans.x
```
For detailed input parameter setting, please refer to `doc/Manual.pdf`

# Test Case
Go to `testcase` directory, and prepare the input potential by
```
gfortran potential_QPC.f90 -o potential_QPC.x
./potential_QPC.x
```
Then, you should find a file named `POTENTIAL` in current directory, and you can run the code by submitting the job script. After the calculation done, there should be a file named `TRANSMISSION`. You can then visualize it by either
```
xmgrace TRANSMISSION
```
or the python script (you may need to do some environment setting before running python on Taiwania 3)
```
python plot_trans.py
```
