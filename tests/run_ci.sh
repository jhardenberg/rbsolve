#!/bin/bash
set -e

echo "Starting CI/CD Integrated Test for rbsolve"

cd src/

# 1. Modify param.h for low resolution: Nx=32, Ny=16 (vertical), Nz=1
# Usually kx is about 2/3 of Nx to dealias. Let's set kx=21, kz=0
echo "Configuring low-resolution grid in param.h"
sed -i.bak 's/parameter *(Nx=.*)/parameter (Nx=32,Ny=16,Nz=1)/g' param.h
sed -i.bak 's/parameter *(kx=.*)/parameter (kx=21,kz=0)/g' param.h
sed -i.bak 's/parameter  *(NPROC=4)/parameter  (NPROC=2)/g' param.h


# 2. Modify param_namelist for Ra=10.4 and short runtime
echo "Configuring param_namelist for Ra=10.4"
sed -i.bak 's/Ra=.*!/Ra=10.4d0                  !/g' param_namelist
sed -i.bak 's/ttot=.*/ttot=100/g' param_namelist
sed -i.bak 's/nsave=.*/nsave=50/g' param_namelist

# 3. Clean and Compile
echo "Enabling MPI in Makefile"
sed -i.bak 's/NOMPI = 1/#NOMPI = 1/g' Makefile

echo "Compiling inicond and rb"
make clean
make inicond
make

# 4. Cleanup old files and Run
echo "Executing inicond"
rm -f rb_*.nc nrec.d coord
mpirun -np 2 ./inicond

echo "Executing main rb simulation"
mpirun -np 2 ./rb

echo "CI/CD Test Completed Successfully!"
