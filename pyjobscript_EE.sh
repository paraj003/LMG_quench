#!/bin/bash
#SBATCH -p debug
#SBATCH -t 15:00
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -A physics-hi
#SBATCH --mail-user=paraj@umd.edu
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=6999
##SBATCH --exclusive
#SBATCH --array=1-1
. ~/.profile
export MKL_NUM_THREADS=$SLURM_NTASKS
export FFTW_NUM_THREADS=$SLURM_NTASKS

module unload python
module unload java
module load python/3.5.1


echo "

import numpy as np
from scipy import linalg as LA
import scipy.sparse.linalg as spla
import os
import mod_LMG_v1 as LMG
import h5py

ΔL=10
L=2000#np.concatenate([np.arange(100,1000,100),np.arange(1000,2000,1000)])#choose even
part_szarr=np.concatenate([np.linspace(10,100,10)]),np.linspace(200,1000,9)],axis=0)#np.arange(ΔL,int(L/2)+ΔL,ΔL)
dt=1.
tarr=np.linspace(1,50,int((50-1)/dt)+1)#np.arange(1,50+dt,dt)
Stot=L/2
paramvals0=LMG.Ham_params(N=L,S=Stot,J=1.,γz=1.,γy=0.,Γ=4.)
paramvalsf=LMG.Ham_params(N=L,S=Stot,J=1.,γz=1.0,γy=0.0,Γ=1.)
">>LMG-params-EE-$SLURM_JOB_ID.py


cat LMG-params-$SLURM_JOB_ID.py cluster_EE_LMGsim.py > temp-cluster_EE_LMGsim-$SLURM_JOB_ID.py


python3 temp-cluster_EE_LMGsim-$SLURM_JOB_ID.py

rm *$SLURM_JOB_ID.py

hostname
date
