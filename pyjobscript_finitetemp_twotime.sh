#!/bin/bash
##SBATCH -p debug
#SBATCH -t 48:00:00
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -A physics-hi
#SBATCH --mail-user=paraj@umd.edu
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=6999
##SBATCH --exclusive
#SBATCH --array=1-12
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

Larr=np.concatenate([np.linspace(100,1000,10),np.linspace(2000,9000,8)],axis=0)
L=Larr[int(os.environ[\"SLURM_ARRAY_TASK_ID\"])-1]  #Set system size.
paramvals=LMG.Ham_params(N=L,S=L/2,J=1.,γz=1.,γy=0.,Γ=0.7)
β=LMG.finitetemp_criticalβ(paramvals)
dt=1. #time step
t1arr=np.array([20.])
t2arr=np.linspace(10,50,int((50-10)/dt)+1)
Az=1.
Ay=0.
">>LMG-params-$SLURM_JOB_ID.py


cat LMG-params-$SLURM_JOB_ID.py cluster_finitetemp_twotime_LMGsim.py > temp-cluster_finitetemp_twotime_LMGsim-$SLURM_JOB_ID.py


python3 temp-cluster_finitetemp_twotime_LMGsim-$SLURM_JOB_ID.py

rm *$SLURM_JOB_ID.py

hostname
date
