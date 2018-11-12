#!/bin/bash
#SBATCH -t 7:00:00
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

Larr=[100]#np.concatenate([np.linspace(100,1000,10),np.linspace(2000,5000,4)],axis=0)
L=Larr[int(os.environ[\"SLURM_ARRAY_TASK_ID\"])-1]  #Set system size.
paramvals0=LMG.Ham_params(N=L,S=L/2,J=1.,γz=1.,γy=0.,Γ=4)
paramvalsf=LMG.Ham_params(N=L,S=L/2,J=1.,γz=1.0,γy=0.0,Γ=1)
dt=1. #time step
t1arr=np.array([15.])
t2arr=np.arange(10,50+dt,dt)
">>LMG-params-$SLURM_JOB_ID.py


cat LMG-params-$SLURM_JOB_ID.py cluster_twotime_LMGsim.py > temp-cluster_twotime_LMGsim-$SLURM_JOB_ID.py


python3 temp-cluster_twotime_LMGsim-$SLURM_JOB_ID.py

rm *$SLURM_JOB_ID.py

hostname
date
