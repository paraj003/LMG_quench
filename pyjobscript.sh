#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -n 5
#SBATCH -N 1
#SBATCH -A physics-hi
#SBATCH --mail-user=paraj@umd.edu
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=6999
##SBATCH --exclusive
#SBATCH --array=1-2
. ~/.profile
export MKL_NUM_THREADS=$SLURM_NTASKS
export FFTW_NUM_THREADS=$SLURM_NTASKS

module unload python
module load python/3.5.1

echo "

import numpy as np
from scipy import linalg as LA
import scipy.sparse.linalg as spla
import os
import mod_LMG as LMG
import h5py

Larr=np.linspace(100,1000,10)
L=Larr[int(os.environ[\"SLURM_ARRAY_TASK_ID\"])]  #Set system size.
paramvals0=LMG.Ham_params(N=L,S=L/2,J=1,γ=0.1,Γ=4)
paramvalsf=LMG.Ham_params(N=L,S=L/2,J=1,γ=0.1,Γ=1.02)
dt=0.2 #time step
Tf=20 # final time step
Nsteps=int(Tf/dt) 
tarr=np.arange(dt,Tf+dt,dt)
">>LMG-params-$SLURM_JOB_ID.py


cat LMG-params-$SLURM_JOB_ID.py cluster_LMGsim.py > temp-cluster_LMGsim-$SLURM_JOB_ID.py


python3 temp-cluster_LMGsim-$SLURM_JOB_ID.py

rm *$SLURM_JOB_ID.py

hostname
date
