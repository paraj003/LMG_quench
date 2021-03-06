#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -n 4
#SBATCH -N 1
#SBATCH -A physics-hi
#SBATCH --mail-user=paraj@umd.edu
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=6999
##SBATCH --exclusive
#SBATCH --array=1-14
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

L=600
LAarr=np.concatenate([np.linspace(10,100,10),np.linspace(150,300,4)],axis=0)#np.concatenate([np.linspace(10,100,10),np.linspace(200,1000,9)],axis=0)
La=LAarr[int(os.environ[\"SLURM_ARRAY_TASK_ID\"])-1]  #Set system size.
S=L/2
SA=La/2
SB=S-SA
directory='data/CGmats'
LMG.CGmatrix_to_file(SA,SB,S,directory)
 ">>CGmatscript-$SLURM_JOB_ID.py


python3 CGmatscript-$SLURM_JOB_ID.py

rm *$SLURM_JOB_ID.py

hostname
date
