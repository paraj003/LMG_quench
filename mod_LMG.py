# Module containing all definitions necessary to run a quench for the LMG model. Use python 3.
import numpy as np
import os
import h5py



# Class definition to define Hamiltonian

#define Hamiltonian parameters
class Ham_params:
    def __init__(self, N:int,S:float,J:float,γ:float,Γ:float):
        self.N=N #number of spins, keep it even
        self.S=S #spin sector
        self.J=J #Ising hopping
        self.γ=γ #anisotropy factor
        self.Γ=Γ #Transverse field
    def paramstr(self):
        #returns a string that contains the parameters of the Hamiltonian
        return'L_'+str(self.N)+',S_'+str(self.S)+',J_'+str(self.J)+',Γ_'+str(self.Γ)+',γ_'+str(self.γ)


#function definitions
def LMG_matrixelement(X:Ham_params,M:float,Mprime:float):
    #computes the matrix element <S,M|H|S,M'>
    value=0 
    if abs(M-Mprime)<10**-5:
        value= (X.J/2)*(1+X.γ*(1-2*X.S*(X.S+1)/X.N))-(M**2)*X.J*(2-X.γ)/X.N
    elif abs(M-(Mprime-2))<10**-5:
        value= X.J*X.γ/(2*X.N)*np.sqrt((X.S*(X.S+1)-(M+2)*(M+1))*(X.S*(X.S+1)-M*(M+1)))
    elif abs(M-(Mprime+2))<10**-5:
        value= X.J*X.γ/(2*X.N)*np.sqrt((X.S*(X.S+1)-(M-2)*(M-1))*(X.S*(X.S+1)-M*(M-1)))
    elif abs(M-(Mprime-1))<10**-5:
        value=-X.Γ*np.sqrt(X.S*(X.S+1)-M*(M+1))
    elif abs(M-(Mprime+1))<10**-5:
        value=-X.Γ*np.sqrt(X.S*(X.S+1)-M*(M-1))
    return value         
def LMG_generateHam(X:Ham_params):
    #Generate (2*S+1,2*S+1) matrix.
    Ham=np.zeros((2*X.N+1,2*X.N+1))
    Marr=np.linspace(-X.N//2,X.N//2,2*X.N+1)
    for p in range(np.size(Marr)):
        for q in range(np.size(Marr)):
            Ham[p,q]=LMG_matrixelement(X,Marr[p],Marr[q])
    return Ham
def magnetization2(state,N):
    #takes in a column vector denoting wavefunction, and calcuates average magnetization as defined for ising spins (See lyx file for def)
    Marr=np.linspace(-N//2,N//2,2*N+1)
    magsq=4/N**2*np.sum(np.square(np.abs(state)*Marr))
    return magsq
def Sz2(state,N):
    #takes in a column vector denoting wavefunction, and calcuates average Sz squared.
    Marr=np.linspace(-N//2,N//2,2*N+1)
    Szsq=np.sum(np.square(np.abs(state)*Marr))
    return Szsq 

def time_evolved_Sz2(InitState,Nsteps,U_dt,N):
    #returns an array with the Sz^2 calculated at intervals of dt for Nsteps
    Sz2arr=np.zeros(Nsteps)
    ψ_t=np.copy(InitState)
    for p in np.arange(Nsteps):
        print(p) #print(p, end='\r', flush=True)
        ψ_t=np.dot(U_dt,ψ_t)
        Sz2arr[p]=Sz2(ψ_t,N)
    return Sz2arr

def save_data_Sz2t(paramvals0:Ham_params,paramvalsf:Ham_params,Sz2arr,initstate,Nsteps,dt):
    # saves data in a h5py dictionary
    directory='data/Sz2t/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename=directory+'Sz2t_[0_'+str(dt)+'_'+str(dt*Nsteps)+']_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'
    print(filename)
    with h5py.File(filename, "w") as f:
        f.create_dataset("Sz2arr", Sz2arr.shape, dtype=Sz2arr.dtype, data=Sz2arr)
        f.create_dataset("InitState", initstate.shape, dtype=initstate.dtype, data=initstate)
        f.close()
    with open("list_of_Sz2t.txt", "a") as myfile:
        myfile.write(filename)
    

