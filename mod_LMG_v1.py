# Module containing all definitions necessary to run a quench for the LMG model. Use python 3.
#v1.0 with new defined Hamiltonian with parameters γz and γy
import numpy as np
import os
import h5py



# Class definition to define Hamiltonian

#define Hamiltonian parameters
class Ham_params:
    def __init__(self, N:int,S:float,J:float,γz:float,γy:float,Γ:float):
        self.N=int(N) #number of spins, keep it even
        if S<=N/2:
            self.S=float(S) #spin sector
        else:
            raise Exception("Total spin S must be smaller than N/2")
        self.J=float(J) #Ising hopping
        self.γz=float(γz) #z-direction interaction  
        self.γy=float(γy) #y direction interaction
        self.Γ=float(Γ) #Transverse field
    def paramstr(self):
        #returns a string that contains the parameters of the Hamiltonian
        return'L_'+str(int(self.N))+',S_'+str(float(self.S))+',J_'+str(float(self.J))+',Γ_'+str(float(self.Γ))+',γz_'+str(float(self.γz))+',γy_'+str(float(self.γy))


#function definitions
def LMG_matrixelement(X:Ham_params,M:float,Mprime:float):
    #computes the matrix element <S,M|H|S,M'>
    value=0 
    if abs(M-Mprime)<10**-5:
        value= (X.J/2)*(X.γz+X.γy*(1-2*X.S*(X.S+1)/X.N))-(M**2)*X.J*(2*X.γz-X.γy)/X.N
    elif abs(M-(Mprime-2))<10**-5:
        value= X.J*X.γy/(2*X.N)*np.sqrt((X.S*(X.S+1)-(M+2)*(M+1))*(X.S*(X.S+1)-M*(M+1)))
    elif abs(M-(Mprime+2))<10**-5:
        value= X.J*X.γy/(2*X.N)*np.sqrt((X.S*(X.S+1)-(M-2)*(M-1))*(X.S*(X.S+1)-M*(M-1)))
    elif abs(M-(Mprime-1))<10**-5:
        value=-X.Γ*np.sqrt(X.S*(X.S+1)-M*(M+1))
    elif abs(M-(Mprime+1))<10**-5:
        value=-X.Γ*np.sqrt(X.S*(X.S+1)-M*(M-1))
    return value         
def LMG_generateHam(X:Ham_params):
    #Generate (2*S+1,2*S+1) matrix.
    Ham=np.zeros((int(2*X.S+1),int(2*X.S+1)))
    Marr=np.linspace(-X.S,X.S,int(2*X.S+1))
    for p in range(np.size(Marr)):
        for q in range(np.size(Marr)):
            Ham[p,q]=LMG_matrixelement(X,Marr[p],Marr[q])
    return Ham
def magnetizationz2(state,X:Ham_params):
    #takes in a column vector denoting wavefunction, and calcuates average magnetization along z direction as defined for ising spins (See lyx file for def)
    Marr=np.linspace(-X.S,X.S,int(2*X.S+1))
    magsq=4/X.N**2*np.sum(np.square(np.abs(state)*Marr))
    return magsq
def magnetizationϕ2(state,X:Ham_params,Az:complex,Ay:complex):
    #calculates magnetization squared along 'ϕ'-direction: Sϕ=(Az*Sz+Ay*Sy) (see lyx file for def)
    Marr=np.linspace(-X.S,X.S,int(2*X.S+1))
    A_Marr=np.zeros(np.size(Marr),dtype=complex)
    A_Marr[0]=(Az*X.S*state[0]-np.sqrt(2*X.S)*state[1]*Ay/(2*1j))
    A_Marr[-1]=(Az*X.S*state[-1]+np.sqrt(2*X.S)*state[-2]*Ay/(2*1j))
    for m in range(1,np.size(Marr)-1):
        A_Marr[m]=(Az*Marr[m]*state[m]+Ay/(2*1j)*(np.sqrt(X.S*(X.S+1)-Marr[m]*(Marr[m]-1))*state[m-1]-np.sqrt(X.S*(X.S+1)-Marr[m]*(Marr[m]+1))*state[m+1]))
    magsq=4/X.N**2*np.sum(np.square(np.abs(A_Marr)))
    return magsq 

def Sz2(state,X:Ham_params):
    #takes in a column vector denoting wavefunction, and calcuates average Sz squared.
    Marr=np.linspace(-X.S,X.S,(2*X.S+1))
    Szsq=np.sum(np.square(np.abs(state)*Marr))
    return Szsq 

def Sϕ2(state,X:Ham_params,Az:complex,Ay:complex):
    #takes in a column vector denoting wavefunction, and calcuates average <Sϕdag*Sϕ> Sϕ=Az*Sz+Ay*Sy.
    Marr=np.linspace(-X.S,X.S,(2*X.S+1))
    Sϕsq=X.N**2/4*magnetizationϕ2(state,X,Az,Ay)
    return Sϕsq 

def time_evolved_Sz2(InitState,Nsteps,U_dt,X:Ham_params):
    #returns an array with the Sz^2 calculated at intervals of dt for Nsteps
    Sz2arr=np.zeros(Nsteps)
    ψ_t=np.copy(InitState)
    for p in np.arange(Nsteps):
        #print(p) #print(p, end='\r', flush=True)
        ψ_t=np.dot(U_dt,ψ_t)
        Sz2arr[p]=Sz2(ψ_t,X)
    return Sz2arr

def time_evolved_Sϕ2(InitState,Nsteps,U_dt,X:Ham_params,Az:complex,Ay:complex):
    #returns an array with the Sz^2 calculated at intervals of dt for Nsteps
    Sϕ2arr=np.zeros(Nsteps)
    ψ_t=np.copy(InitState)
    for p in np.arange(Nsteps):
        #print(p) #print(p, end='\r', flush=True)
        ψ_t=np.dot(U_dt,ψ_t)
        Sϕ2arr[p]=Sϕ2(ψ_t,X,Az,Ay)
    return Sϕ2arr

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
 
def save_data_Sϕ2t(paramvals0:Ham_params,paramvalsf:Ham_params,Sϕ2arr,Az,Ay,initstate,Nsteps,dt):
    # saves data in a h5py dictionary
    directory='data/Sϕ2t/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    filename=directory+'Sϕ2t_Az_'+str(float(Az))+'_Ay_'+str(float(Ay))+'_[0_'+str(dt)+'_'+str(dt*Nsteps)+']_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'
    print(filename)
    with h5py.File(filename, "w") as f:
        f.create_dataset("Sϕ2arr", Sϕ2arr.shape, dtype=Sϕ2arr.dtype, data=Sϕ2arr)
        f.create_dataset("InitState", initstate.shape, dtype=initstate.dtype, data=initstate)
        f.close()
    with open("list_of_Sϕ2t.txt", "a") as myfile:
        myfile.write(filename)
    
#def Sz_2_t_tprime(InitState,Nsteps,U_dt,N):
