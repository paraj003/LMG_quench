directory1='data/Sz2t/'

filename=directory1+'Sz2t_[0_'+str(dt)+'_'+str(dt*Nsteps)+']_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'
if not os.path.exists(filename):
    #Initial state Ground state of the paramagnetic phase
    Ham0=LMG.LMG_generateHam(paramvals0)
    (Ham0.transpose() == Ham0).all() #check hermitian
    GSenergy,vec=spla.eigs(Ham0,k=1,which="SR",tol=10**(-8))
    InitState=vec[:,0]
    print("Obtained Initial state.")
    #quench hamiltonian
    Hamf=LMG.LMG_generateHam(paramvalsf)
    U_dt=LA.expm(-1j*Hamf*dt)
    #time-evolved magnetization squared
    Sz2arr=np.zeros(Nsteps)
    Sz2arr=LMG.time_evolved_Sz2(InitState,Nsteps,U_dt,L)
    LMG.save_data_Sz2t(paramvals0,paramvalsf,Sz2arr,InitState,Nsteps,dt)
else:
    print("Simulation unnecessary. File:"+filename)

