directory1='data/Sϕ2t/'

filename1=directory1+'Twotimecorrelator_Az_'+str(1.)+'_Ay_'+str(0.)+'_t1_['+str(t1arr[0])+'_'+str((t1arr[-1]-t1arr[0])/np.size(t1arr))+'_'+str(t1arr[-1])+']_t2_['+str(t2arr[0])+'_'+str((t2arr[-1]-t2arr[0])/np.size(t2arr))+'_'+str(t2arr[-1])+']_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'
filename2=directory1+'Twotimecorrelator_Az_'+str(0.)+'_Ay_'+str(1.)+'_t1_['+str(t1arr[0])+'_'+str((t1arr[-1]-t1arr[0])/np.size(t1arr))+'_'+str(t1arr[-1])+']_t2_['+str(t2arr[0])+'_'+str((t2arr[-1]-t2arr[0])/np.size(t2arr))+'_'+str(t2arr[-1])+']_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'
if not os.path.exists(filename1) or not os.path.exists(filename2):
    #Initial state Ground state of the paramagnetic phase
    Ham0=LMG.LMG_generateHam(paramvals0)
    #(Ham0.transpose() == Ham0).all() #check hermitian
    GSenergy,vec=spla.eigs(Ham0,k=1,which="SR",tol=10**(-6))
    InitState=vec[:,0]
    print("Obtained Initial state.")
    #quench hamiltonian
    Hamf=LMG.LMG_generateHam(paramvalsf)
    energyf,vecf=LA.eig(Hamf)
    #twotime correlator
    U_t1=LMG.LMG_Ut(1.2,energyf,vecf)
    testmat=LMG.Sϕ_on_state(paramvals0,np.dot(U_t1,vec),1,0)
    #correlationmatz=LMG.twotimecorrelation(paramvals0,t1arr,t2arr,vec,energyf,vecf,1,0)
    #LMG.save_data_twotimecorrelation(paramvals0,paramvalsf,correlationmatz,t1arr,t2arr,1,0)
    #correlationmaty=LMG.twotimecorrelation(paramvals0,t1arr,t2arr,vec,energyf,vecf,0,1)
    #LMG.save_data_twotimecorrelation(paramvals0,paramvalsf,correlationmaty,t1arr,t2arr,0,1)
else:
    print("Simulation unnecessary. File:\n"+filename1+"\n"+filename2)

