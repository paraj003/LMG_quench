directory1='data/Twotimecorrelation/'

filename1=directory1+'Twotimecorrelator_Az_'+str(1.)+'_Ay_'+str(0.)+'_t1_'+LMG.arrtostr(t1arr)+'_t2_'+LMG.arrtostr(t2arr)+'_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'
filename2=directory1+'Twotimecorrelator_Az_'+str(0.)+'_Ay_'+str(1.)+'_t1_'+LMG.arrtostr(t1arr)+'_t2_'+LMG.arrtostr(t2arr)+'_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'

if os.path.exists(filename1) and os.path.exists(filename2):
    print("Simulation unnecessary. File:\n"+filename1+"\n"+filename2)
else:
    #Initial state Ground state of the paramagnetic phase
    Ham0=LMG.LMG_generateHam(paramvals0)
    #(Ham0.transpose() == Ham0).all() #check hermitian
    GSenergy,vec=spla.eigs(Ham0,k=1,which="SR",maxiter=200000)
    InitState=vec[:,0]
    print("Obtained Initial state.")
    #quench hamiltonian
    Hamf=LMG.LMG_generateHam(paramvalsf)
    energyf,vecf=LA.eig(Hamf)
    
    #twotime correlator    
    if not os.path.exists(filename1):
        correlationmatz=LMG.twotimecorrelation(paramvals0,t1arr,t2arr,InitState,energyf,vecf,1,0)
        LMG.save_data_twotimecorrelation(paramvals0,paramvalsf,correlationmatz,t1arr,t2arr,1,0)
        print("Sz correlation done")
    if not os.path.exists(filename2):
        correlationmaty=LMG.twotimecorrelation(paramvals0,t1arr,t2arr,InitState,energyf,vecf,0,1)
        LMG.save_data_twotimecorrelation(paramvals0,paramvalsf,correlationmaty,t1arr,t2arr,0,1)
        print("Sy correlation done")

