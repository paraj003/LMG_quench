directory1='data/FiniteTempTwotimecorrelation/'

filename1=directory1+'FiniteTempTwotimecorrelator_β_'+str(β)+'_Az_'+str(float(Az))+'_Ay_'+str(float(Ay))+'_t1_'+LMG.arrtostr(t1arr)+'_t2_'+LMG.arrtostr(t2arr)+'_from_'+paramvals.paramstr()+'.hdf5'    

if os.path.exists(filename1):
    print("Simulation unnecessary. File:\n"+filename1)
else:
    #twotime correlator    
    if not os.path.exists(filename1):
        finitetempcorr=LMG.finitetemp_twotimecorrelation(paramvals,t1arr,t2arr,β,Az,Ay)
        LMG.save_data_finitetemp_twotimecorrelation(β,paramvals,finitetempcorr,t1arr,t2arr,Az,Ay) 
