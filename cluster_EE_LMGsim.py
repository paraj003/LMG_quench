directory='data/CGmats/'
directory1='data/EE/'
filename=directory1+'EE_LA_['+str(part_szarr[0])+'_'+str(part_szarr[-1])+']_t_'+LMG.arrtostr(tarr)+'_from_'+paramvals0.paramstr()+'_to_'+paramvalsf.paramstr()+'.hdf5'
         
entropyarr=np.zeros((np.size(tarr),np.size(part_szarr)))
if not os.path.exists(filename):
    print("Running...")
    Ham=LMG.LMG_generateHam(paramvals0)
    GSenergy,vec=spla.eigs(Ham,k=1,which="SR")
    GState=vec[:,0]
    Hamquench=LMG.LMG_generateHam(paramvalsf)
    energyf,vecf=LA.eig(Hamquench)
    for t,q in zip(tarr,range(np.size(tarr))):
        print(q)
        U_t=LMG.LMG_Ut(t,energyf,vecf)
        ψ_t=np.dot(U_t,GState)
        for p,part_sz in zip(range(np.size(part_szarr)),part_szarr):
            #print([p,q], end='\r', flush=True)
            SA=part_sz/2 #Size of A subsystem
            SB=Stot-SA #Size of B subsystem
            cgmat=LMG.CGmatrix(SA,SB,Stot,directory)
            ψ_tAB=np.matmul(cgmat,ψ_t)
            ρA=LMG.Reduced_ρ(ψ_tAB,SA,SB)
            entropyarr[q,p]=LMG.EEntropy_VN(ρA)
    LMG.save_data_EE(paramvals0,paramvalsf,entropyarr,tarr,GState,part_szarr)
else:
    print("Simulation unnecessary. File:\n"+filename)

