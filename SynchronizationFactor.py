import numpy as np

def SynchronizationFactor(mutant, WTinput,Nmax=10000):
    mutant  = np.asarray(mutant).reshape(-1,).astype(np.double)
    WTinput = np.asarray(WTinput).reshape(-1,).astype(np.double)

    mutant = mutant/np.sum(mutant)

    #computation of N : (only two decimals after the coma are considered for WTinput)
    #N is the number of tested windows
    diff1 = np.round(100*WTinput)-100*WTinput
    diff2 = np.round(10*WTinput)-10*WTinput
    diff3 = np.round(WTinput)-WTinput
    if np.sum(np.abs(diff3))==0:
        k=1
    elif np.sum(np.abs(diff2))==0:
        k=10
    else :
        k=100
    
    WT = np.round(k*WTinput) 
    s = np.sum(WT)
    N = s
    toto = np.floor(Nmax/N)
    if toto>0:
        N = toto * N
    DWT = (np.round(WT/s * N)).astype(np.int64)
    #for information for the user.
    DWTDisplay = s/k/N*DWT.astype(np.double)
    
    difference = np.sum(np.abs(DWTDisplay-WTinput))
    if difference>0.0001:
        print(f"WT {WTinput} is approximated by {DWTDisplay}.")
    if np.min(DWT)==0:
        print("We cannot run this algorithm when  a stage is 0 for the WT")
        print("However, the algorithm can be modified easily. I suggest to you to write to the developer")
        return 0;
    Ntot  =  np.sum(DWT)
 
    
    cycle = np.zeros((Ntot,))
    debut = 0
    for i in range(mutant.shape[0]):
        cycle[debut:debut+DWT[i]] =  mutant[i]/DWT[i] * np.ones((int(DWT[i]),))
        debut = debut + DWT[i] 

    res=0
    for i in range(cycle.shape[0]):
        cumcycle= np.cumsum(cycle)
        indice = np.where(cumcycle >= 0.6826)[0][0]
        tmp = np.double(cycle.shape[0]) * cumcycle[indice] / np.double(indice+1)
        #cumcycle[indice] should be close to 0.6826 is N is large.
        #Note that we use the same percentiles for the WT and for the mutant so that the result is not biased.
        if tmp>res:
            res=tmp
        cycle = np.roll(cycle,-1)

    
    print("The synchronization factor is {:.2f}".format(res))
    return res
